# Created by: Erica Behrens
# Send any questions to eabehrens97@gmail.com
# Revised 12/23/21

# Script to extract molecualr line integrated intensity image pixel values averaged over a user-defined area for entire image
# See 00README.txt file for instructions on how to run code





from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.wcs import WCS
import regions, os, sys, math
from astropy.table import Table
from astropy.io import ascii
import os
import pandas as pd


### CHANGE OBJECT DISTANCE ###
# distance in pc
dist = 3.34e6 # NGC 253 distance in pc


############################## DEFINED FUNCTIONS #############################

# MAKE HEXAGON REGIONS AND SAVE TO .REG FILE 

def build_hex_reg(h,l,csvrow,makeCirc):

    impath = csvrow['intimage']
    hdr = fits.getheader(impath)
    data = fits.getdata(impath)
    w = WCS(hdr)
    pix_val = hdr['CDELT2']  #pixel size in deg
    xdim = hdr['NAXIS1']
    ydim = hdr['NAXIS2']
    num_cols = math.floor((xdim*pix_val - l/2) / (3*l/2))
    num_rows = math.floor(ydim * pix_val / h)


    # make hexagon .reg file--cover image
    # first hexagon at origin (bottom left corner)

    p0 = w.wcs_pix2world(0,0,0)      # origin
    p1 = [p0[0] - l/2, p0[1]]        # bottom left corner
    p2 = [p1[0] + l/2, p1[1] + h/2]  # middle left corner
    p3 = [p2[0] - l/2, p2[1] + h/2]  # top left corner
    p4 = [p3[0] - l, p3[1]]          # top right corner
    p5 = [p4[0] - l/2, p4[1] - h/2]  # middle right corner
    p6 = [p5[0] + l/2, p5[1] - h/2]  # bottom right corner

    ptop = w.wcs_pix2world(0,ydim,0)    # coordinates of top row of image
    pright = w.wcs_pix2world(xdim,0,0)  # coordinates of right-most column of image

    h_pc = h*u.deg.to(u.radian) * dist
    if makeCirc:
        shape = 'circ'
    else:
        shape = 'hex'
    f = open('output/' + csvrow['species']+csvrow['transition']+ '_' + shape +  str(round(h_pc)) + 'pc' +'_full.reg','w')
    f.write('fk5\n')

    col = 0

    # Starting from bottom left corner, build hexagons in a column upwards until you hit the top row. If another full hexagon can't fit onto the image, start the next column at the bottom of the image again and work your way toward the right until another column of hexagons cannot fit on the image.
   
    while abs(pright[0] - p6[0]) > 2*l:   


        for row in range(num_rows):

            if row == 0 and col != 0:
                if col % 2 == 0:  # check if column is even or odd--even columns start at y = 0, odd start at y = h/2
                    p1 = [p0[0] - l * (0.5 + 3./2 * col), p0[1]]
                    p2 = [p1[0] + l/2, p1[1] + h/2]
                    p3 = [p2[0] - l/2, p2[1] + h/2]
                    p4 = [p3[0] - l, p3[1]]
                    p5 = [p4[0] - l/2, p4[1] - h/2]
                    p6 = [p5[0] + l/2, p5[1] - h/2]

                else:

                    p1 = [p0[0] - 2*l - 3*l * (col // 2), p0[1] + h/2]
                    p2 = [p1[0] + l/2, p1[1] + h/2]
                    p3 = [p2[0] - l/2, p2[1] + h/2]
                    p4 = [p3[0] - l, p3[1]]
                    p5 = [p4[0] - l/2, p4[1] - h/2]
                    p6 = [p5[0] + l/2, p5[1] - h/2]

            elif row == (num_rows - 1):  # make sure there's enough room at the top for another hex
                if abs(ptop[1] -  p3[1]) < h:
                    continue

                p1 = [p1[0], p1[1] + h]
                p2 = [p2[0], p2[1] + h]
                p3 = [p3[0], p3[1] + h]
                p4 = [p4[0], p4[1] + h]
                p5 = [p5[0], p5[1] + h]
                p6 = [p6[0], p6[1] + h]

            elif row != 0:
                p1 = [p1[0], p1[1] + h]
                p2 = [p2[0], p2[1] + h]
                p3 = [p3[0], p3[1] + h]
                p4 = [p4[0], p4[1] + h]
                p5 = [p5[0], p5[1] + h]
                p6 = [p6[0], p6[1] + h]


            # write region to file
            if not makeCirc:
                f.write('polygon(' + str(p1[0]) + ',' + str(p1[1]) + ',' + str(p2[0]) +  ',' + str(p2[1]) + ',' + str(p3[0]) + ',' + str(p3[1]) + ',' + str(p4[0]) + ',' + str(p4[1]) +  ',' + str(p5[0]) + ',' + str(p5[1]) + ',' + str(p6[0]) + ',' + str(p6[1]) + ') # color=green\n')
                
            else:
                cra = np.mean([p1[0],p6[0]])
                cdec = p2[1]
                f.write('circle(' + str(cra) + ',' + str(cdec) + ',' + str(h/2) + ') # color=green\n') 
                

        col += 1

    f.close()


        

# CALCULATE REGION STATISTICS 

def calc_stats(outtab,regfile,switch,makeCirc,nanpct,csvrow):

    impath = csvrow['intimage']
    hdr = fits.getheader(impath)
    data = fits.getdata(impath)
    w = WCS(hdr)
    nan_mask = np.ma.masked_array(data,np.isnan(data))  # mask image NaN vals
    reglist = []

    #i = 0
    for reg in regfile:
        pixreg = reg.to_pixel(w)
        reg_mask = pixreg.to_mask()
        reg_arr = reg_mask.cutout(data)    # select region pixels

        # check if percentage of nan vals is below cutoff
        if not switch :
            if (len(reg_arr[np.isnan(reg_arr)]) / reg_arr.size) > nanpct:
                continue
            else:  # no nanpct variable when user enters their own .reg file
                reglist.append(reg)
                if makeCirc:
                    ra_cen = reg.center.ra.value
                    dec_cen = reg.center.dec.value
                else:
                    ra_cen = np.mean([reg.vertices[1].ra.value,reg.vertices[4].ra.value])
		    #ra_cen.append(np.mean([reg.vertices[1].ra.value,reg.vertices[4].ra.value]))
		    #dec_cen.append(reg.vertices[1].dec.value)
                    dec_cen = reg.vertices[1].dec.value
        else:
            if reg_type == 'circ':
                ra_cen = np.float64(reg.center.ra.value)
                dec_cen = np.float64(reg.center.dec.value)
            elif reg_type == 'hex':
                ra_cen = np.mean([reg.vertices[1].ra.value,reg.vertices[4].ra.value])
                dec_cen = reg.vertices[1].dec.value
            reglist.append(reg)


        # if region falls off edge of image:
        if nan_mask[reg_mask.bbox.slices].shape != reg_mask.shape:
            continue 
        
        reg_mean = np.ma.average(nan_mask[reg_mask.bbox.slices],weights=reg_mask)
        var = np.ma.average((nan_mask[reg_mask.bbox.slices]-reg_mean)**2,weights=reg_mask)
        rmsaver = math.sqrt(var)
        #rmsint = csvrow['chanrms']*0.001*float(chanwidth)*math.sqrt(csvrow['fwzi']/float(chanwidth))

        rmsint = np.sqrt((0.001*csvrow['chanrms']*float(csvrow['chanwidth']))**2*(csvrow['fwzi']/float(csvrow['chanwidth'])) + (reg_mean*float(csvrow['fluxcaluncert']))**2)
        rmsint_measnoiseportion = (0.001*csvrow['chanrms']*float(csvrow['chanwidth']))**2*(csvrow['fwzi']/float(csvrow['chanwidth']))
        rmsint_fluxcalportion = (reg_mean*float(csvrow['fluxcaluncert']))**2
        #i += 1

        outtab.add_row([csvrow['species'], csvrow['transition'], csvrow['intimage'], csvrow['chanrms'], csvrow['chanwidth'],csvrow['fluxcaluncert'],ra_cen, dec_cen, reg_mean, rmsaver, rmsint,rmsint_measnoiseportion,rmsint_fluxcalportion])

    return outtab,reglist



###################################################################

########################     MAIN     #############################

###################################################################

               
# Allow running this thing from anywhere...

path = os.getcwd()
sys.path.append(path)
files = os.listdir(path)

# Set floating-point error handling for arithmetic

np.seterr(all='warn')



while os.path.isfile:
    csvinput = input('Enter csv file name (see README for input format): ')
    if os.path.isfile(csvinput) == True:
        csvtable = pd.read_csv(csvinput)
        break
    else:
        print(csvinput,' does not exist...try again')




y_n = input('\nWould you like to enter your own ds9 .reg file (only circular and hexagonal regions are accepted)? ')

################# IF USER WANTS TO ENTER THEIR OWN .REG FILE ###################

if y_n == 'yes' or y_n == 'y':
    
    switch = True
    makeCirc = False
    nan_pct = 0

    while os.path.isfile:
        reg_path = input('Enter .reg file path: ')
        if os.path.isfile(reg_path) == True:
            regfile = regions.read_ds9(reg_path)
            break
        else:
            print(reg_path,' does not exist...try again')

    if type(regfile[0]) == regions.shapes.circle.CircleSkyRegion:
        reg_type = 'circ'
    elif type(regfile[0]) == regions.shapes.polygon.PolygonSkyRegion:
        reg_type = 'hex'
        
    out_tab = Table(names=('species','transition','intimage','chanrms','chanwidth','fluxcaluncert','regionra','regiondec','averint','rmsaver','rmsint','rmsint_measnoiseportion','rmsint_fluxcalportion'), dtype=('U32','U32','U128','float64','float64','float64','float64','float64','float64','float64','float64','float64','float64'))

    for i in range(len(csvtable)):

        csvrow = csvtable.iloc[i]

        print('\nWorking on ' + csvrow['species'] + ' ' + csvrow['transition'])

        impath = csvrow['intimage']
        out_tab,reg_list = calc_stats(out_tab,regfile,switch,makeCirc,nan_pct,csvrow)
        out_tab.write('region_avg_stats.csv',format='csv',overwrite=True)




############# IF USER WANTS CODE-GENERATED REGIONS #################
        
else:
                             
    switch = False

    # ask if user wants hexagonal or circular regions
    shape = input('\nWhat shape regions would you like to generate (enter hex or circle): ')
    if shape != 'hex' and shape != 'hexagon' and shape != 'h' and shape != 'circle' and shape != 'circ' and shape != 'c':
        shape = input('Shape not recognized. Enter hex or circle: ')
        
    unit = input('\nWhat unit are you using to specify region size (enter pc or arcsec): ')
    if unit != 'pc' and unit != 'arcsec':
            unit = input('Unit not recognized. Please enter pc or arcsec: ')
            
    if shape == 'hex' or shape  == 'hexagon' or shape == 'h':
        makeCirc = False
        shape = 'hex'
        s = float(input('\nEnter desired hexagon size (length from top to bottom): '))
    else:
        makeCirc = True
        shape = 'circ'
        s = float(input('\nEnter desired circle diameter: '))
        
        
    if unit == 'pc':
        s = (s / dist)*u.radian.to(u.deg)   # change reg size to degrees
    else:                           
        s = s*u.arcsec.to(u.deg)            # change reg size to degrees
    l = s / np.sqrt(3)                      # length of hexagon side
    s_pc = s*u.deg.to(u.radian) * dist


    nan_pct = float(input('\nSome regions will include masked pixels with nan values. Enter the maximum fraction of a region that can be filled with nan values: '))
    if nan_pct > 1:
        nan_pct = input('Input must be a fraction less than one (e.g. 0.25): ')
    out_tab = Table(names=('species','transition','intimage','chanrms','chanwidth','fluxcaluncert','regionra','regiondec','averint','rmsaver','rmsint','rmsint_measnoiseportion','rmsint_fluxcalportion'), dtype=('U32','U32','U128','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8'))

    if not os.path.isdir('output'):
        os.mkdir('output')

    for i in range(len(csvtable)):

        csvrow = csvtable.iloc[i]

        print('\nWorking on ' + csvrow['species'] + ' ' + csvrow['transition'])

        impath = csvrow['intimage']
        build_hex_reg(s,l,csvrow,makeCirc)

        regfile = regions.read_ds9('output/' + csvrow['species']+csvrow['transition']+ '_' + shape + str(round(s_pc)) + 'pc' + '_full.reg')
        out_tab,reg_list = calc_stats(out_tab,regfile,switch,makeCirc,nan_pct,csvrow)

        regions.write_ds9(reg_list,csvrow['species']+csvrow['transition']+'_' + shape + str(round(s_pc)) + 'pc_' + str(nan_pct)+ '_final.reg')
        out_tab.write('reg_avg_'+ shape + str(round(s_pc)) + 'pc_' + str(nan_pct) + '.csv', format='csv', overwrite=True)
