README file for region_average.py
Direct any questions to Erica Behrens at eabehrens97@gmail.com



-----------------   BEFORE RUNNING CODE --------------------


Before running region_average.py, you will need to create a CSV file with one row per molecular transition and with the following columns:
Species, Transition, Single-channel RMS (Jy), Path to integrated intensity images, Full Width at Zero Intensity (km/s), Spectral channel width (km/s), Flux calibration uncertainty (fraction from 0 to 1)

The column names that are referenced in the code and that should appear in the CSV file are:
'species','transition','chanrms','intimage','fwzi', 'chanwidth','fluxcaluncert'

An example input CSV file row:

HC15N	1-0	0.18	/Users/eb7he/research_uva/HCNHNC/CubeLineMoment/moment0/NGC253_HC15N_10_moment0_widthscale1.0_sncut3.0_widthcutscale1.0.fits	100	10.0	0.15



------------------  USER INPUT ------------------------

Run the code, and when prompted, enter the full path to the CSV file (unless it is in your working directory, then you can just use the file name).

You will be asked whether you would like to input your own ds9 .reg file. If you would like to specify your own regions to average over, enter 'yes' and then enter the path to your ds9 .reg file.

If you would like the code to define regions for you, enter 'no' and then follow the prompts to choose your desired region shape (hexagon or circle) and size. Note that using hexagonal regions will sample the entire image, while circular regions result in some wasted space between regions. If using code-defined regions, you will be asked to enter a maximum percentage of the region that can be occupied by Not a Number (NaN) values (if images are masked). Any regions with a higher percentage of NaN values than the user-defined fraction will be discarded when doing statistics and averaging. Enter a fraction less than 1.


------------------- CODE OUTPUT -------------------------

If you input your own ds9 .reg file, the only output will be the output CSV file 'region_avg_stats.csv'. The output CSV table will have the following column names:
'species','transition','intimage','chanrms','chanwidth','fluxcaluncert','regionra','regiondec','averint','rmsaver','rmsint','rmsint_measnoiseportion','rmsint_fluxcalportion'
with one row per region.

The first 6 columns are the same as the input csv file.
'regionra' and 'regiondec' = RA and Dec of center pixel in that region 
'averint' = averaged integrated intensity over region
'rmsaver' = spatial rms over region 
'rmsint' = rms of integrated intensity
'rmsint_measnoiseportion' = measured noise portion of the rms, used to calculate 'rmsint'
'rmsint_fluxcalportion' = flux calibration portion of the rms, used to calculate 'rmsint'

If you had the code define regions, you will get the same output CSV file but named "reg_avg_{Shape}{Shapesize}_{Nanpct}.csv", with the user-defined region, size, and NaN percentage in the file name (e.g. 'reg_avg_hex50pc_0.8.csv'). The code will also output 2 ds9 .reg files per transition: one with full region coverage over the whole image (before Nan-filled regions are discarded) and one with region coverage only over the key areas in the image (after NaN-filled/background regions are discarded), named respectively "{Species}{Transition}_{Shape}{Size}_full.reg" and "{Species}{Transition}_{Shape}{Size}_{Nanpct}_final.reg" (e.g. 'HNC3-2_hex26pc_full.reg' and 'HNC3-2_hex26pc_0.5_final.reg').
