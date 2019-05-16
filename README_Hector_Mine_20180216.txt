#Demonstration processing or Hector Mine interferometric pair cw 16-Feb-2018
#set the default raster data type to BMP
export GAMMA_RASTER='BMP'

#create CEOS_list file with directories of the raw data

cat CEOS_list
19990915
19991020

#ERS raw data preprocessing script
ERS_pre_proc
*** /Users/cw/gamma_software/MSP/scripts/ERS_pre_proc
*** Copyright 2014, Gamma Remote Sensing, v4.0 7-Jun-2014 clw ***
*** ERS raw data pre-processing ***

usage: /Users/cw/gamma_software/MSP/scripts/ERS_pre_proc <CEOS_list> <path_to_DELFT> <out_dir> <log> <proc_list> <mode> [keyword] [value]

    CEOS_list	(input) parameter file with 1 column entry/record
		  1. directory (including path) containing CEOS data set from ERS
    DELFT_path  path to DELFT orbits (example: /usr/local/data/DELFT)
    out_dir     directory for output fixed raw data files and processing parameter files
    log         (output) ERS_pre_proc processing log file
    proc_list   (output) processing list for use by ERS_proc_par and ERS_proc_all
    mode        processing mode:
		  1: create processing parameter files from CEOS leaders
		  2: fix raw data files (read and scan for missing lines)
		  3: extract and interpolate DELFT state vectors and update MSP processing parameter files
		  4: estimate Doppler ambiguity and centroid
    		  5: set value in the processing parameter files for a keyword:value pair
		  6: generate processing list for use by ERS_proc_all
		  7: gzip fixed raw data files
    keyword 	(mode 5) keyword in the MSP processing parameter file (example: doppler_polynomial)
    value	(mode 5) new value delimited by double quotes (example: "317.0 0. 0. 0.")

#create processing parameter files from CEOS leaders
ERS_pre_proc CEOS_list ./DELFT raw ERS_pre_proc_1.log proc_list 1

#fix raw data files (read and scan for missing lines)
ERS_pre_proc CEOS_list ./DELFT raw ERS_pre_proc_2.log proc_list 2

#extract and interpolate DELFT state vectors and update MSP processing parameter files
ERS_pre_proc CEOS_list ./DELFT raw ERS_pre_proc_3.log proc_list 3

#estimate Doppler ambiguity and centroid
ERS_pre_proc CEOS_list ./DELFT raw ERS_pre_proc_4.log proc_list 4

#generate processing list file proc_list for use by ERS_proc_all
ERS_pre_proc CEOS_list ./DELFT raw ERS_pre_proc_6.log proc_list 6

cat proc_list
19990915 ERS2_ESA.par - - - -     29.5 0.0 0.8
19991020 ERS2_CCRS.par - - - -     88.5 0.0 0.8

# ERS raw data processing to generate SLC images
$ ERS_proc_all
*** /Users/cw/gamma_software/MSP/scripts/ERS_proc_all
*** Copyright 2015, Gamma Remote Sensing, v3.5 22-Mar-2016 clw ***
*** ERS-1 and ERS-2 SAR image processing from RAW data ***

usage: /Users/cw/gamma_software/MSP/scripts/ERS_proc_all <proc_list> <raw_dir> <rc_dir> <SLC_dir> <MLI_dir> <rlks> <azlks> <SLC_format> [az_patch] [autof_min] 
    proc_list	(input) processing list (9 columns): 
    		    1. scene_identifier (example: 19960816)
		    2. MSP sensor parameter file
		    3. offset in echoes to start of processed data (enter - for default)
		    4. number of echoes to process (enter - for default)
		    5. range offset in samples (enter - for default)
		    6. number of range samples to process (enter - for default)
		    7. Doppler centroid for scene (Hz)
		    8. Doppler slope for scene (Hz/m)
		    9. azimuth processing bandwidth fraction for scene 
    raw_dir	 data directory containing fixed ERS raw data files 
    rc_dir	 directory to temporarily store intermediate range compressed data (example: /tmp)
    SLC_dir	 directory to store output SLC data 
    MLI_dir	 directory to store multilook intensity (MLI) files derived from the slc data
    rlks	 number of range looks to generate MLI images   (nominal: 1)
    azlks	 number of azimuth looks to generate MLI images (nominal: 5)
    SLC_format	 SLC image format:
                    0: FCOMPLEX
                    1: SCOMPLEX
    az_patch     azimuth patch size (enter - for default: 8192)
    autof_min    minimum SNR threshold for autofocus, 0.0 for no autofocus (enter - for nominal = 8)
    
    NOTE:  current directory is denoted by .


#process 2 scenes using autofocus
ERS_proc_all proc_list raw . slc mli_2_10 2 10 0 6144 10
##################################################################################################
#Determine SRTM DEM boundaries:

SLC_corners slc/19990915.slc.par
*** Calculate SLC/MLI image corners in geodetic latitude and longitude (deg.) ***
*** Copyright 2017, Gamma Remote Sensing, v1.8 13-Dec-2017 clw/awi/cm ***
latitude (deg.):      35.09902384  longitude (deg.):    -115.51972142
latitude (deg.):      35.28353472  longitude (deg.):    -116.59450850
latitude (deg.):      34.15576794  longitude (deg.):    -115.76996352
latitude (deg.):      34.33974572  longitude (deg.):    -116.83119021

center latitude (deg.):    34.72746240  center longitude (deg.):  -116.21749899
min. latitude (deg.):      34.15576794  max. latitude (deg.):       35.28353472
min. longitude (deg.):   -116.83119021  max. longitude (deg.):    -115.51972142
delta latitude (deg.):      1.12776678  delta longitude (deg.):      1.31146879

upper left corner latitude,  longitude (deg.): 35.32  -116.87
lower right corner latitude, longitude (deg.): 34.12  -115.48

user time (s):         0.000
system time (s):       0.000
elapsed time (s):      0.020

#Download SRTM 1-arcsec DEM that covers this region
#NASA V3.0 1 arc-sec DEM mosaic from USGS Combined 1 arcs SRTMGL1:
#   http://gdex.cr.usgs.gov/gdex/
#information on the SRTMGL1 data set: 
#   https://lpdaac.usgs.gov/dataset_discovery/measures/measures_products_table/srtmgl1
#check the NASA SRTM V3.0, 1 arcsec
#enter corner coordinates using xy icon on the selection toolbar in the gdex web browser:

upper left corner latitude,  longitude (deg.): 35.32  -116.87
lower right corner latitude, longitude (deg.): 34.12  -115.48

#Download DEM as a single GeoTIFF file, click on the folder download icon on the gdex toolbar
#select NASA SRTM combined images V3.0, 1 arcsec, as lat/lon projection

#The SRTM DEM uses a geoidal reference. The srtm2dem script interpolates a geoid correction map to create a DEM with WGS84 as the height datum:
srtm2dem
*** /Users/cw/gamma_software/DIFF/scripts/srtm2dem
*** Copyright 2017, Gamma Remote Sensing, v1.9 16-Aug-2017 clw/oc/cm ***
*** Convert SRTM DEM in GeoTIFF format to Gamma DEM format with optional geoid offset correction ***

usage: /Users/cw/gamma_software/DIFF/scripts/srtm2dem <SRTM_DEM> <DEM> <DEM_par> <gflg> <geoid> [geoid_dem]
    SRTM_DEM   (input) DEM in geographic coordinates (lat.,lon.) in GeoTIFF format
               NOTE: DEM downloaded from USGS http://gdex.cr.usgs.gov or CGIAR http://srtm.csi.cgiar.org
    DEM        (output) DEM data 
    DEM_par    (output) Gamma DEM parameter file
    gflg       geoid offset correction flag:
                 0: no geoid offset correction, replace NODATA with 1.e-20.
                 1: add constant geoid offset relative to the WGS84 datum, NODATA are set to the offset value
                 2: add interpolated geoid offset relative to the WGS84 datum
                    NODATA are set to the interpolated geoid offset
                 3: no geoid offset correction, replace NODATA with 0.0
    geoid      (input) geoid offset relative to the WGS84 datum  (enter - for default):
                 gflg:1 constant offset (m) (default: 0.0)
                 gflg:2 GeoTIFF file containing geoid offset in the WGS84 datum 
                 (default in the Gamma Software: /Users/cw/gamma_software/DIFF/scripts/egm96_wgs84_diff.tif)
    geoid_dem  (output) resampled geoid offset map in the DEM grid (gflg: 2)

cd DEM
srtm2dem 20150914023253_1455875154.tif hector_eqa.dem hector_eqa.dem_par 2 -

#Terrain geocode the scene from 19990915 using geographic (EQA) coordinates using mk_geo_radcal script
mk_geo_radcal2
*** /Users/cw/gamma_software/DIFF/scripts/mk_geo_radcal2
*** Copyright 2018, Gamma Remote Sensing, v6.7 1-Feb-2018 clw ***
*** Terrain geocoding of SAR images using gc_map2, including lookup-table refinement + resampling of the DEM to SAR Range-Doppler Coordinates (RDC) ***

usage: /Users/cw/gamma_software/DIFF/scripts/mk_geo_radcal2 <MLI> <MLI_par> <DEM> <DEM_par> <DEM_seg> <DEM_seg_par> <GEO_dir> <scene_id> <post> <mode> [msk_mode] [r_ovr] [gc_n_ovr] [rad_max] [rlks] [azlks] [thres] [rpos] [azpos] [roff] [azoff] [r_patch] [az_patch] [nr] [naz]

    MLI          (input) MLI SAR image (including path)
    MLI_par      (input) ISP image parameter file of the MLI image (including path)
    DEM          (input) DEM in desired output projection (including path)
    DEM_par      (input) DEM parameter file (including path)
    DEM_seg      (output) DEM segment for output image products (including path)
    DEM_seg_par  (input/output) DEM parameter file for output image products (including path), regenerated each time 
    GEO_dir      directory for output images, lookup tables and DEM products
    scene_id     scene name to identify output files
    post         output image sample spacing in meters or degrees for geographic (EQA) projection
                 NOTE: for the EQA projection, use the -z option to calculate the longitude posting to give square pixels at the
                 center of the frame, otherwise the latitude and longitude posting are equal
    mode         processing mode:    
                   0: generate initial lookup table, simulated SAR image, and DEM segment parameters, erase existing DEM segment parameters
                   1: measure initial offset between simulated SAR image and actual SAR image
                   2: perform refinement of lookup table by offset measurement with respect to the simulated SAR image
                   3: update lookup table and produce terrain geocoded SAR image and DEM in SAR range-Doppler coordinates (RDC) 
                   4: ellipsoid geocoding of the SAR image without a DEM, erase existing DEM segment parameters  
    msk_mode     gc_map2 masking (enter - for default):
                   0: no masking (default)
                   1: mask values outside the swath
                   2: masking layover, shadow and values outside swath
    r_ovr        gc_map2 range oversampling factor (only used in internal calculations, for improved accuracy)(enter - for default: 2)
    gc_n_ovr     interpolation oversampling factor used by geocode (enter - for default: 2)
    rad_max      maximum interpolation search radius used by geocode (enter - for default: 8*gc_n_ovr)
    
    NOTES: 1. gc_n_ovr and rad_max are parameters used by the program geocode for transformation of the simulated image and DEM into SAR geometry.
           2. The parameters rlks, azlks, thres, rpos, azpos, roff, azoff are used for estimation of the initial offset of the SAR image with
              respect to the simulated SAR image. 
		 		 
    rlks         number of range looks for the initial offset estimate (default: 1)
    azlks        number of azimuth looks for the initial offset estimate (default: 1)
    thres        cross-correlation threshold for offset measurements over offset grid (default: 0.15)
    rpos         range position for initial offset (enter - for default)
    azpos        azimuth position for initial offset (enter - for default)
    roff         initial range offset estimate (enter - for current value in DIFF_par file)
    azoff        initial azimuth offset estimate (enter - for current value in DIFF_par file)
    r_patch      range patch size for offset estimation (default: 512 samples)
    az_patch     azimuth patch size for offset estimation (default: 512 lines)
    nr           number of range patches for offset estimation (default: 0)
    naz          number of azimuth patches for offset estimation (default: 0)

    -s scale     (option) set image display scale factor (default: 0.9)   
    -e exp       (option) set image display exponent (default: 0.4)
    -r           (option) use existing DEM segment to determine image bounds, do not erase an existing DEM segment parameter file
    -q           (option) quiet mode, run without displaying images on screen
    -p           (option) use pixel_area program to generate simulated SAR image in Range-Doppler Coordinates (RDC)
                          NOTE: ls_mode command line parameter must equal actual area: 2
    -j           (option) do not use layover-shadow map in pixel_area calculation
    -i map_pwr   (option) intensity image in map coordinates to be used for lookup table refinement instead of a simulated SAR image,
                          same bounds and sample spacing as the DEM segment
    -c           (option) calculate radiometrically calibrated radar backscatter image
    -d           (option) resample DEM to Range-Doppler Coordinates (RDC)
    -x psize     (option) image patch size for initial offset estimate using cross-correlation in mode 1 (default: 1024)
    -h           (option) reference height for ellipsoid geocoding, mode 4 (default: 100)
    -t thres2    (option) cross-correlation threshold for initial offset estimate (default: 0.1)
    -n npoly     (option) number of terms in the range and azimuth offset polynomials (default: 1)
    -o n_ovr     (option) offset estimation oversampling factor for image data (default: 2)
    -z           (option) calculate longitude posting in degrees to get square pixels for EQA projection at the center of the frame
                          NOTE: by default the latitude and longitude postings are equal and specified by the post command line parameter
  
    NOTE: -s , -e and -t parameters must have a number before the decimal point, e.g. 0.3

#geocode to spacing of 2.0833333e-4 degrees (about 30m) in Latitude and Longitude
#Generate initial lookup table, simulated SAR image, and DEM segment parameters, erase existing DEM segment parameters, use pixel_area to simulate image, scale the longitude posting to get approximately square pixels
mk_geo_radcal2 mli_2_10/19990915.mli mli_2_10/19990915.mli.par DEM/hector_eqa.dem DEM/hector_eqa.dem_par geo/hector_eqa_seg.dem geo/hector_eqa_seg.dem_par geo hector 2.77777778e-4 0 -p -z -d -j

#measure initial offset between simulated SAR image and actual SAR image
mk_geo_radcal2 mli_2_10/19990915.mli mli_2_10/19990915.mli.par DEM/hector_eqa.dem DEM/hector_eqa.dem_par geo/hector_eqa_seg.dem geo/hector_eqa_seg.dem_par geo hector 2.77777778e-4 1 -p -d -j -n 3

#perform refinement of lookup table by offset measurement with respect to the simulated SAR image
mk_geo_radcal2 mli_2_10/19990915.mli mli_2_10/19990915.mli.par DEM/hector_eqa.dem DEM/hector_eqa.dem_par geo/hector_eqa_seg.dem geo/hector_eqa_seg.dem_par geo hector 2.77777778e-4 2 -p -d -j -n 3

range offset polynomial:   -0.82844  -3.4217e-05  -4.2742e-05  0.0000e+00  0.0000e+00  0.0000e+00
azimuth offset polynomial: 0.61020  5.0871e-06  1.3561e-05  0.0000e+00  0.0000e+00  0.0000e+00

#update lookup table and produce terrain geocoded SAR image and DEM in SAR range-Doppler coordinates (RDC), genedarte DEM in radar coordinates (RDC)
mk_geo_radcal2 mli_2_10/19990915.mli mli_2_10/19990915.mli.par DEM/hector_eqa.dem DEM/hector_eqa.dem_par geo/hector_eqa_seg.dem geo/hector_eqa_seg.dem_par geo hector 2.77777778e-4 3 -p -d -j -n 3

dis2ras geo/hector_1.pix.ras mli_2_10/19990915.mli.bmp&
disras_dem_par geo/hector_map.mli.ras geo/hector_eqa_seg.dem_par&
disras geo/hector_dem.rdc.bmp&
##################################################################################################
#coregister SLC images using polynomial model for interfereometry
#generate SLC_tab

mk_tab slc slc slc.par SLC_tab

#output tab file: SLC_tab
#slc/19990915.slc  slc/19990915.slc.par
#slc/19991020.slc  slc/19991020.slc.par

SLC_resamp_all
*** /Users/cw/gamma_software/DIFF/scripts/SLC_resamp_all
*** Copyright 2015, Gamma Remote Sensing, v5.1 18-Aug-2016 clw ***
*** Resample a set of SLCs to a common reference SLC using a polynomial offset model ***

usage: /Users/cw/gamma_software/DIFF/scripts/SLC_resamp_all <SLC_tab> <ref_SLC> <ref_par> <rslc_dir> <RSLC_tab> <mode> [rflag] [rlks] [azlks] [rpos] [azpos] [rpatch] [azpatch] [npoly] [n_ovr]

    SLC_tab   (input) two column list of SLC filenames and SLC parameter filenames (including paths) (ascii)
    ref_SLC   (input) reference SLC (including path)
    ref_par   (input) ISP image parameter file of the reference SLC (including path)
    rslc_dir  directory to receive the resampled SLCs and ISP image parameter files   
    RSLC_tab  (output) RSLC_tab file for the resampled SLC files
    mode      processing mode:
                0: create offset parameter files    
                1: estimate initial range and azimuth offsets using orbit state vectors
                2: measure initial range and azimuth offsets using SLC images
                3: estimate range and azimuth offset models using correlation of image intensities
	        4: resample SLC images using offset models
    rflag     flag for interactive refinement of the resampled SLC:
                0: resample and measure residual range and azimuth offsets  (no interation, default)
		1: interactively improve range and azimuth offset polynomials then measure residual offsets
    rlks      number of range looks for initial offset estimation  (enter - for default)
    azlks     number of azimuth looks for initial offset estimation  (enter - for default)
    rpos      center of patch in range (samples) (enter - for default: image center)
    azpos     center of patch in azimuth (lines) (enter - for default: image center)
    rpatch    range patch size for offset estimation (enter - for default: 128)
    azpatch   azimuth patch size for offset estimation (enter - for default: 256)
    npoly     number of terms in the polynomial fit estimated from offsets (1,3,4,6 default: 6)
    n_ovr     SLC oversampling factor for offset estimation (integer 2**N (1,2,4) default: 2)

    -i "roff azoff"  (option) initial values for range and azimuth offsets, example: "-10 30", default: "- -"
    -t thres         (option) cross-correlation threshold to accept an SLC offset measurement (default: 0.2) 
    -p psz           (option) size of the patch for initial offset estimation (default: 256)
    -c clip          (option) intensity clipping threshold for offset measurement data (default: );
    -u               (option) update flag, coregister only new SLC scenes that have not yet been processed
   
#create offset parameter files, copy reference SLC to rslc directory
SLC_resamp_all SLC_tab slc/19990915.slc slc/19990915.slc.par rslc RSLC_tab 0 

# estimate initial range and azimuth offsets using orbit state vectors
SLC_resamp_all SLC_tab slc/19990915.slc slc/19990915.slc.par rslc RSLC_tab 1 

#estimate initial range and azimuth offsets using SLC images
SLC_resamp_all SLC_tab slc/19990915.slc slc/19990915.slc.par rslc RSLC_tab 2 

#range_offset_polynomial: -172.63754 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00
#azimuth_offset_polynomial: -115.83609 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00

#estimate range and azimuth offset models using correlation of image intensities
SLC_resamp_all SLC_tab slc/19990915.slc slc/19990915.slc.par rslc RSLC_tab 3 

# resample SLC images using offset models
SLC_resamp_all SLC_tab slc/19990915.slc slc/19990915.slc.par rslc RSLC_tab 4

##################################################################################################
#coregister SLC images using DEM to resample, uses rdc_trans to generate the initial lookup table for interferometry
#generate SLC_tab

mk_tab slc slc slc.par SLC_tab
output tab file: SLC_tab
slc/19990915.slc  slc/19990915.slc.par
slc/19991020.slc  slc/19991020.slc.par

$ SLC_resamp_lt_all
*** /Users/cw/gamma_software/DIFF/scripts/SLC_resamp_lt_all
*** Copyright 2016, Gamma Remote Sensing, v6.6 22-Sep-2016 clw ***
*** Resample a set of SLCs to a common reference SLC using a lookup table and DEM ***

usage: /Users/cw/gamma_software/DIFF/scripts/SLC_resamp_lt_all <SLC_tab> <ref_SLC> <ref_SLC_par> <ref_MLI_par> <DEM_rdc> <MLI_dir> <RSLC_dir> <RSLC_tab> <mode> [rflag] [rlks] [azlks] [rpos] [azpos] [patch] [thres_mli]

    SLC_tab      (input) two column list of SLC filenames and SLC parameter filenames (including paths) (ascii)
    ref_SLC      (input) reference SLC (including path)
    ref_SLC_par  (input) ISP image parameter file of the reference SLC (including path)
    ref_MLI_par  (input) ISP image parameter file of the reference MLI (including path)
    DEM_rdc      (input) DEM in Range-Doppler Coordinates with width and height as described in the MLI_par
    MLI_dir      (input) directory containing MLI parameter files, same the same width and height as DEM_rdc.
    RSLC_dir     directory to contain the resampled SLCs, lookup tables, ISP image parameter files, and log files   
    RSLC_tab     (output) RSLC_tab file for the resampled SLC files  
    mode         processing modes:
                   0: generate lookup table and resample MLI-1 image into the geometry of MLI-2
                   1: refine lookup table based on measured offsets between MLI-1 resampled to the geometry of MLI-2 (optional)
                   2: resample SLC-2 to the geometry of the reference SLC using lookup table    
                   3: create offset parameter files for SLC resampling refinement, measure offsets, and calculate the SLC offset fit polynomials
	           4: resample SLC images using lookup table offsets determined in mode 3 and generate RSLC_tab
    rflag        flag for interative refinement of the resampled SLC:
                   0: resample and measure residual range and azimuth offsets  (no interation, default)
		   1: interatively improve range and azimuth offset polynomials then measure residual offsets
    rlks         number of range looks for initial MLI offset estimation  (enter - for default: 1)
    azlks        number of azimuth looks for initial MLI offset estimation  (enter - for default: 1)
    rpos         center of patch for initial MLI offset estimation range (samples) (enter - for default: image center)
    azpos        center of patch for initial MLI offset estimation azimuth (lines) (enter - for default: image center)
    patch        patch size for initial MLI offset estimate (enter - for default: 1024)
    thres_mli    cross-correlation threshold to accept the MLI offset estimates (enter - for default: 0.2)

    -m           (option) generate MLI-1 in geometry of MLI-2, required for optional refinement of lookup table using MLI-1 and MLI-2
    -b blksz     (option) number of lines/block, only required for spotlight-mode SAR where there are along-track changes in the Doppler centroid, minimum value: 64
    -n npoly     (option) number of terms in the range and azimuth offset polynomials (npoly: 1,3,4,6) (default: 3)
    -t thres     (option) cross-correlation threshold to accept an SLC offset measurement (default: 0.2) 
    -s scale     (option) set image display scale factor (default: 0.8)   
    -e exp       (option) set image display exponent (default: 0.4)
    -c           (option) clean up and delete intermediate data files
    -u           (option) update flag, coregister only new SLC scenes that have not yet been processed

#generate lookup table and resample MLI-1 image into the geometry of MLI-2
SLC_resamp_lt_all SLC_tab slc/19990915.slc slc/19990915.slc.par mli_2_10/19990915.mli.par geo/hector_dem.rdc mli_2_10 rslc_lt RSLC_lt_tab 0 

#refine lookup table based on measured offsets between MLI-1 resampled to the geometry of MLI-2 
SLC_resamp_lt_all SLC_tab slc/19990915.slc slc/19990915.slc.par mli_2_10/19990915.mli.par geo/hector_dem.rdc mli_2_10 rslc_lt RSLC_lt_tab 1 

#resample SLC-2 to the geometry of the reference SLC using lookup table
SLC_resamp_lt_all SLC_tab slc/19990915.slc slc/19990915.slc.par mli_2_10/19990915.mli.par geo/hector_dem.rdc mli_2_10 rslc_lt RSLC_lt_tab 2

#create offset parameter files for SLC resampling refinement, measure offsets, and calculate the SLC offset fit polynomials
SLC_resamp_lt_all SLC_tab slc/19990915.slc slc/19990915.slc.par mli_2_10/19990915.mli.par geo/hector_dem.rdc mli_2_10 rslc_lt RSLC_lt_tab 3

#resample SLC images using lookup table offsets determined in mode 3 and generate RSLC_tab, perform second refinement!
SLC_resamp_lt_all SLC_tab slc/19990915.slc slc/19990915.slc.par mli_2_10/19990915.mli.par geo/hector_dem.rdc mli_2_10 rslc_lt RSLC_lt_tab 4 1

range_offset_polynomial: -0.00030 2.6402e-08 3.8420e-09 0.0000e+00 0.0000e+00 0.0000e+00
azimuth_offset_polynomial: 0.00002 -2.0569e-08 -3.9128e-09 0.0000e+00 0.0000e+00 0.0000e+00
##################################################################################################
#generate MLI images of resampled data, either using offset polynomial, or lookup table
mk_mli_all
*** Copyright 2015, Gamma Remote Sensing, v3.0 23-Jul-2015 clw ***
*** Calculate MLI images for a stack of SLCs and average as an option ***

usage: /home/cw/gamma_software/DIFF/scripts/mk_mli_all <SLC_tab> <MLI_dir> <rlks> <azlks> [sflag] [scale] [exp] [MLI_ave] 

    SLC_tab     (input) two column list of SLC filenames and SLC parameter filenames (including paths) (ascii)
    MLI_dir     directory to contain MLI images and MLI image parameter files
    rlks        range looks for the MLI images
    azlks       azimuth looks for the MLI images
    sflag       MLI image average flag:
                  0: no   (default)
                  1: yes 
    scale       relative intensity display scale factor (default: 0.8)
    exp         display exponent (default: 0.35)
    MLI_ave     MLI average image filename (without path)
	
    -s scale    MLI radiometric scale factor, nominal value for calibrated SCOMPLEX HH and VV data: 1.e-6
    -u          (option) update flag, calculate MLI images only if they do not yet exist

#polynomial offset model SLC imges
mk_mli_all RSLC_tab rmli_2_10 2 10 1 .8 .35 hector_rmli.ave

#calculate MLI of resampled SLC images from lookup-table 
mk_mli_all RSLC_lt_tab rmli_lt_2_10 2 10 1 .8 .35 hector_rmli.ave

##################################################################################################
#Generate differential interferogram, first estimate baselines, either from
#polynomial model or lookup-table resampled SLC data using base_calc script

base_calc
*** /Users/cw/gamma_software/DIFF/scripts/base_calc
*** Copyright 2016, Gamma Remote Sensing, v4.0 26-Nov-2016 clw/uw ***
*** Generate baseline plot and output file with perpendicular baselines and delta_T values ***
*** Generate interferogram table (itab) file specifying SLCs for each interferogram  ***

usage: /Users/cw/gamma_software/DIFF/scripts/base_calc <SLC_tab> <SLC_par> <bperp_file> <itab> <itab_type> [plt_flg] [bperp_min] [bperp_max] [delta_T_min] [delta_T_max] [delta_n_max]

    SLC_tab      (input) two column list of SLC filenames and SLC parameter filenames (including paths) (text)
                   1. SLC filename  (includes path)
                   2. SLC parameter filename (includes path)
    SLC_par      (input) reference SLC parameter filename (include path)
    bperp_file   (output) list of dates, bperp and delta_T for interferogram pairs in the itab (text)
    itab         (output) interferogram table with 4 column format:
                   1. row number in SLC_tab of the reference SLC
                   2. row number in SLC_tab of SLC-2 of the interferogram
                   3. line number in the itab
                   4. flag used to indicate if this interferogram is used in IPTA processing (0:not used  1:used)
    itab_type    itab type (enter - for default):
                   0: single reference (default)
                   1: all pairs
    pltflg       bperp plotting flag (enter - for default):
                   0: none (default)
                   1: output plot in PNG format
                   2: screen output
    bperp_min    minimum magnitude of bperp (m) (default = all, enter - for default)
    bperp_max    maximum magnitude of bperp (m) (default = all, enter - for default)
    delta_T_min  minimum number of days between passes (default = 0, enter - for default)
    delta_T_max  maximum number of days between passes
    delta_n_max  maximum scene number difference between passes

#SLC data interpolated using polynomial
base_calc RSLC_tab rslc/19990915.rslc.par hector.bperp itab 0 1 1

#SLC data interpolated using lookup table
base_calc RSLC_lt_tab rslc_lt/19990915.rslc.par hector.bperp itab 0 1 1

visras.py hector.bmp.png    #Gamma raster image display program for all platforms
eog hector.bperp.png&	    #Linux
preview hector.bperp.png&   #macOS
#

#bperp: 20.9 meters
cat itab
   1    2    1  1

##################################################################################################
#calculate differential interferogram, apply adf filter, and unwrap using mcf, and geocode the output raster image
#use phase_sim to simulate the phase, if you want to use phase_sim_orb set option -o 
mk_diff_2d
*** /Users/cw/gamma_software/DIFF/scripts/mk_diff_2d
*** Copyright 2016, Gamma Remote Sensing, v4.8 3-Jun-2016 clw ***
*** Calculate 2D diff. interferograms using RSLC_tab, itab, DEM and deformation rate in Range-Doppler Coordinates (RDC) ***

usage: mk_diff_2d <RSLC_tab> <itab> <bflag> <DEM_rdc> <def> <mli> <mli_dir> <diff_dir> <rlks> <azlks> <cc_win> [rsflg] [azflg] [mflg] [r_samp] [az_line] [cc_min] [cc_max] 

    RSLC_tab  (input) two column list of resampled SLC filenames and SLC parameter filenames (including paths) (text)
                1. SLC filename  (includes path)
                2. SLC parameter filename (includes path)
    itab      (input) table associating interferogram stack records with pairs of SLC stack records (text)
                1. row number in SLC_tab of the reference SLC 
                2. row number in SLC_tab of SfLC-2 of the interferogram
                3. line number in the itab
                4. flag used to indicate if this interferogram is used for IPTA processing (0:not used  1:used)
    bflag     baseline flag:
                0: use initial baseline in the baseline file
                1: use precision baseline in the baseline file
    DEM_rdc   (input) terrain height in radar range-Doppler coordinates (meters, float, enter - for none)
    def       (input) deformation rate (m/year) (enter - for none)
    mli       (input) reference MLI image with same rlks and azlks as the interferogram used for background
    mli_dir   directory containing MLI images of the coregistered SLCs                      
    diff_dir  differential interferograms after subtraction of simulated unwrapped phase   
    rlks      range looks for interferogram generation
    azlks     azimuth looks for interferogram generation
    cc_win    correlation estimation window size in pixels with linear weighting (default: 3)
    rsflg     range spectral shift filtering flag:
                0: off
                1: on (default)
    azflg     azimuth common-band filtering flag:
                0: off (default)
                1: on
    mflg      initial baseline estimation and refinement flag (enter - for default):
                0: orbit state vector data (default)
                1: orbit state vector data + baseline refinement using 2-D FFT
                2: use existing baseline file
    r_samp    range pixel offset to center of FFT window for baseline refinement using 2-D FFT (enter - for default: image center)
    az_line   azimuth line offset to center of FFT window for baseline refinement using 2-D FFT (enter - for default: image center)
    cc_min    minimum correlation threshold for display (enter - for default: 0.1)
    cc_max    maximum correlation threshold for display (enter - for default: 0.9)

    -s scale  (option) set image display scale factor (default: 1)   
    -e exp    (option) set image display exponent (default: 0.35)
    -c        (option) use of cc_ad rather than cc_wave to estimate correlation
    -o        (option) simulate interferogram phase using orbit state vectors with phase_sim_orb
    -a        (option) when using orbit state vectors for simulation of the phase, add phase calculated from residual baseline obtained by base_ls
    -r        (option) SLC parameter file from the scene used for coregistration, required by phase_sim_orb
    -t        (option) Tandem-X single-pass interferometry mode
    -u        (option) update flag, calculate interferograms only if they do not yet exist

#create shaded relief map from DEM
mapshd geo/hector_dem.rdc 2456 40. 40. - - geo/hector_dem_rdc.shd

#calculate interferogram from SLCs resampled without lookup table
# mk_diff_2d RSLC_tab itab 0 geo/hector_dem.rdc - rmli_2_10/hector_rmli.ave rmli_2_10 diff0_2d 2 10 3 1 0 0

#calculate interferogram from SLCs resampled using lookup table from rdc_trans()
mk_diff_2d RSLC_lt_tab itab 0 geo/hector_dem.rdc - geo/hector_dem_rdc.shd rmli_lt_2_10 diff0_lt_2d 2 10 3 1 0 0 -o

preview diff0_lt_2d/*.bmp

#Interferogram spectral filter with adf, exponent 0.4, window size: 32
mk_adf_2d
*** /Users/cw/gamma_software/DIFF/scripts/mk_adf_2d
*** Copyright 2017, Gamma Remote Sensing, v2.9 31-Aug-2017 clw ***
*** Adaptive filtering (adf) of a set of differential interferograms ***

usage: /Users/cw/gamma_software/DIFF/scripts/mk_adf_2d <RSLC_tab> <itab> <mli> <diff_dir> [cc_win] [adf_exp] [adf_win] [adf_step] 

    RSLC_tab  (input) two column list of coregistered SLC filenames and SLC parameter filenames
             	1. SLC filename  (includes path)
                2. SLC parameter filename (includes path)
    itab      (input) table associating interferogram stack records with pairs of RSLC_tab records (text)		
                1. row number in RSLC_tab of the reference SLC 
             	2. row number in RSLC_tab of SLC-2 of the interferogram
             	3. line number in the itab
                4. flag used to indicate if this interferogram is used for IPTA processing (0:not used  1:used)
    mli       (input) background image with the same dimensions as the interferogram (MLI or shaded relief)
    diff_dir  differential interferogram directory containing *.diff differential interferogram files    
    cc_win    correlation estimation window range and azimuth size (linear weighting) in pixels, (enter - for default: 9)
    adf_exp   exponent parameter for adf interferogram filter, nominal range 0.2-->1.0, (enter - for default: 0.4)
    adf_win   window size for adf filter (enter - for default: 64)
    adf_step  range and azimuth filter step size (1 < adf_step <= adf_win/2, enter - for default: )
    cc_min    minimum correlation threshold for display (enter - for default: 0.1)
    cc_max    maximum correlation threshold for display (enter - for default: 0.9)
    
    -m MLI_dir  (option) use MLI2 as the background image rather than mli image specified on the command line
                MLI_dir is the directory containing the mli images
    -s scale    (option) set image display scale factor (default: 0.7)   
    -e exp      (option) set image display exponent (default: 0.35)
    -u          (option) update flag, filter only interferograms that have not yet been filtered

#Apply adf interferogram phase filter, use MLI image as background
#mk_adf_2d RSLC_tab itab rmli_2_10/hector_rmli.ave diff0_2d 7 .4 32 
#mk_adf_2d RSLC_lt_tab itab rmli_lt_2_10/hector_rmli.ave diff0_lt_2d 7 .4 32

#use shaded relief image as background
mk_adf_2d RSLC_lt_tab itab geo/hector_dem_rdc.shd diff0_lt_2d 7 .4 32 -s 1.2

#unwrap the phase using minimum cost flow (mcf)
mk_unw_2d
*** /Users/cw/gamma_software/DIFF/scripts/mk_unw_2d
*** Copyright 2017, Gamma Remote Sensing, v5.0 31-Aug-2017 clw ***
*** Unwrap the phase of a stack of interferograms using Minimum Cost Flow (MCF) algorithm ***

usage: /Users/cw/gamma_software/DIFF/scripts/mk_unw_2d <RSLC_tab> <itab> <rmli> <diff_dir> [cc_thres] [pwr_thres] [nlks] [npat_r] [npat_az] [mode] [r_init] [az_init] [tri_mode] [unw_mask] [roff] [loff] [nr] [nlines]

    RSLC_tab    (input) two column list of coregistered SLC filenames and SLC parameter filenames (including paths)
                  1. SLC filename  (includes path)
                  2. SLC parameter filename (includes path)
    itab        (input) table associating interferograms with pairs of SLCs listed in the RSLC_tab
                  1. row number in SLC_tab of the reference SLC 
                  2. row number in SLC_tab of SLC-2 of the interferogram
                  3. line number in the itab
                  4. flag used to indicate if this interferogram is to be considered in time-series processing (0:not used  1:used)
    rmli        (input) MLI image derived from the coregistered SLCs with the same pixel spacing and dimensions as the interferogram
    diff_dir    differential interferogram directory containing *.diff differential interferogram files 
    cc_thres    threshold for correlation for creating the unwrapping mask (0.0 --> 1.0) (default: 0.4)
    pwr_thres   threshold for relative intensity for creating the unwrapping mask (0.0 --> 1.0) (default: 0)
    nlks        number of looks in range and azimuth to scale before unwrapping (default: 1)
    npat_r      number of patches in range (default: 1)
    npat_az     number of patches in azimuth (default: 1)
    mode        processing mode:
                  0: unwrap unfiltered data (*.diff)
	          1: (default) unwrap adf filtered data (*.adf.diff) using adf correlation (*.adf.cc)
    r_init      phase reference range offset (default: -)
    az_init     phase reference azimuth offset (default: -)
    tri_mode    MCF triangulation mode:
                  0: filled triangular mesh 
                  1: Delaunay triangulation (default)
    unw_mask1   mask file to specify unwrap region rather than using a mask generated using cc_thres and pwr_thres
                parameters (enter - for none) (8-bit/pixel Sun raster, BMP, or TIFF format) 
    roff        offset to starting range of region to unwrap (default: 0)
    loff        offset to starting line of region to unwrap (default: 0)
    nr          number of range samples of region to unwrap (default(-): width - roff)
    nlines      number of lines of region to unwrap (default(-): total number of lines - loff)
    
    -b rmli_dir  (option) use MLI2 as the background image for display rather than MLI image specified on the command line
    -d diff_tab (option) output a DIFF_tab file containing 2 column list of unwrapped diff. interferograms and delta_T values in decimal days 
    -n          (option) generate DIFF_tab only, no interferometric processing 
    -s scale    (option) set image display scale factor (default: 0.7)   
    -e exp      (option) set image display exponent (default: 0.35)
    -p pscale   (option) set phase scaling for output display, 1 cycle = 2PI/pscale (default: 0.5)
    -m mask     (option) mask to apply to the unwrapped phase (Sun raster or BMP format)
    -u          (option) update flag, unwrap only interferograms that have not yet been unwrapped

#unwrap with CC threshold at 0.1, intensity threshold at 0.01, single patch, reference phase at 100,100, unwrap region covered by the differential interferogram
#leave no extra frame of 0 values to avoid phase unwrapping errors
#mk_unw_2d RSLC_tab itab rmli_2_10/hector_rmli.ave diff0_2d 0.1 0.01 1 1 1 1 100 100 1 -d diff_tab -s .9
mk_unw_2d RSLC_lt_tab itab geo/hector_dem_rdc.shd diff0_lt_2d 0.1 0.01 1 1 1 1 100 100 1 - 90 12 2360 -d diff_tab -s 1.2
disras diff0_lt_2d/19990915_19991020.adf.unw.bmp&

#############################################################################################################
#Terrain geocode the raster image of the unwrapped phase just for visulisation

mk_geo_data
*** /Users/cw/gamma_software/DIFF/scripts/mk_geo_data
*** Copyright 2017, Gamma Remote Sensing, 1.7 29-Jan-2017 clw ***
*** Terrain geocoding of data sets in SAR slant-range coordinates ***

usage: /Users/cw/gamma_software/DIFF/scripts/mk_geo_data <MLI_par> <DEM_seg_par> <gc> <data> <data_geo> <interp_mode> <format> <log>
  MLI_par      (input) MLI image parameter file with same dimensions as input data
  DEM_seg_par  (input) DEM parameter file for output geocoded product 
  gc           (input) geocoding lookup table (map coordinates --> RDC), use refined lookup table if available
  data         (input) input data set with same dimensions as the MLI image
  data_geo     (output) output terrain geocoded data set 
  interp_mode  interpolation mode:
                 0: nearest-neighbor
                 1: bicubic spline (default)
                 2: bicubic-log spline, interpolates log(data)
                 3: bicubic-sqrt spline, interpolates sqrt(data)
                 4: B-spline interpolation (default B-spline degree: 5)
                 5: B-spline interpolation log(x)  (default B-spline degree: 5)
                 6: B-spline interpolation sqrt(x) (default B-spline degree: 5)
                 7: Lanczos interpolation (default Lanczos function order: 5)
                 8: Lanczos interpolation log(x)  (default Lanczos function order: 5)
                 9: Lanczos interpolation sqrt(x) (default Lanczos function order: 5)
  format       input data format:
                 0: FLOAT single precision (4 bytes/value)
		 1: FCOMPLEX (float real, float imaginary) 
		 2: SUN/BMP/TIFF 8 or 24-bit raster image
		 3: UNSIGNED CHAR (1 byte/value)
		 4: SHORT integer (2 bytes/value)
                 5: DOUBLE precision floating point (8 bytes/value)
  log          log file name

  -o           (option) Lanczos function order or B-spline degree (2->9) (enter - default: 5)
	   
#mk_geo_data rmli_2_10/19990915.rmli.par geo/hector_eqa_seg.dem_par geo/hector_1.map_to_rdc diff0_2d/19990915_19991020.adf.unw.bmp geo/19990915_19991020.adf.unw.geo.bmp 0 2 geo/mk_geo_data.log
mk_geo_data rmli_lt_2_10/19990915.rmli.par geo/hector_eqa_seg.dem_par geo/hector_1.map_to_rdc diff0_lt_2d/19990915_19991020.adf.unw.bmp geo/19990915_19991020.adf.unw.geo.bmp 0 2 geo/mk_geo_data.log

#generate KML file for the unwrapped data:
*** /Users/cw/gamma_software/DIFF/scripts/mk_kml
*** Copyright 2016, Gamma Remote Sensing, 1.2 27-Sep-2016 clw ***
*** Generate a Google Earth KML file for an JPEG/PNG/BMP/TIFF rasterformat image in EQA projection ***

usage: /Users/cw/gamma_software/DIFF/scripts/mk_kml <DEM_par> <image> <kml>
  DEM_par (input) DEM parameter file with image dimensions
  image   (input) image in geometry described by DEM_par (JPEG, BMP format)
  kml     (output) kml format output file
  
mk_kml geo/hector_eqa_seg.dem_par geo/19990915_19991020.adf.unw.geo.bmp geo/hector_unw.kml

#convert to displacment along the LOS
mk_dispmap_2d diff_tab geo/hector_dem_rdc.shd rmli_lt_2_10/19990915.rmli.par diff0_lt_2d 1.0 disp_tab 0 -s 1.1 -e .4

#generate look_vector components
look_vector rslc_lt/19990915.rslc.par - geo/hector_eqa_seg.dem_par geo/hector_eqa_seg.dem geo/lv_theta geo/lv_phi
disdt_pwr24 geo/lv_phi geo/hector_map.mli 6284 1 1 0 .01
disdt_pwr24 geo/lv_theta geo/hector_map.mli 6284 1 1 0 .01

#####################################################################
# Command summary: no terrain corrrection for resampling of the SLC

#create processing parameter files from CEOS leaders
ERS_pre_proc CEOS_list ./DELFT raw ERS_pre_proc_1.log proc_list 1
ERS_pre_proc CEOS_list ./DELFT raw ERS_pre_proc_2.log proc_list 2
ERS_pre_proc CEOS_list ./DELFT raw ERS_pre_proc_3.log proc_list 3
ERS_pre_proc CEOS_list ./DELFT raw ERS_pre_proc_4.log proc_list 4

#generate processing list file proc_list for use by ERS_proc_all
ERS_pre_proc CEOS_list ./DELFT raw ERS_pre_proc_6.log proc_list 6

#geocode to space int 2.0833333e-4 degrees (about 22.5m)
mk_geo_radcal mli_2_10/19990915.mli mli_2_10/19990915.mli.par DEM/hector_eqa.dem DEM/hector_eqa.dem_par geo/hector_eqa_seg.dem geo/hector_eqa_seg.dem_par geo hector 2.083333333e-4 0 -p -q -d -z -j
mk_geo_radcal mli_2_10/19990915.mli mli_2_10/19990915.mli.par DEM/hector_eqa.dem DEM/hector_eqa.dem_par geo/hector_eqa_seg.dem geo/hector_eqa_seg.dem_par geo hector 2.083333333e-4 1 -p -q -d -z -j -n 3
mk_geo_radcal mli_2_10/19990915.mli mli_2_10/19990915.mli.par DEM/hector_eqa.dem DEM/hector_eqa.dem_par geo/hector_eqa_seg.dem geo/hector_eqa_seg.dem_par geo hector 2.083333333e-4 2 -p -q -d -z -j -n 3
mk_geo_radcal mli_2_10/19990915.mli mli_2_10/19990915.mli.par DEM/hector_eqa.dem DEM/hector_eqa.dem_par geo/hector_eqa_seg.dem geo/hector_eqa_seg.dem_par geo hector 2.083333333e-4 3 -p -q -d -z -j -n 3

mk_geo_radcal2 mli_2_10/19990915.mli mli_2_10/19990915.mli.par DEM/hector_eqa.dem DEM/hector_eqa.dem_par geo/hector_eqa_seg.dem geo/hector_eqa_seg.dem_par geo hector 2.77777778e-4 0 -p -z -d -j
mk_geo_radcal2 mli_2_10/19990915.mli mli_2_10/19990915.mli.par DEM/hector_eqa.dem DEM/hector_eqa.dem_par geo/hector_eqa_seg.dem geo/hector_eqa_seg.dem_par geo hector 2.77777778e-4 1 -p -d -j -n 3
mk_geo_radcal2 mli_2_10/19990915.mli mli_2_10/19990915.mli.par DEM/hector_eqa.dem DEM/hector_eqa.dem_par geo/hector_eqa_seg.dem geo/hector_eqa_seg.dem_par geo hector 2.77777778e-4 2 -p -d -j -n 3
mk_geo_radcal2 mli_2_10/19990915.mli mli_2_10/19990915.mli.par DEM/hector_eqa.dem DEM/hector_eqa.dem_par geo/hector_eqa_seg.dem geo/hector_eqa_seg.dem_par geo hector 2.77777778e-4 3 -p -d -j -n 3

#create shaded relief map from DEM
mapshd geo/hector_dem.rdc 2456 40. 40. - -  geo/hector_dem_rdc.shd 0

#coregister and resample SLC data 
SLC_resamp_all SLC_tab slc/19990915.slc slc/19990915.slc.par rslc RSLC_tab 0 
SLC_resamp_all SLC_tab slc/19990915.slc slc/19990915.slc.par rslc RSLC_tab 1 
SLC_resamp_all SLC_tab slc/19990915.slc slc/19990915.slc.par rslc RSLC_tab 2 
SLC_resamp_all SLC_tab slc/19990915.slc slc/19990915.slc.par rslc RSLC_tab 3 
SLC_resamp_all SLC_tab slc/19990915.slc slc/19990915.slc.par rslc RSLC_tab 4

#generate MLI images
mk_mli_all RSLC_tab rmli_2_10 2 10 1 .8 .35 hector_rmli.ave

#estimate baseline and generate itab to generate interferograms
base_calc RSLC_tab rslc/19990915.rslc.par hector.bperp itab 0 1 1

#generate differential interferogram
mk_diff_2d RSLC_tab itab 0 geo/hector_dem.rdc - geo/hector_dem_rdc.shd rmli_2_10 diff0_2d 2 10 3 1 0 0

#filter with adaptive non-linear spectral filter: adf
mk_adf_2d RSLC_tab itab geo/hector_dem_rdc.shd diff0_2d 3 .4 32 -u

#unwrap the phase using Minimum Cost Flow (rascc_mask, mcf)
mk_unw_2d RSLC_lt_tab itab geo/hector_dem_rdc.shd diff0_lt_2d 0.1 0.01 1 1 1 1 100 100 1 - 90 12 2360 -d diff_tab -s 1.2
#disras diff0_lt_2d/19990915_19991020.adf.unw.bmp&

#convert unwrapped phase to displacement along the line of sight (LOS), 1 meter/color cycle (dispmap, rasdt_pwr24)
mk_dispmap_2d diff_tab geo/hector_dem_rdc.shd  rmli_lt_2_10/19990915.rmli.par diff0_lt_2d 1.0 disp_tab 0 -s 1.1 -e .4
#output list of displacement data: disp_tab
disdt_pwr24 diff0_lt_2d/19990915_19991020.adf.disp geo/hector_dem_rdc.shd 2456 1 1 0 1.0 1.1 0.35 &

#geocode raster image of the unwrapped phase
mk_geo_data rmli_2_10/19990915.rmli.par geo/hector_eqa_seg.dem_par geo/hector_1.map_to_rdc diff0_2d/19990915_19991020.adf.unw.bmp geo/19990915_19991020.adf.unw.geo.bmp 0 2 geo/mk_geo_data.log

#convert image to BMP format for display with Google Earth
#convert geo/19990915_19991020.adf.unw.geo.ras geo/19990915_19991020.adf.unw.geo.bmp

#Generate KML file
mk_kml geo/hector_eqa_seg.dem_par geo/19990915_19991020.adf.unw.geo.bmp geo/hector_unw.kml

#####################################################################
#####################################################################
# Command summary: Use terrain corrrection for SLC resampling
#####################################################################
#####################################################################
ERS_pre_proc CEOS_list ./DELFT raw ERS_pre_proc_1.log proc_list 1
ERS_pre_proc CEOS_list ./DELFT raw ERS_pre_proc_2.log proc_list 2
ERS_pre_proc CEOS_list ./DELFT raw ERS_pre_proc_3.log proc_list 3
ERS_pre_proc CEOS_list ./DELFT raw ERS_pre_proc_4.log proc_list 4

#Generate processing list file proc_list for use by ERS_proc_all
ERS_pre_proc CEOS_list ./DELFT raw ERS_pre_proc_6.log proc_list 6

#process all the images
ERS_proc_all proc_list raw . slc mli_2_10 2 10 0 6144 10

#Generate table of SLC images using script mk_tab. The SLC_tab is used by the resampling script: SLC_resamp_lt_all
mk_tab slc slc slc.par SLC_tab

#Geocode MLI image toi lat/lon geographic coordinates, pixel spacing is 2.77777778e-4 degrees (about 30m)
#programs: gc_map,  geocode, geocode_back, offset_pwrm, offset_fitm, create_dem_par, gc_map_fine
mk_geo_radcal mli_2_10/19990915.mli mli_2_10/19990915.mli.par DEM/hector_eqa.dem DEM/hector_eqa.dem_par geo/hector_eqa_seg.dem geo/hector_eqa_seg.dem_par geo hector 2.77777778e-4 0 -p -d -j -z -n 3
mk_geo_radcal mli_2_10/19990915.mli mli_2_10/19990915.mli.par DEM/hector_eqa.dem DEM/hector_eqa.dem_par geo/hector_eqa_seg.dem geo/hector_eqa_seg.dem_par geo hector 2.77777778e-4 1 -p -q -d -z -j -n 3
mk_geo_radcal mli_2_10/19990915.mli mli_2_10/19990915.mli.par DEM/hector_eqa.dem DEM/hector_eqa.dem_par geo/hector_eqa_seg.dem geo/hector_eqa_seg.dem_par geo hector 2.77777778e-4 2 -p -q -d -z -j -n 3
mk_geo_radcal mli_2_10/19990915.mli mli_2_10/19990915.mli.par DEM/hector_eqa.dem DEM/hector_eqa.dem_par geo/hector_eqa_seg.dem geo/hector_eqa_seg.dem_par geo hector 2.77777778e-4 3 -p -d -j -z -n 3

mk_geo_radcal2 mli_2_10/19990915.mli mli_2_10/19990915.mli.par DEM/hector_eqa.dem DEM/hector_eqa.dem_par geo/hector_eqa_seg.dem geo/hector_eqa_seg.dem_par geo hector 2.77777778e-4 0 -p -z -d -j
mk_geo_radcal2 mli_2_10/19990915.mli mli_2_10/19990915.mli.par DEM/hector_eqa.dem DEM/hector_eqa.dem_par geo/hector_eqa_seg.dem geo/hector_eqa_seg.dem_par geo hector 2.77777778e-4 1 -p -d -j -n 3
mk_geo_radcal2 mli_2_10/19990915.mli mli_2_10/19990915.mli.par DEM/hector_eqa.dem DEM/hector_eqa.dem_par geo/hector_eqa_seg.dem geo/hector_eqa_seg.dem_par geo hector 2.77777778e-4 2 -p -d -j -n 3
mk_geo_radcal2 mli_2_10/19990915.mli mli_2_10/19990915.mli.par DEM/hector_eqa.dem DEM/hector_eqa.dem_par geo/hector_eqa_seg.dem geo/hector_eqa_seg.dem_par geo hector 2.77777778e-4 3 -p -d -j -n 3

#Create shaded relief map from DEM, program mapshd.
mapshd geo/hector_dem.rdc 2456 40. 40. - -  geo/hector_dem_rdc.shd

#Resample SLC(s) to create differential interferograms
#reference SLC  slc/19990915.slc
#reference scene MLI parameter file used to determine DEM geometry parameters: mli_2_10/19990915.mli.par
#DEM resampled into MLI radar geometry: geo/hector_dem.rdc
#resampled SLC directory: rslc_lt
#resampled SLC_tab: RSLC_lt_tab
#programs: rdc_trans, offset_pwr, offset_pwrm, offset_fitm, offset_fit, SLC_interp_lt, multi_look
SLC_resamp_lt_all SLC_tab slc/19990915.slc slc/19990915.slc.par mli_2_10/19990915.mli.par geo/hector_dem.rdc mli_2_10 rslc_lt RSLC_lt_tab 0 -m
SLC_resamp_lt_all SLC_tab slc/19990915.slc slc/19990915.slc.par mli_2_10/19990915.mli.par geo/hector_dem.rdc mli_2_10 rslc_lt RSLC_lt_tab 1
SLC_resamp_lt_all SLC_tab slc/19990915.slc slc/19990915.slc.par mli_2_10/19990915.mli.par geo/hector_dem.rdc mli_2_10 rslc_lt RSLC_lt_tab 2
SLC_resamp_lt_all SLC_tab slc/19990915.slc slc/19990915.slc.par mli_2_10/19990915.mli.par geo/hector_dem.rdc mli_2_10 rslc_lt RSLC_lt_tab 3
SLC_resamp_lt_all SLC_tab slc/19990915.slc slc/19990915.slc.par mli_2_10/19990915.mli.par geo/hector_dem.rdc mli_2_10 rslc_lt RSLC_lt_tab 4 1

#generate shaded relief image of DEM in RDC coordiantes
mapshd geo/hector_dem.rdc 2456 40. 40. - - geo/hector_dem_rdc.shd 0

#Generate MLI images from the stack of resampled SLCs:  RSLC_lt_tab
#output MLI directory: rmli_lt_2_10
#range looks: 2  azimuth looks: 10
#sum images mode selected. sum image file: hector_rmli.ave
#intensity scale: 0.9, exponent: 0.35
mk_mli_all RSLC_lt_tab rmli_lt_2_10 2 10 1 .9 .35 hector_rmli.ave

#Estimate baseline and generate itab to generate interferograms. The itab file is used by mk_diff_2d 
#Single reference mode selected
#baseline plot generated (PNG format): hector.bperp.png
#Minimum baseline length (m): 1.0
base_calc RSLC_lt_tab rslc_lt/19990915.rslc.par hector.bperp itab 0 1 1

#Generate differential interferogram using DEM resampled into MLI image coordinates: phase_sim_orb, SLC_diff_intf
#Option -o means that orbits rather than polynomial baseline model is used to estimate the topographic phase. 
#The baseline is estimated for each point in the scene.
mk_diff_2d RSLC_lt_tab itab 0 geo/hector_dem.rdc - geo/hector_dem_rdc.shd rmli_lt_2_10 diff0_lt_2d 2 10 3 1 0 0 -o
rasmph_pwr diff0_lt_2d/19990915_19991020.diff rmli_lt_2_10/hector_rmli.ave 2456 1 1 0 1 1 .8 0.35 1 diff0_lt_2d/19990915_19991020_mli.diff.bmp
dis2ras diff0_lt_2d/19990915_19991020.diff.bmp diff0_lt_2d/19990915_19991020_mli.diff.bmp&

#Filter differential interferogram with ISP program adf
#adf exponent paramter set to 0.4  window size: 32x32
#adf correlation estimation window size: 7
mk_adf_2d RSLC_lt_tab itab geo/hector_dem_rdc.shd diff0_lt_2d 7 .4 32 -s 1.2

#Unwrap the phase of the differential interferogram using minimum cost flow phase unwrapping (ISP programs: rascc_mask, mcf)
#reference point is at range pixel: 100  azimuth line: 100
#rascc_mask threshold set to 0.01 power and 0.1 correlation coefficient 
#output list of unwrapped interferograms: diff_tab
mk_unw_2d RSLC_lt_tab itab geo/hector_dem_rdc.shd diff0_lt_2d 0.3 0.005 1 1 1 1 100 100 1 -d diff_tab -s 1.2

#Calculate LOS displacement map from the Unwrapped phase and generate a 24-bit raster image using the displacement data
#1.0 meter/color cycle. Output list of displacements in the disp_tab
mk_dispmap_2d diff_tab geo/hector_dem_rdc.shd rmli_lt_2_10/19990915.rmli.par diff0_lt_2d 1.0 disp_tab 0 -s 1.1 -e .4

#Geocode raster image of the unwrapped phase into Lat/Lon geometric coordinates
mk_geo_data rmli_lt_2_10/19990915.rmli.par geo/hector_eqa_seg.dem_par geo/hector_1.map_to_rdc diff0_lt_2d/19990915_19991020.adf.unw.bmp geo/19990915_19991020.adf.unw.geo.bmp 0 2 geo/mk_geo_data.log

#Generate KML file for the unwrapped phase
mk_kml geo/hector_eqa_seg.dem_par geo/19990915_19991020.adf.unw.geo.bmp geo/hector_unw.kml

#To display in Google Earth, open the KML file located in the geo directory: geo/hector_unw.kml

#generate look_vector components
look_vector rslc_lt/19990915.rslc.par - geo/hector_eqa_seg.dem_par geo/hector_eqa_seg.dem geo/lv_theta geo/lv_phi
disdt_pwr24 geo/lv_phi geo/hector_map.mli 3874 1 1 0 .01
disdt_pwr24 geo/lv_theta geo/hector_map.mli 3874 1 1 0 .01
