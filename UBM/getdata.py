"""
These are data input download and prep scripts. They download and massage the data for the UBM calculations (calc.py)

"""
from __future__ import absolute_import, division, print_function, unicode_literals
import time
import urllib
import urllib2
import re
import glob
import os
import arcpy
from arcpy.sa import *


def get_modis(tiles, save_path, months='', years=''):
    """The following script automatically retrieves monthly MODIS16 hdf file from the ntsg website.

    :param tiles: Tile number in format h##v##; based on grid from https://modis-land.gsfc.nasa.gov/MODLAND_grid.html
    :param save_path: name of output file name
    :param months: months of interest; defaults to [1,12]
    :param years: years of interest; defaults to [2000,2015]
    :return: saves files in outpath
    """


    from bs4 import BeautifulSoup
    if months == '':
        months = [1, 12]
    if years == '':
        years = [2000, 2015]

    mons = [str(i).zfill(2) for i in range(months[0], months[1] + 1)]
    yrs = [str(i) for i in range(years[0], years[1] + 1)]

    for tile in tiles:
        for yr in yrs:
            for m in mons:
                base_url = "http://files.ntsg.umt.edu/data/NTSG_Products/MOD16/MOD16A2_MONTHLY.MERRA_GMAO_1kmALB/"

                dir_path = "Y{:}/M{:}/".format(yr, m)
                url = base_url + dir_path
                soup = BeautifulSoup(urllib2.urlopen(url), "lxml")
                hdf_name = soup.find_all('', {
                    'href': re.compile('MOD16A2.A{:}M{:}.{:}.105'.format(yr, m, tile), re.IGNORECASE)})
                files = urllib.urlretrieve(url + hdf_name[0].text, save_path + hdf_name[0].text)
                print(save_path + hdf_name[0].text)
                time.sleep(0.5)


def get_file_list(save_path):
    return glob.glob(os.path.join(save_path, '*.105*.hdf'))


def reproject_modis(files, save_path, data_type, proj=102003):
    """Iterates through MODIS files in a folder reprojecting them.

    Takes the crazy MODIS sinusoidal projection to a user defined projection.

    Args:
        files: list of file paths of MODIS hdf files; created using files = glob.glob(os.path.join(save_path, '*.105*.hdf'))
        save_path: folder to store the reprojected files
        data_type: type of MODIS16 data being reprojected; options are 'ET','PET','LE', and 'PLE'
        proj: projection of output data by epsg number; default is nad83 zone 12
    Returns:
        Reprojected MODIS files

    ..notes:
    The EPSG code for NAD83 Zone 12 is 26912.
    The EPSG code for Albers Equal Area is 102003
    http://files.ntsg.umt.edu/data/NTSG_Products/MOD16/MOD16_global_evapotranspiration_description.pdf
    https://modis-land.gsfc.nasa.gov/MODLAND_grid.html
    """
    import pymodis
    # dictionary to designate a directory
    datadir = {'ET': '/ET/', 'PET': '/PET/', 'LE': '/LE/', 'PLE': '/PLE/'}
    # dictionary to select layer from hdf file that contains the datatype
    matrdir = {'ET': [1, 0, 0, 0], 'LE': [0, 1, 0, 0], 'PET': [0, 0, 1, 0], 'PLE': [0, 0, 0, 1]}

    # check for file folder and make it if it doesn't exist
    if not os.path.exists(save_path + datadir[data_type]):
        os.makedirs(save_path + datadir[data_type])
        print('created {:}'.format(save_path + datadir[data_type]))

    for f in files:
        year = f.split('\\')[1].split('.')[1][1:5]  # parse year from hdf filename
        month = f.split('\\')[1].split('.')[1][-2:]  # parse month from hdf filename
        v = f.split('\\')[1].split('.')[2][-2:]  # parse v (cell coordinate) from hdf filename
        h = f.split('\\')[1].split('.')[2][1:3]  # parse h (cell coordinate) from hdf filename
        pref = os.path.join(save_path + datadir[data_type] + 'A' + year + 'M' + month + 'h' + h + 'v' + v)
        convertsingle = pymodis.convertmodis_gdal.convertModisGDAL(hdfname=f, prefix=pref,
                                                                   subset=matrdir[data_type],
                                                                   res=1000, epsg=proj)
        # [ET,LE,PET,PLE]
        try:
            convertsingle.run()
        except:
            print('A' + year + 'M' + month + 'h' + h + 'v' + v + ' failed!')
            pass


def clip_and_fix(path, outpath, data_type, area=''):
    """Clips raster to Utah's Watersheds and makes exception values null.

    Args:
        path: folder of the reprojected MODIS files
        outpath: ESRI gdb to store the clipped files
        data_type: type of MODIS16 data being reprojected; options are 'ET','PET','LE', and 'PLE'
        area: path to polygon used to clip tiles

    """
    # Check out the ArcGIS Spatial Analyst extension license
    arcpy.CheckOutExtension("Spatial")

    arcpy.env.workspace = path
    arcpy.env.overwriteOutput = True

    if area == '':
        area = 'H:/GIS/Calc.gdb/WBD_UT'

    arcpy.env.mask = area
    arcpy.CheckOutExtension("spatial")
    for rast in arcpy.ListRasters():
        calc = SetNull(arcpy.Raster(rast) > 32760, arcpy.Raster(rast))
        calc.save(outpath + data_type + rast[1:5] + rast[6:8] + 'h' + rast[10:11] + 'v' + rast[13:14])
        print(outpath + data_type + rast[1:5] + rast[6:8] + 'h' + rast[10:11] + 'v' + rast[13:14])


def merge_rasts(path, data_type='ET', monthRange='', yearRange='', outpath=''):
    """Mosaics (merges) different MODIS cells into one layer.


    """
    if monthRange == '':
        monthRange = [1, 12]
    if yearRange == '':
        yearRange = [2000, 2015]
    if outpath == '':
        outpath = path

    arcpy.env.workspace = path
    outCS = arcpy.SpatialReference('NAD 1983 UTM Zone 12N')
    for y in range(yearRange[0], yearRange[-1] + 1):  # set years converted here
        for m in range(monthRange[0], monthRange[-1] + 1):  # set months converted here
            nm = data_type + str(y) + str(m).zfill(2)
            rlist = []
            for rast in arcpy.ListRasters(nm + '*'):
                rlist.append(rast)
            try:
                arcpy.MosaicToNewRaster_management(rlist, outpath, nm + 'c', outCS, \
                                                   "16_BIT_UNSIGNED", "1000", "1", "LAST", "LAST")

                print(path + nm + 'c')
            except:
                print(nm + ' failed!')
                pass


def scale_modis(path, out_path, scaleby=10000.0, data_type='ET', monthRange=[1, 12], yearRange=[2000, 2014]):
    """

    :param path: directory to unconverted modis tiles
    :param out_path: directory to put output in
    :param scaleby: scaling factor for MODIS data; default converts to meters/month
    :param data_type: type of MODIS16 data being scaled; used for file name; options are 'ET','PET','LE', and 'PLE'
    :param monthRange: range of months to process data
    :param yearRange: range of years to process data
    :return:
    """
    arcpy.CheckOutExtension("spatial")

    for y in range(yearRange[0], yearRange[-1] + 1):  # set years converted here
        for m in range(monthRange[0], monthRange[-1] + 1):  # set months converted here
            nm = data_type + str(y) + str(m).zfill(2)
            calc = Divide(nm + 'c', scaleby)
            calc.save(out_path + nm)


def untar(filepath, outfoldername='.', compression='r', deletesource=False):
    """
    Given an input tar archive filepath, extracts the files.
    Required: filepath -- the path to the tar archive
    Optional: outfoldername -- the output directory for the files; DEFAULT is directory with tar archive
              compression -- the type of compression used in the archive; DEFAULT is 'r'; use "r:gz" for gzipped archives
              deletesource -- a boolean argument determining whether to remove the archive after extraction; DEFAULT is false
    Output:   filelist -- the list of all extract files
    """
    import tarfile

    with tarfile.open(filepath, compression) as tfile:
        filelist = tfile.getnames()
        tfile.extractall(path=outfoldername)

    if deletesource:
        try:
            os.remove(filepath)
        except:
            raise Exception("Could not delete tar archive {0}.".format(filepath))

    return filelist


def ungz(filepath, compression='rb', deletesource=False):
    """
    Given an input gz archive filepath, extracts the files.
    Required: filepath -- the path to the tar archive
    Optional: outfoldername -- the output directory for the files; DEFAULT is directory with tar archive
              compression -- the type of compression used in the archive; DEFAULT is 'r'; use "r:gz" for gzipped archives
              deletesource -- a boolean argument determining whether to remove the archive after extraction; DEFAULT is false
    Output:   filelist -- the list of all extract files
    """

    import gzip

    with gzip.open(filepath, compression) as f:
        outF = open(filepath[:-3], 'wb')
        outF.write(f.read())
        f.close()
        outF.close()
    if deletesource:
        try:
            os.remove(filepath)
        except:
            raise Exception("Could not delete gz archive {0}.".format(filepath))

    return filepath[:-3]


def replace_hdr_file(hdrfile):
    """
    Replace the .hdr file for a .bil raster with the correct data for Arc processing
    Required: hdrfile -- filepath for .hdr file to replace/create
    Output:   None
    """
    # hdr file replacment string
    HDRFILE_STRING = "byteorder M\nlayout bil\nnbands 1\nnbits 16\nncols 6935\nnrows 3351\n\
    ulxmap -124.729583333331703\nulymap 52.871249516804028\nxdim 0.00833333333\nydim 0.00833333333\n"
    with open(hdrfile, 'w') as o:
        o.write(HDRFILE_STRING)


def get_snodas(out_dir, months='', years=''):
    """Downloads daily SNODAS data from ftp.  This is slow.

    :param out_dir: directory to store downloaded SNODAS zip files
    :param months: months desired for download
    :param years: years desired for download
    :return: saved zip files in out_dir

    .. note:
    Use polaris: http://nsidc.org/data/polaris/
    """
    import ftplib

    if months == '':
        months = [1, 12]
    if years == '':
        years = [2000, 2015]

    monnames = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    mons = [str(i).zfill(2) + "_" + monnames[i - 1] for i in range(months[0], months[1] + 1)]

    yrs = [str(i) for i in range(years[0], years[1] + 1)]

    for yr in yrs:
        for m in mons:
            ftp_addr = "sidads.colorado.edu"
            ftp = ftplib.FTP(ftp_addr)
            ftp.login()

            dir_path = "pub/DATASETS/NOAA/G02158/masked/" + yr + "/" + m + "/"
            ftp.cwd(dir_path)
            files = ftp.nlst()

            for f in files:
                if len(f) > 4:
                    save_file = open(out_dir + "/" + f, 'wb')
                    ftp.retrbinary("RETR " + f, save_file.write)
                    save_file.close()
                    print(f)
            ftp.close()


def rename_polaris_snodas(path):
    prodcode = {'us_ssmv11038wS__A': 'SPAT', 'us_ssmv11044bS__T': 'SNML', 'us_ssmv11050lL00T': 'SPSB',
                'us_ssmv11034tS__T': 'SWEQ', 'us_ssmv01025SlL00': 'RAIN', 'us_ssmv01025SlL01': 'SNOW',
                'us_ssmv11036tS__T': 'SNOD', 'us_ssmv11039lL00T': 'BSSB'}

    for filename in os.listdir(path):
        if filename.startswith("us_ssmv"):
            code = prodcode[filename[0:17]]
            yrsrt = filename.find('TNATS') + 5
            yr = filename[yrsrt:yrsrt + 4]
            mo = filename[yrsrt + 4:yrsrt + 6]
            dy = filename[yrsrt + 6:yrsrt + 8]
            try:
                os.rename(os.path.join(path, filename), os.path.join(path, code + yr + mo + dy + filename[-4:]))
            except:
                pass


def snow_summary(code, scalingFactor, statistics="SUM", outcellsize='1000', monthRange='', yearRange='',
                 path="H:/GIS/SNODAS/SNWDS/", outpath="H:/GIS/SNODAS.gdb/", area=''):
    """
    summarizes daily SNODAS data to monthly values

    INPUT
    -----
    code = text; prefix of dataset to use; choices are 'RAIN','SWEQ','SNOD','SPAT','BSSB','SNML', or 'SPSB'
    scalingFactor = float; table 1 at http://nsidc.org/data/docs/noaa/g02158_snodas_snow_cover_model/
    statistics = text; from arcpy sa CellStatistics; choices are MEAN, MAJORITY, MAXIMUM, MEDIAN, MINIMUM, MINORITY,
                    RANGE, STD, SUM, or VARIETY
    monthRange = len 2 list; begin and end month of data you wish to analyze
    yearRange = len 2 list; bengin and end year of data you wish to analyze
    path = directory where raw geoTiffs are located
    outpath = directory where final data will be stored

    OUTPUT
    ------
    projected and scaled monthly rasters

    """
    if monthRange == '':
        months = [1, 12]
    if yearRange == '':
        years = [2000, 2015]

    g = {}
    arcpy.env.workspace = path
    arcpy.env.overwriteOutput = True
    if area == '':
        area = 'H:/GIS/Calc.gdb/WBD_UT'
    # arcpy.env.mask = area

    statstype = {'MEAN': 'AVG', 'MAJORITY': 'MAJ', 'MAXIMUM': 'MAX', 'MEDIAN': 'MED', 'MINIMUM': 'MIN',
                 'MINORITY': 'MNR',
                 'RANGE': 'RNG', 'STD': 'STD', 'SUM': 'SUM', 'VARIETY': 'VAR'}

    for y in range(yearRange[0], yearRange[1] + 1):  # set years converted here
        for m in range(monthRange[0], monthRange[1] + 1):  # set months converted here
            g[code + str(y) + str(m).zfill(2)] = []  # this defines the dictionary key based on data type month and year
            for name in sorted(
                    glob.glob(path + code + '*.tif')):  # pick all tiff files from raw data folder of a data type
                rast = os.path.basename(name)
                if rast[0:4] == code and int(rast[4:8]) == y and int(rast[8:10]) == m:
                    g[code + str(y) + str(m).zfill(2)].append(rast)  # create a list of rasters for each month
                else:
                    pass
            if len(g[code + str(y) + str(m).zfill(2)]) > 0:
                # print(g[code+str(y)+str(m).zfill(2)])
                # ifnull = 'in_memory/ifnull'
                # arcpy sa functions that summarize the daily data to monthly data
                cellstats = CellStatistics(g[code + str(y) + str(m).zfill(2)], statistics_type=statistics,
                                           ignore_nodata="DATA")
                div = Divide(cellstats, scalingFactor)  # scale factor, converts to kg/m2 10 then to m 0.001
                calc = Con(div < 0.0, 0.0, div)  # remove negative and null values
                ifnull = Con(IsNull(calc), 0, calc)  # remove null
                # WKID 102039
                outCS = arcpy.SpatialReference(102039)  # change coordinate units to m for spatial analysis
                # define save path for file
                outnm = outpath + rast[0:4] + str(y).zfill(2) + str(m).zfill(2) + statstype[statistics]
                memoryFeature = "in_memory/myMemoryFeature"
                # memoryFeature = outnm
                arcpy.ProjectRaster_management(ifnull, memoryFeature, outCS, 'BILINEAR', outcellsize,
                                               'WGS_1984_(ITRF00)_To_NAD_1983', '#', '#')
                # Execute ExtractByMask to clip snodas data to Utah watersheds
                extrc = arcpy.sa.ExtractByMask(memoryFeature, area)
                extrc.save(outnm)
                print(outnm)
                arcpy.Delete_management("in_memory")


def totalavg(code, statistics="MEAN", monthRange=[1, 12], yearRange=[2003, 2016],
             path="H:/GIS/SNODAS/SNODASproj.gdb/", outpath="H:/GIS/SNODAS/SNODASproj.gdb/"):
    """Summarizes daily raster data into monthly data.

    INPUT
    -----
        code = string with four letters represting data type to summarize (example 'BSSB')
        statistics = how data will be summarized; defaults to monthly averages; options are
            ['MEAN','MAJORITY','MAXIMUM','MEDIAN','MINIMUM','MINORITY','RANGE','STD','SUM','VARIETY']
            Most common are 'MEAN','MEDIAN', and 'SUM'
            These are inputs that will be used in the ArcPy CellStatistics function.
            See http://pro.arcgis.com/en/pro-app/tool-reference/spatial-analyst/cell-statistics.htm for documentation
        monthRange = beginning and end months of summary statistics
        yearRange = beginning and end years of summary statistics
        path = location of geodatabase of data to summarize
        outpath = location of geodatabase where output data should be stored
    OUTPUT
    ------
        summary raster(s) stored in outpath

    """
    g = {}
    statstype = {'MEAN': 'AVG', 'MAJORITY': 'MAJ', 'MAXIMUM': 'MAX', 'MEDIAN': 'MED', 'MINIMUM': 'MIN',
                 'MINORITY': 'MNR',
                 'RANGE': 'RNG', 'STD': 'STD', 'SUM': 'SUM', 'VARIETY': 'VAR'}
    arcpy.env.workspace = path
    arcpy.env.overwriteOutput = True

    # iterate over month range set here; default is 1 to 12 (Jan to Dec)
    for m in range(monthRange[0], monthRange[1] + 1):

        # this defines the dictionary key based on data type, month, and year
        g[code + '0000' + str(m).zfill(2)] = []

        # pick all tiff files from raw data folder of a data type
        for rast in arcpy.ListRasters():
            yrrng = range(yearRange[0], yearRange[1] + 1)  # set years converted here

            # create a list of rasters with the right code and month and year
            if rast[0:4] == code and int(rast[4:8]) in yrrng and int(rast[8:10]) == m:
                g[code + '0000' + str(m).zfill(2)].append(rast)  # create a list of rasters for each month
            else:
                pass
        if len(g[code + '0000' + str(m).zfill(2)]) > 0:
            # arcpy sa functions that summarize the daily data to monthly data
            calc = CellStatistics(g[code + '0000' + str(m).zfill(2)], statistics_type=statistics, ignore_nodata="DATA")
            calc.save(code + '0000' + str(m).zfill(2) + statstype[statistics])
            print(code + '0000' + str(m).zfill(2) + statstype[statistics])

if __name__ == '__main__':
    main()
