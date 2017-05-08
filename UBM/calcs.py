"""
These are where the water balance calcuations take place
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import arcpy
from arcpy.sa import *



def calc_avail_water(path, path2, months='',years=''):

    arcpy.env.workspace = path
    arcpy.env.overwriteOutput = True

    if months == '':
        months = [1,12]
    if years == '':
        years = [2004,2014]

    for y in range(years[0], years[1]+1): #set years converted here
        for m in range(months[0], months[1]+1): #set months converted here
            my = str(y)+str(m).zfill(2)
            newdn = 'AVWT' + my
            rain = Raster('RAIN'+ my +'SUM')
            melt = Raster('SNML'+ my +'SUM')
            calc = rain + melt
            calc.save(path2+newdn)
            print(newdn)


def UBM_calc(results, field_cap, wilt_point, T_soil_water, geol_k, avail_water, pet, months='', years=''):

    arcpy.env.overwriteOutput = True

    if months == '':
        months = [1, 12]
    if years == '':
        years = [2004, 2014]

    av_soil_water = field_cap
    for y in range(years[0], years[1] + 1):  # set years converted here
        for m in range(months[0], months[1] + 1):  # set months converted here
            my = str(y) + str(m).zfill(2)
            av_wtr = Raster(avail_water + 'AVWT' + my)
            pet_rast = Raster(pet + 'PET' + my)

            av_soil_water = av_wtr + av_soil_water

            # Eq 1
            av_recharge1 = "T_soil_water - field_cap"
            recharge1 = "Con(eval(av_recharge1) > geol_k, geol_k, eval(av_recharge1))"
            runoff1 = "(av_soil_water - T_soil_water) + Con(eval(av_recharge1) > geol_k, eval(av_recharge1) - geol_k, 0)"
            # Eq2
            av_recharge2 = "av_soil_water - field_cap"
            recharge2 = "Con(eval(av_recharge2) > geol_k, geol_k, eval(av_recharge2))"
            runoff2 = "Con(eval(av_recharge2) > geol_k, eval(av_recharge2) - geol_k, 0)"
            
            # Eq3 recharge3 = 0 runoff3 = 0 aet = pet_rast
            av_evap = "av_soil_water - pet_rast"
            aet3 = "Con(eval(av_evap) <= wilt_point, wilt_point, pet_rast)"
            # Eq4 recharge3 = 0 runoff3 = 0 aet = 0

            # Order of if/then is Eq 1, Eq 4, Eq 2, Eq 3
            recharge = Con(av_soil_water > T_soil_water, eval(recharge1),
                           Con(av_soil_water < wilt_point, 0,
                               Con((av_soil_water < T_soil_water) & (av_soil_water > field_cap), eval(recharge2),
                                   Con((av_soil_water > wilt_point) & (av_soil_water < field_cap), 0, 0))))
            recharge.save(results + 'rec' + my)

            runoff = Con(av_soil_water > T_soil_water, eval(runoff1),
                         Con(av_soil_water < wilt_point, 0,
                             Con((av_soil_water < T_soil_water) & (av_soil_water > field_cap), eval(runoff2),
                                 Con((av_soil_water > wilt_point) & (av_soil_water < field_cap), 0, 0))))
            runoff.save(results + 'run' + my)

            aet = Con(av_soil_water > T_soil_water, pet_rast,
                      Con(av_soil_water < wilt_point, 0,
                          Con((av_soil_water < T_soil_water) & (av_soil_water > field_cap), eval(aet3),
                              Con((av_soil_water > wilt_point) & (av_soil_water < field_cap), pet_rast, 0))))
            aet.save(results + 'aet' + my)

            av_soil_water = av_soil_water - runoff - recharge - aet
            av_soil_water = Con(av_soil_water > 0.0, av_soil_water, 0.0)
            av_soil_water.save(results + 'asw' + my)
            print(my)

def summarize(path, code, statistics='MEAN'):
    arcpy.env.workspace = path
    arcpy.env.overwriteOutput = True

    rlist = arcpy.ListRasters(code+"*") #pick all files from raw data folder of a data type
    # arcpy sa functions that summarize the daily data to monthly data
    calc = CellStatistics(rlist, statistics_type = statistics, ignore_nodata="DATA")
    outnm = code+"_"+statistics
    calc.save(outnm)
    print(outnm)


def monthly_to_yearly(path, code, yearRange='', statistics='SUM'):
    if yearRange=='':
        yearRange = [2004,2014]
    arcpy.env.workspace = path
    arcpy.env.overwriteOutput = True
    for y in range(yearRange[0],yearRange[1]+1): #set years converted here
        ylist = arcpy.ListRasters(code+str(y)+"*") #pick all files from raw data folder of a data type
        calc = CellStatistics(ylist, statistics_type = statistics, ignore_nodata="DATA")
        outnm = 'y'+code+str(y)
        desc = arcpy.Describe(calc)
        print(outnm)
        calc.save(outnm)

if __name__ == '__main__':
    main()
