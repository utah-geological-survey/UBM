import arcpy

if arcpy.CheckExtension("Spatial") == "Available":
    arcpy.AddMessage("Checking out Spatial")
    arcpy.CheckOutExtension("Spatial")
else:
    arcpy.AddError("Unable to get spatial analyst extension")
    arcpy.AddMessage(arcpy.GetMessages(0))
    sys.exit(0)

from arcpy.sa import *


results = arcpy.GetParameterAsText(0)
field_cap = arcpy.GetParameterAsText(1)
wilt_point = Raster(arcpy.GetParameterAsText(2))
T_soil_water = Raster(arcpy.GetParameterAsText(3))
geol_k = Raster(arcpy.GetParameterAsText(4))
avail_water = arcpy.GetParameterAsText(5)
pet = arcpy.GetParameterAsText(6)
starting_soil_water = Raster(arcpy.GetParameterAsText(7))
begin_month = int(arcpy.GetParameterAsText(8))
end_month = int(arcpy.GetParameterAsText(9))
begin_year = int(arcpy.GetParameterAsText(10))
end_year = int(arcpy.GetParameterAsText(11))

months = [begin_month, end_month]
years = [begin_year, end_year]

arcpy.env.workspace = results
arcpy.env.overwriteOutput = True


av_soil_water = starting_soil_water


for y in range(years[0], years[1] + 1):  # set years converted here
    for m in range(months[0], months[1] + 1):  # set months converted here
        my = str(y) + str(m).zfill(2)
        av_wtr = Raster(avail_water + '/AVWT' + my)
        pet_rast = Raster(pet + '/PET' + my)

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
        av_evap = "av_soil_water - wilt_point"
        aet3 = "Con(eval(av_evap) <= wilt_point, eval(av_evap), pet_rast)"
        # Eq4 recharge3 = 0 runoff3 = 0 aet = 0
        # Order of if/then is Eq 1, Eq 4, Eq 2, Eq 3

        recharge = Con(av_soil_water > T_soil_water, eval(recharge1),
                       Con(av_soil_water < wilt_point, 0,
                           Con((av_soil_water < T_soil_water) & (av_soil_water > field_cap), eval(recharge2),
                               Con((av_soil_water > wilt_point) & (av_soil_water < field_cap), 0, 0))))
        recharge.save(results + '/rec' + my)

        runoff = Con(av_soil_water > T_soil_water, eval(runoff1),
					 Con(av_soil_water < wilt_point, 0,
						 Con((av_soil_water < T_soil_water) & (av_soil_water > field_cap), eval(runoff2),
							 Con((av_soil_water > wilt_point) & (av_soil_water < field_cap), 0, 0))))
        runoff.save(results + '/run' + my)

        aet = Con(av_soil_water > T_soil_water, pet_rast,
				  Con(av_soil_water < wilt_point, 0,
					  Con((av_soil_water < T_soil_water) & (av_soil_water > field_cap), eval(aet3),
						  Con((av_soil_water > wilt_point) & (av_soil_water < field_cap), pet_rast, 0))))
        aet.save(results + '/aet' + my)

        av_soil_water = av_soil_water - runoff - recharge - aet
        av_soil_water = Con(av_soil_water > 0.0, av_soil_water, 0.0)
        av_soil_water.save(results + '/asw' + my)
        arcpy.AddMessage('You made a BM for ' + my + '!')
