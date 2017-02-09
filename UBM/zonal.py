import arcpy
import pandas as pd
import wellapplication as wa
import numpy as np
import matplotlib.pyplot as plt


def calcvols(tablegdb, searchstr, source, variable, mult = 1.0):
    """Calculates volume of water per zone in ac-ft. Uses output from zone_gdb. Created pandas DataFrame.

    :param tablegdb: Path to file geodatabase in which tables are stored
    :param searchstr: Search wildcard to select a subset of input tables; use astrix (*) for any string after search string
    :param source: Designate name of data source field
    :param variable: Designate name of data variable; ex. 'runoff'
    :param mult: Multiplier to adjust values of input
    :return: pandas DataFrame of zonal values in ac-ft
    """
    arcpy.env.workspace = tablegdb
    tables = arcpy.ListTables(searchstr)
    fields = arcpy.ListFields(tables[0])

    fieldlist = [field.name for field in fields]

    f = {}
    for table in tables:
        f[table] = pd.DataFrame(arcpy.da.TableToNumPyArray(table ,fieldlist))
    g = pd.concat(f)
    g.reset_index(inplace=True)

    g['YearMonth'] = g['level_0'].apply(lambda x: str(x)[-9:-3] if str(x)[-3:] == 'SUM' else str(x)[-6:] ,1)
    g['HUC_10'] = g['HUC_12'].apply(lambda x: str(x)[:-2] ,1)
    g['HUC_08'] = g['HUC_12'].apply(lambda x: str(x)[:-4] ,1)
    g.drop(['level_0', 'level_1', 'OBJECTID', 'ZONE_CODE'], axis=1, inplace=True)
    g['SOURCE'] = source
    g['variable'] = variable
    g['volume_m_cubed'] = g['MEAN'] * g['AREA'] * mult
    g['volume_acft'] = g['volume_m_cubed'] * 0.000810714

    return g

def get_zone(rast, z_Name, Zonal_HUCS, Zone_field):
    dsc = arcpy.Describe(rast)
    nm = dsc.baseName
    arcpy.sa.ZonalStatisticsAsTable(Zonal_HUCS, Zone_field, rast, z_Name + "/z_" + nm, "DATA", "ALL")
    print("z_" + nm)

def zone_gdb(indata, z_Name, Zonal_HUCS, Zone_field, wildcard='*'):
    """Creates geodatabase of tables summarizing zonal statistics from rasters.

    :param indata: Geodatabase containing rasters to summarize with zonal statistics
    :param z_Name: Output geodatabase for zonal statistics tables
    :param Zonal_HUCS: Polygons used to conduct zonal statistics
    :param Zone_field: Field in Zonal_HUCS to summarize zonal statistics
    :param wildcard: Wildcard to select subset of rasters
    :return: zonal statistics tables (1 for each raster)
    """
    arcpy.env.workspace = indata
    arcpy.CheckOutExtension("Spatial")
    arcpy.env.overwriteOutput = True

    for rast in arcpy.ListRasters(wildcard):
        dsc = arcpy.Describe(rast)
        nm = dsc.baseName
        arcpy.sa.ZonalStatisticsAsTable(Zonal_HUCS, Zone_field, rast, z_Name + "/z_" + nm, "DATA", "ALL")
        print("z_" + nm)


def runModel(mrg, geo_k=''):
    huc12list = mrg['HUC_12'].unique()

    mrg['avail_water'] = np.nan
    mrg['avail_rech'] = np.nan
    mrg['aet'] = np.nan
    mrg['runoff'] = np.nan
    mrg['recharge'] = np.nan
    mrg['eqt'] = np.nan

    grp = {}

    for h in huc12list:

        grp[h] = mrg[mrg['HUC_12'] == h]
        grp[h] = grp[h][~grp[h].index.duplicated(keep='first')]

        soil_max = grp[h]['total soil moisture'].mean()
        field_cap = grp[h]['field capacity'].mean()
        if geo_k == '':
            geo_k = grp[h]['conductivity'].mean()

        wilt_pnt = grp[h]['wilting point'].mean()
        dates = pd.date_range(start=grp[h].index.min(), end=grp[h].index.max(), freq='MS')


        for i in dates:

            inwater = grp[h].ix[i, 'snow and rain']
            pet = grp[h].ix[i, 'PET']

            if i == dates[0]:
                avail_water = inwater + field_cap
            elif i.month == 1:
                avail_water = inwater + grp[h].ix[pd.datetime(i.year - 1, 12, 1), 'avail_water']
            else:
                avail_water = inwater + grp[h].ix[pd.datetime(i.year, i.month - 1, 1), 'avail_water']

            if avail_water > soil_max:

                avail_rech = soil_max - field_cap
                grp[h].ix[i, 'aet'] = pet

                if avail_rech > geo_k:
                    grp[h].ix[i, 'eqt'] = 1.1
                    grp[h].ix[i, 'runoff'] = (avail_water - soil_max) + (avail_rech - geo_k)
                    grp[h].ix[i, 'recharge'] = geo_k
                else:
                    grp[h].ix[i, 'eqt'] = 1.2
                    grp[h].ix[i, 'runoff'] = avail_water - soil_max
                    grp[h].ix[i, 'recharge'] = avail_rech
            elif (avail_water < soil_max) and (avail_water > field_cap):

                avail_rech = avail_water - field_cap
                grp[h].ix[i, 'aet'] = pet

                if avail_rech > geo_k:
                    grp[h].ix[i, 'eqt'] = 2.1
                    grp[h].ix[i, 'runoff'] = avail_rech - geo_k
                    grp[h].ix[i, 'recharge'] = geo_k
                else:
                    grp[h].ix[i, 'eqt'] = 2.2
                    grp[h].ix[i, 'runoff'] = 0
                    grp[h].ix[i, 'recharge'] = avail_rech
            elif (avail_water > wilt_pnt) and (avail_water < field_cap):
                avail_rech = 0
                grp[h].ix[i, 'eqt'] = 3
                grp[h].ix[i, 'runoff'] = 0
                grp[h].ix[i, 'recharge'] = 0
                grp[h].ix[i, 'aet'] = pet
            elif avail_water < wilt_pnt:
                avail_rech = 0
                grp[h].ix[i, 'eqt'] = 4
                grp[h].ix[i, 'runoff'] = 0
                grp[h].ix[i, 'recharge'] = 0
                grp[h].ix[i, 'aet'] = 0
            else:
                pass

            grp[h].ix[i, 'avail_rech'] = avail_rech

            av_water_0 = avail_water - grp[h].ix[i, 'runoff'] - grp[h].ix[i, 'recharge'] - grp[h].ix[i, 'aet']
            if av_water_0 > 0:
                grp[h].ix[i, 'avail_water'] = av_water_0
            else:
                grp[h].ix[i, 'avail_water'] = 0


    if len(huc12list) > 1:
        df = pd.concat([grp[h] for h in huc12list])
    else:
        df = grp[h]

    return df

def process_huc(HUC):
    HUC = [str(i) for i in HUC]
    huc10 = ["'{:}'".format(i[:-2]) for i in HUC]
    huc10 = list(set(huc10))
    return HUC, huc10

def get_usgs_data(SITE, strtDT=''):
    # Grab USGS data for comparison
    if strtDT == '':
        strtDT = '2003-01-01'

    nw = wa.nwis('dv', SITE, 'sites', startDT=strtDT)
    data = nw.data
    label = nw.sites.station_nm[0].title()
    data['cfd'] = data.value * 86400  # cfs to cfd
    if isinstance(data.index, pd.core.index.MultiIndex):
        data.index = data.index.droplevel(0)
    # convert to acft by converting from monthly sum of cfd to acft-mo
    acft = data['cfd'].groupby(pd.TimeGrouper('MS', label='left')).sum() * 2.29569E-5
    acft = acft.to_frame()
    acft['month'] = acft.index.month

    return acft, label

def get_UBM_data(HUC, engine, table):
    HUC, huc10 = process_huc(HUC)
    quer = "SELECT HUC_12,YearMonth,volume_acft FROM ubm.{:} WHERE SOURCE = '{:}' AND HUC_10 IN({:}) AND variable IN('{:}')"
    dataset = 'UBM'
    variable = 'recharge'
    Urc = pd.read_sql_query(sql=quer.format(table, dataset, ','.join(huc10), variable), con=engine)
    variable = 'runoff'
    Urn = pd.read_sql_query(sql=quer.format(table, dataset, ','.join(huc10), variable), con=engine)

    Urun = Urn[Urn['HUC_12'].isin(HUC)]
    Urun = Urun.rename(columns={'volume_acft': 'runoff_acft'})
    Urun.set_index(['HUC_12', 'YearMonth'], inplace=True)
    Urec = Urc[Urc['HUC_12'].isin(HUC)]
    Urec = Urec.rename(columns={'volume_acft': 'recharge_acft'})
    Urec.set_index(['HUC_12', 'YearMonth'], inplace=True)

    UBM = pd.concat([Urun, Urec], axis=1)
    UBM['dt'] = pd.to_datetime(UBM.index.get_level_values(1), errors='coerce', format='%Y%m')

    UBM.reset_index(inplace=True)
    UBM.dropna(inplace=True)

    UBMgrp = UBM.groupby(['dt']).sum()
    UBMgrp['month'] = UBMgrp.index.month

    return UBM, UBMgrp

def get_model_inputs(HUC,engine,table):
    HUC, huc10 = process_huc(HUC)
    UBM = get_UBM_data(HUC, engine, table)[0]
    # This section pulls the individual model inputs and plots them (third figure)
    quer = "SELECT HUC_12,YearMonth,volume_acft,SOURCE,AREA,variable FROM ubm.{:} WHERE HUC_10 IN({:}) AND SOURCE IN({:})"
    sources = "'Surrgo','State Geologic Maps','SNODAS','MODIS16'"

    chk = pd.read_sql_query(sql=quer.format(table,','.join(huc10), sources), con=engine)
    chk = chk[chk['HUC_12'].isin(HUC)]
    chk['dt'] = pd.to_datetime(chk.YearMonth, errors='coerce', format='%Y%m')

    piv = pd.pivot_table(chk, index=['HUC_12', 'dt'], columns='variable', values='volume_acft')
    pv = pd.pivot_table(chk, index=['HUC_12'], columns='variable', values='volume_acft')
    pv.drop([u'evapotranspiration', u'precip as rain', u'snowmelt'], inplace=True, axis=1)
    piv.reset_index(inplace=True)
    pv.reset_index(inplace=True)
    mrg1 = pd.merge(piv, pv, on='HUC_12')
    mrg1.reset_index(inplace=True)

    areas = chk.drop(['YearMonth', 'SOURCE', 'variable', 'volume_acft', 'dt'], axis=1)
    areas.drop_duplicates(inplace=True)
    mrg = pd.merge(mrg1, areas, on='HUC_12')
    mrg.dropna(inplace=True)

    # mrg.set_index(['HUC_12','dt'],inplace=True)

    mrg['incoming_water'] = mrg[u'precip as rain'] + mrg[u'snowmelt']

    mrg.set_index(['dt'], inplace=True)
    mrg.drop('index', inplace=True, axis=1)
    mrg['YearMonth'] = [str(x.year) + str(x.month).zfill(2) for x in mrg.index]
    mrg1 = pd.merge(UBM, mrg, on=['HUC_12', 'YearMonth'])
    mrg1.set_index(['dt'], inplace=True)
    return mrg1

def merge_huc(mrg):
    mrg_grp = mrg.groupby(level=0).agg([np.sum, np.mean])
    # mrg_grp = mrg.groupby(['dt']).agg(np.sum)

    mrg_grp.drop(['AREA', u'precip as rain', u'snowmelt', 'porosity'], inplace=True, axis=1)
    if isinstance(mrg_grp.index, pd.core.index.MultiIndex):
        mrg_grp.index = mrg_grp.index.droplevel(0)
    mrg_grp.drop([('evapotranspiration', 'mean'), ('incoming_water', 'mean'),
                  ('field capacity', 'sum'), ('total soil moisture', 'sum'),
                  ('wilting point', 'sum'), ('conductivity', 'sum')], axis=1, inplace=True)
    mrg_grp.columns = mrg_grp.columns.droplevel(-1)
    return mrg_grp

def plotfits(HUC, SITE, engine, fileloc, table):
    from matplotlib.backends.backend_pdf import PdfPages

    HUC = [str(i) for i in HUC]
    huc10 = ["'{:}'".format(i[:-2]) for i in HUC]
    huc10 = list(set(huc10))

    pdf = PdfPages(fileloc + str(HUC[0])[:-2] + '.pdf')

    ubm, UBMgrp = get_UBM_data(HUC, engine, table)
    acft, label = get_usgs_data(SITE)

    ubm.reset_index(inplace=True)
    ubm.dropna(inplace=True)

    # get monthly average values (plot 2)
    UBMmon = UBMgrp.groupby('month').mean()
    acgp = acft.groupby('month').mean()

    plt.figure()
    plt.plot(UBMgrp.index, UBMgrp.recharge_acft, 'o-', label='Modeled Recharge')
    plt.plot(UBMgrp.index, UBMgrp.runoff_acft, 'o-', label='Modeled Runoff')
    plt.plot(acft.index, acft['cfd'], 'o-', label='Measured Surface Discharge')
    plt.ylabel('Discharge (ac-ft)')
    plt.legend()
    plt.title(label + ' ' + str(SITE))
    plt.grid()
    pdf.attach_note(label)
    pdf.savefig()
    plt.close()

    plt.figure()
    plt.plot(UBMmon.index, UBMmon.recharge_acft, 'o-', label='Modeled Recharge')
    plt.plot(UBMmon.index, UBMmon.runoff_acft, 'o-', label='Modeled Runoff')
    plt.plot(acgp.index, acgp, 'o-', label='Measured Surface Discharge')
    plt.ylabel('Discharge (ac-ft)')
    plt.legend()
    plt.title(label + ' ' + str(SITE))
    plt.grid()
    pdf.attach_note(label)
    pdf.savefig()
    plt.close()

    mrg = get_model_inputs(HUC, engine, table)
    mrg_grp = merge_huc(mrg)


    mrg_grp.plot()
    plt.text('1/1/2009', mrg_grp['total soil moisture'].median(), 'Max Soil Capacity')
    plt.text('1/1/2009', mrg_grp['wilting point'].median(), 'Wilting Point')
    plt.text('1/1/2009', mrg_grp['field capacity'].median(), 'Field Capacity')
    plt.text('1/1/2009', mrg_grp['conductivity'].median(), 'Geologic K')
    plt.legend(loc=1)
    plt.ylabel('ac-ft')
    plt.title(label + ' ' + str(SITE))
    pdf.attach_note(label)
    pdf.savefig()
    plt.close()
    pdf.close()

    return ubm, acft, UBMmon, acgp, mrg


