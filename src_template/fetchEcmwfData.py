# Nov. 11, 2021, Xiaole Zhang
# Download meteo data from ecmwf 
# Dec, 29, 2021, Xiaole Zhang
# fix the date bugs when combine the files
# Apr. 28, 2023, Xiaole Zhang
# Get meteo data for Beijing ARG
# May 07, 2023, Xiaole Zhang
# fix the misalignment of the data
import cdsapi
from datetime import date, timedelta, datetime
from dateutil.relativedelta import relativedelta
import os
from pathlib import Path
import calendar 

# covert number to formated string
def number2String(f):
    return '%02d' % f

c = cdsapi.Client() 

# [up, left, down, right]
regionName = 'City_temp'
regionBound = {'Beijing': [50, 105, 30, 127],\
               'Xian': [44.25, 97.5, 24.25, 119.5], \
              'Wuhan': [40.5, 104.5, 20.5, 124.5], \
              'Chongqing':[39.8, 96.5 ,19.8, 116.5], \
              'Shijiazhuang':[48, 104.5, 28, 124.5], \
	      'Langfang': [50, 105, 30, 127]}

# data for Wuhan
initialYear = initialYear_temp
initialMonth = initialMonth_temp
numberOfMonths = numberOfMonths_temp
initialDate = datetime(initialYear, initialMonth, 1, 0) 

#targetYear = 2021
#targetMonths = [1,2]
#monthsList = list(map(number2String, targetMonths))

targetDays = range(1,32)
daysList = list(map(number2String, targetDays))

# to download the meteo day by month to meet the volume requirement
for currentMonthNum in range(numberOfMonths):
    currentDay = initialDate+ relativedelta(months=currentMonthNum)#date(targetYear, currentMonth, 1) 
    initalDay = currentDay 
    
    year = '%4d' % currentDay.year
    month = '%02d' % currentDay.month
    day = '%02d' % currentDay.day
    print('Current Day:' + currentDay.isoformat())
    
    # download single level data
    fileName = regionName + '/'+'SingleLevel_'+currentDay.strftime('%Y%m%d')+'.grib'
    filePath = Path(regionName + '/')
    filePath.mkdir(parents=True, exist_ok=True)
    
    c.retrieve( 
        'reanalysis-era5-single-levels', 
        { 
            'product_type': 'reanalysis', 
            'format': 'grib', 
            'variable': [ 
                '2m_temperature', 'boundary_layer_height', 'eastward_turbulent_surface_stress', 
                'evaporation', 'high_cloud_cover', 'medium_cloud_cover', 
                'northward_turbulent_surface_stress', 'skin_temperature', 'surface_sensible_heat_flux', 
                'surface_solar_radiation_downwards', 'surface_pressure', 'convective_precipitation','large_scale_precipitation','volumetric_soil_water_layer_1', 
            ], 
            'year': year, 
            'month': month,
            'day': daysList, 
            'time': [ 
                '00:00', '01:00', '02:00', 
                '03:00', '04:00', '05:00', 
                '06:00', '07:00', '08:00', 
                '09:00', '10:00', '11:00',
                '12:00', '13:00', '14:00',
                '15:00', '16:00', '17:00',
                '18:00', '19:00', '20:00',
                '21:00', '22:00', '23:00',
            ], 
            'area': regionBound[regionName], 
        }, 
        fileName) 
    
    # split the data to daily files
    ruleFile =regionName+'/rules_file'
    os.system("grib_filter "+ruleFile+" "+fileName)
    
    ###############################################
    # download 3d data
    fileName3d =  regionName + '/'+'PressureLevel_'+currentDay.strftime('%Y%m%d')+'.grib'
    
    c.retrieve( 
        'reanalysis-era5-pressure-levels', 
        { 
        'product_type': 'reanalysis', 
        'variable': [ 
        'specific_cloud_ice_water_content', 'specific_cloud_liquid_water_content', 'specific_humidity', 
        'temperature', 'u_component_of_wind', 'v_component_of_wind', 
        ], 
        'pressure_level': [ 
        '500', '550', '600', 
        '650', '700', '750', 
        '775', '800', '825', 
        '850', '875', '900', 
        '925', '950', '975', 
        '1000', 
        ], 
        'year': year, 
        'month': month,
        'day': daysList, 
        'time': [ 
        '00:00', '01:00', '02:00', 
        '03:00', '04:00', '05:00', 
        '06:00', '07:00', '08:00', 
        '09:00', '10:00', '11:00',
        '12:00', '13:00', '14:00',
        '15:00', '16:00', '17:00',
        '18:00', '19:00', '20:00',
        '21:00', '22:00', '23:00', 
        ], 
        'area': regionBound[regionName], 
        'format': 'grib', 
        }, 
        fileName3d) 
    
    # split the data to daily files
    ruleFile =regionName+'/rules_file_pressure'
    os.system("grib_filter "+ruleFile+" "+fileName3d)
    
    ###########################
    # combine the daily 2d and 3d data files into single files and delete the splitted data
    # need to install the grib api tools
    # sudo apt-get install libgrib-api-tools
    if(currentMonthNum!=0):      
        currentDay = initalDay + timedelta(days=-1)
        fileNameAll = regionName + '/'+'ECMWF-'+currentDay.strftime('%Y%m%d')+'.grb'
        currentDay = initalDay + timedelta(hours=-6)
        fileName2dFc = regionName + '/' + 'ecmf_fc_'+currentDay.strftime('%Y%m%d_%-H')+'.grib'
        # ">>" should be used instead of ">", otherwise the input file will be erased
        os.system("cat "+fileNameAll+" "+fileName2dFc+" >> "+fileNameAll)

    dayRange = calendar.monthrange(initalDay.year, initalDay.month)
    for day in range(0, dayRange[1]):
        currentDay = initalDay + timedelta(days=day) 
        print('Current Day:' + currentDay.isoformat())
        fileName3d = regionName + '/' + 'ecmf_'+currentDay.strftime('%Y%m%d')+'_pressureLevels.grib'
        fileNameAll = regionName + '/'+'ECMWF-'+currentDay.strftime('%Y%m%d')+'.grb'
        if(currentMonthNum==0 or day!=0):
            os.system("cat "+fileNameAll+" > "+fileNameAll)
        
        for hour in [-6, 6, 18]:
            currentDay = initalDay + timedelta(days=day) + timedelta(hours = hour)
            fileName2dFc = regionName + '/' + 'ecmf_fc_'+currentDay.strftime('%Y%m%d_%-H')+'.grib'
            os.system("cat "+fileNameAll+" "+fileName2dFc+" >> "+fileNameAll)
            if(hour!=18):
               os.system("rm "+fileName2dFc)  
        
        for hour in range(0, 24):
            currentDay = initalDay + timedelta(days=day) + timedelta(hours = hour)
            fileName2dAn = regionName + '/' + 'ecmf_an_'+currentDay.strftime('%Y%m%d_%-H')+'.grib'
            os.system("cat "+fileNameAll+" "+fileName2dAn+" >> "+fileNameAll)
            os.system("rm "+fileName2dAn)
            
        os.system("cat "+fileNameAll+" "+fileName3d+" >> "+fileNameAll)
        os.system("rm "+fileName3d) 
        
    currentDay = initalDay + timedelta(days=day+1)
    fileNameAll = regionName + '/'+'ECMWF-'+currentDay.strftime('%Y%m%d')+'.grb'
    os.system("cat "+fileNameAll+" > "+fileNameAll)
    currentDay = initalDay + timedelta(days=day+1) + timedelta(hours = -6)
    fileName2dFc = regionName + '/' + 'ecmf_fc_'+currentDay.strftime('%Y%m%d_%-H')+'.grib'
    os.system("cat "+fileNameAll+" "+fileName2dFc+" >> "+fileNameAll)
    os.system("rm "+fileName2dFc)
        
     
   #fileNameAll = regionName + '/'+'ECMWF-'+currentDay.strftime('%Y%m%d')+'.grb'
   #os.system("cat "+fileName+" "+fileName3d+" > "+fileNameAll)
   #os.system("rm "+fileName+" "+fileName3d)
