# Nov. 12, 2021, Xiaole Zhang
# Run the Polair3D model

import os
from pathlib import Path
from datetime import date, timedelta

startDay = date(2021, 2, 1)

# Landuse data
# check whether the ground data dir exists or not
# if not, create one
groundFile = Path("./data/ground/LUC-glcf.bin")
if not groundFile.is_file():
    groundPath = Path("./data/ground/")
    groundPath.mkdir(parents=True, exist_ok=True)
    os.system("luc-glcf config/general.cfg config/luc-glcf.cfg")
    
# covert the landuse data
groundZhangFile = Path("./data/ground/LUC-glcf-zhang.bin")
if not groundZhangFile.is_file():
    os.system("luc-convert  config/general.cfg config/glcf_to_zhang.cfg")

# roughness
roughnessFile = Path("./data/ground/Roughness-glcf.bin")
if not roughnessFile.is_file():
    os.system("roughness config/general.cfg config/roughness.cfg")

# Meteo data
meteoPath = Path("./data/meteo/")
meteoPath.mkdir(parents=True, exist_ok=True)

# Kz_TM data folder
kzTMPath = Path("./data/meteo/Kz_TM")
kzTMPath.mkdir(parents=True, exist_ok=True)

# dep data folder
depPath = Path("./data/dep")
depPath.mkdir(parents=True, exist_ok=True)

# emission data folder
emissionPath = Path("./data/emissions")
emissionPath.mkdir(parents=True, exist_ok=True)

# results folder
resultsPath = Path("./results")
resultsPath.mkdir(parents=True, exist_ok=True)

for daysN in range(0,12):
    durations = timedelta(days=daysN)
    currentDay = durations+startDay
    nextDay = currentDay+timedelta(days=1)
    print('Current Day:' + currentDay.isoformat())
    
    # process the meteo data
    os.system("meteo config/general.cfg config/meteo.cfg "+currentDay.strftime('%Y-%m-%d'))
    
    # estimate Kz
    os.system("Kz config/general.cfg config/meteo.cfg "+currentDay.strftime('%Y-%m-%d')+" "+nextDay.strftime('%Y-%m-%d'))
    
    # estimate Kz_TM
    os.system("Kz_TM config/general.cfg config/meteo.cfg "+currentDay.strftime('%Y-%m-%d')+" "+nextDay.strftime('%Y-%m-%d'))

    # estimate dep
    os.system("dep_aerosol config/general.cfg config/dep.cfg "+currentDay.strftime('%Y-%m-%d')+" "+nextDay.strftime('%Y-%m-%d'))
