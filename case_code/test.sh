#!/bin/bash

# 算例数据(每次需修改)
cityName="Shijiazhuang" 
initialYear="2021"
initialMonth="10"
numberOfMonths="6"
days="170"

siteX="114.5358"
siteY="38.0247" 

Species="SJZ"

# 存储静态数据
Date_min="$initialYear-$initialMonth-02_00-00-00"
preProcessing_date="$initialYear, $initialMonth, 2"
rawDataDirectory="rawData\/meteo\/${cityName}\/"

declare -A coords_x
declare -A coords_y
coords_x["Shijiazhuang"]="x_min = 104.5	Delta_x = 0.10	Nx = 200"
coords_x_ecmwf["Shijiazhuang"]="x_min = 104.5	Delta_x = 0.25	Nx = 81"
coords_x["Beijing"]="x_min = 105	Delta_x = 0.10	Nx = 200"
coords_x["Chongqing"]="x_min = 96.5	Delta_x = 0.10	Nx = 200"

coords_y["Shijiazhuang"]="y_min = 28	Delta_y = 0.10	Ny = 200"
coords_y_ecmwf["Shijiazhuang"]="y_min = 28	Delta_y = 0.25	Ny = 81"
coords_y["Beijing"]="y_min = 30	Delta_y = 0.10	Ny = 200"
coords_y["Chongqing"]="y_min = 19.8	Delta_y = 0.10	Ny = 200"


mkdir -p 'forwardRun/rawData/meteo/'$cityName
cp ../src_template/fetchEcmwfData.ipynb forwardRun/rawData/meteo/fetchEcmwfData.ipynb
cp ../src_template/fetchEcmwfData.py forwardRun/rawData/meteo/fetchEcmwfData.py
cp ../src_template/hybrid_coefficients.dat 'forwardRun/rawData/meteo/'$cityName'/hybrid_coefficients.dat'
cp ../src_template/rules_file 'forwardRun/rawData/meteo/'$cityName'/rules_file'
cp ../src_template/rules_file_pressure 'forwardRun/rawData/meteo/'$cityName'/rules_file_pressure'

sed -i 's/City_temp/'$cityName'/g' 'forwardRun/rawData/meteo/'$cityName'/rules_file'
sed -i 's/City_temp/'$cityName'/g' 'forwardRun/rawData/meteo/'$cityName'/rules_file_pressure'

sed -i 's/City_temp/'$cityName'/g' 'forwardRun/rawData/meteo/fetchEcmwfData.py'
sed -i 's/initialYear_temp/'$initialYear'/g' 'forwardRun/rawData/meteo/fetchEcmwfData.py'
sed -i 's/initialMonth_temp/'$initialMonth'/g' 'forwardRun/rawData/meteo/fetchEcmwfData.py'
sed -i 's/numberOfMonths_temp/'$numberOfMonths'/g' 'forwardRun/rawData/meteo/fetchEcmwfData.py'

cp ../src_template/forwardRun/preProcessing.py forwardRun/preProcessing.py
sed -i 's/date_default/'"${preProcessing_date}"'/g' 'forwardRun/preProcessing.py'
sed -i 's/days_default/'$days'/g' 'forwardRun/preProcessing.py'

cp -r ../src_template/forwardRun/config forwardRun/config
sed -i 's/Date_min_default/'$Date_min'/g' 'forwardRun/config/caseConfig.cfg'
sed -i 's/coords_x_default/'"${coords_x[${cityName}]}"'/g' 'forwardRun/config/caseConfig.cfg'
sed -i 's/coords_y_default/'"${coords_y[${cityName}]}"'/g' 'forwardRun/config/caseConfig.cfg'

sed -i 's/Date_min_default/'$Date_min'/g' 'forwardRun/config/case-data.cfg'

sed -i 's/date_default/'"$Date_min"'/g' 'forwardRun/config/general.cfg'
sed -i 's/coords_x_default/'"${coords_x[${cityName}]}"'/g' 'forwardRun/config/general.cfg'
sed -i 's/coords_y_default/'"${coords_y[${cityName}]}"'/g' 'forwardRun/config/general.cfg'

sed -i 's/Database_meteo_default/'"$rawDataDirectory"'/g' 'forwardRun/config/meteo.cfg'
sed -i 's/coords_x_ecmwf_default/'"${coords_x_ecmwf[${cityName}]}"'/g' 'forwardRun/config/meteo.cfg'
sed -i 's/coords_y_ecmwf_default/'"${coords_y_ecmwf[${cityName}]}"'/g' 'forwardRun/config/meteo.cfg'

sed -i 's/siteX_default/'"${siteX}"'/g' 'forwardRun/config/source.dat'
sed -i 's/siteY_default/'"${siteY}"'/g' 'forwardRun/config/source.dat'
sed -i 's/Species_default/'"${Species}"'/g' 'forwardRun/config/source.dat'

cd ./forwardRun
python3 preProcessing.py > preProcessing.log 2>&1
