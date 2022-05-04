# adjointModel
This is a case study to combine the source oriented and receptor source apportionment methods based on the measurements during Feb.7 to Feb.11, 2021.

## Preparation
Now we have a long queue to download the meteo data from ECMWF. The precessing tools in Polyphemus need the meteo data saved by days. To shorten the queueing time, we download all the data for the calculation period once. Then use the [grib_filter](https://confluence.ecmwf.int/display/ECC/grib_filter) in the grib tools to split the whole data into days.

Install the grib api tools
> sudo apt-get install libgrib-api-tools

## Meteorological data preparation
* Meteorological data: the adjoint model should run backward in time. Considering the measurements duration, we decided to prepare the meteorological data covering whole January until Feb. 12, 2021. 

To split the file, in "rules_file"
> write "[centre]_[date].grib";

Then, run the grib_filter to split the data file
> grib_filter rules_file SingleLevel_20210101.grib 

where SingleLevel_20210101.grib is the downloaded file

[Xygrib](https://opengribs.org/en/) can be used to overview the generated data files, and check the file information, including area and numbers of grid points.

* [preProcessing.py](preProcessing.py): generate the model-ready ground data and meteorological data for Polair3d model.

## Emission data preparation
* Temporal profile: use the newly developed [high resolution emission profiles in EDGAR](https://edgar.jrc.ec.europa.eu/dataset_temp_profile). The hourly profile is generated by the combination of weekly factors and daily factors.
* [profiles_MEIC.m](src/preprocessing/profiles_MEIC.m): generate the model ready emissions for five sectors (i.e. industry, residential combustion, traffic, power and agriculture) by combining three emission inventories, Beijing-Tianjing-Hebei fine emission inventory, MEIC emission inventory for China and MIX for Asia, which can be obtained from the [MEIC database](http://meicmodel.org/).   

## Run adjoint model
* [adjointFullRun.m](src/preprocessing/adjoint/adjointFullRun.m): generate the meteorological data for the adjoint model, and run the adjoint model based on the time period of the [measurement data](forwardRun/rawData/measurements/Xian.csv)

## Polyphemus modifications
* [meteo.cpp](src/modifiedPolyphemus/meteo.cpp): The original meteo.cpp in Polyphemus has been modified in order to use ERA5 data, which only 1 hour accumulation
https://confluence.ecmwf.int/pages/viewpage.action?pageId=155337784 
> scons meteo_parallel openmp=yes flag_openmp=mp

* [polair3dAdjoint.cpp](src/modifiedPolyphemus/polair3dAjoint/polair3dAdjoint.cpp): modified code for adjoint run, and it should be put in the 'processing/photochemistry' folder of the original Polyphemus codes.
* [included files](src/modifiedPolyphemus/polair3dAjoint/includeModels/): the modified codes in this folder should be put in the "/include/models" folder o the original Polyphemus codes.
