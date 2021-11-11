# adjointModel
This is a case study to combine the source oriented and receptor source apportionment methods based on the measurements durng Feb.7 to Feb.11, 2021.

## Preparation
Now we have a long queue to download the meteo data from ECMWF. The precessing tools in Polyphemus need the meteo data saved by days. To shorten the queueing time, we download all the data for the calculation period once. Then use the [grib_filter](https://confluence.ecmwf.int/display/ECC/grib_filter) in the grib tools to split the whole data into days.

Install the grib api tools
> sudo apt-get install libgrib-api-tools

In "rules_file"
> write "[centre]_[date].grib";

Then, run the grib_filter to split the data file
> grib_filter rules_file SingleLevel_20210101.grib 
where SingleLevel_20210101.grib is the downloaded file

## Data preparation
* Meteorological data: the adjoint model should run backward in time. Considering the measurements duration, we decided to prepare the meteorological data covering whole January until Feb. 12, 2021.  

