function [dataTmp] = readMEIC(dataFolder,sector)
%read the MEIC data
%  Mean value between 2013 and 2015
%% meic grid
xllcorner = 70.0;
yllcorner = 10.0;
ncols = 320;
nrows = 200; 
cellsize = 0.25;

xMeic = xllcorner:cellsize:(xllcorner+(ncols-1)*cellsize);
yMeic = yllcorner:cellsize:(yllcorner+(nrows-1)*cellsize);
[xMeic, yMeic] = meshgrid(xMeic, yMeic);

 r = 6371229; % meter
 areas = (cos(yMeic/180*pi)*r.*cellsize/180*pi)*(cellsize/180*pi.*r);
 %%
switch sector
    case {'INDUSTRY'}
        sectorMEIC = 'industry';
    case {'RESIDENTIAL'}
        sectorMEIC = 'residential';
    case {'TRANSPORT'}
        sectorMEIC = 'transportation';
    case {'POWER'}
        sectorMEIC = 'power';
    case {'agriculture'}
        sectorMEIC = 'agriculture';
end

emi = 0;
for years = 2013:2015
    filename = fullfile(dataFolder, num2str(years), [num2str(years) '_0_' sectorMEIC '_PM25.asc'] );
    dataTmp =  importdata(filename,' ',6);
    emi = emi + dataTmp.data;
end

% convert from ton/year to microg/m2/s
dataTmp = flipud(emi)/3./areas*10^3*10^3*10^6/(365*24*3600);
end

