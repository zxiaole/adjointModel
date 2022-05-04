%% Nov. 10, 2021
% process the MEIC data for adjoint model
clear
% close all
addpath('../../utils')

%% common parameters
year = 2017;
sectors = {'residential'
    'industry' 
    'agriculture'
    'transportation'
    'power'};

if(size(sectors,1)~=1)
    sectors = sectors';
end
pollutant = 'PM25';

load constants
r = earthR; % 6371229  meter

% configurations for Polair3d
configFile = '../../../forwardRun/config/caseConfig.cfg';

universalRawDataFolder = '../../../../universalRawData/';
% mix emission inventory
mixFolder = fullfile(universalRawDataFolder, 'emission/MIX/PM2.5/');
mixFileName = 'MICS_Asia_PM2.5_2010_0.25x0.25.nc';

% MEIC emission inventory
meicFolder = fullfile(universalRawDataFolder, ['emission/MEIC_' num2str(year)]);

% % BTH emission inventory
% bthFolder = ['../../forwardRun/rawData/emission/BTH_' num2str(year)];

% emission data folder
output = '../../../forwardRun/data/emissions';
currentMonth = 0;

%% get configurations
[configs] = readConfigs(configFile);

x_min = configs.domain.x_min;
Delta_x = configs.domain.Delta_x;
Nx = configs.domain.Nx;

y_min = configs.domain.y_min;
Delta_y = configs.domain.Delta_y;
Ny = configs.domain.Ny;

xgrid = (x_min+Delta_x/2):Delta_x:(x_min+Nx*Delta_x);
ygrid = (y_min+Delta_x/2):Delta_y:(y_min+Ny*Delta_y);

negIndex = xgrid<0;
xgrid(negIndex) = xgrid(negIndex)+360;
xgrid= max(xgrid, 0.05+0.01);

[xgrid, ygrid] = meshgrid(xgrid, ygrid);

% initial time
inputFormat = 'yyyy-MM-dd_HH-mm-ss';
initialTime = datetime(configs.domain.Date_min,'InputFormat',inputFormat);

% % negative sign is needed because timezone covert local time to utc, but we
% % want utc to local
% timezones = -timezone(xgrid(:));
% timezones = reshape(timezones, size(xgrid, 1), size(xgrid, 2));

% timezone should be in days /24
timeZones =ones(size(xgrid))*8/24;

totalHours = configs.domain.Delta_t*configs.domain.Nt/3600;

%% load mask
load Asia_Mask.mat
countryId = interp2(xMask, yMask, mask, xgrid, ygrid, 'nearest');

%% Read MEIC estimate monthly factors based on MEIC
for sector = sectors
    data = zeros(size(xgrid,1), size(xgrid,2),12);
    for monthId = 1:12
        filenameMEIC = fullfile(meicFolder, [num2str(year) '_' num2str(monthId) '_' sector{1} '_' pollutant '.nc']);
        xrangeMEIC = ncread(filenameMEIC, 'x_range');
        yrangeMEIC = ncread(filenameMEIC, 'y_range');
        dims = ncread(filenameMEIC, 'dimension');
        spacings = ncread(filenameMEIC, 'spacing');
        
        xMEIC = (xrangeMEIC(1)+spacings(1)/2):spacings(1):(xrangeMEIC(2));
        yMEIC = (yrangeMEIC(1)+spacings(2)/2):spacings(2):(yrangeMEIC(2));
        
        if(length(xMEIC)~=dims(1) || length(yMEIC)~=dims(2))
            error('MEIC dimension error')
        end
        
        [xMeic, yMeic] = meshgrid(xMEIC, yMEIC);
        areas = (cos(yMeic/180*pi)*r.*spacings(1)/180*pi)*(spacings(2)/180*pi.*r);
        
        dataTmp = ncread(filenameMEIC, 'z');
        dataTmp = reshape(dataTmp, dims(1), dims(2));
        dataTmp = flipud(dataTmp');
        
        % convert from ton/year to microg/m2/s
        daysN = eomday(year,monthId);
        dataTmp = dataTmp./areas*10^3*10^3*10^6/(daysN*24*3600);
        
        % interpolate MEIC onto the polair3d grid
        data(:,:, monthId) = interp2(xMeic, yMeic, dataTmp, xgrid, ygrid);
    end
    meic.(sector{1}).emission = data;
    meic.(sector{1}).monthlyFactors = data./mean(data,3);
end

%% read MIX
filenameMix = fullfile(mixFolder, mixFileName );
latMix = ncread(filenameMix, 'lat');
lonMix = ncread(filenameMix, 'lon');
cellsize = [mean(diff(lonMix)) mean(diff(latMix))];
[xMix, yMix] = meshgrid(lonMix, latMix);
areasMIX = (cos(yMix/180*pi)*r.*cellsize(1)/180*pi)*(cellsize(2)/180*pi.*r);

for sector = sectors
    sectorName = convertSectorName(sector{1},'MIX');
    if(~isempty(sectorName))
        data = zeros(size(xgrid,1), size(xgrid,2),12);
        datatmp = ncread(filenameMix, ['PM2.5_' sectorName]);
        
        for monthId = 1:12
            datatmpMonth = squeeze(datatmp(:,:,monthId));
            daysN = eomday(year,monthId);
            datatmpMonth = datatmpMonth'./areasMIX*10^3*10^3*10^6/(daysN*24*3600);
            data(:,:,monthId) = interp2(xMix, yMix, datatmpMonth, xgrid, ygrid);
        end
        mix.(sector{1}).emission = data;
    else
        mix.(sector{1}).emission = ones(size(xgrid,1),size(xgrid,2),12)*-1;
    end
end

% %% read BTH
% for sectorName = sectors
%     fileName = ['BTH_' num2str(year) '_' sectorName{1} '_' lower(pollutant) '_0.1x0.1.asc'];
%     fileName = fullfile(bthFolder, sectorName{1}, fileName);
%     [rawData, configsTmp, configs2] = importdata(fileName,' ',6);
%     
%     dataTmp = flipud(rawData.data);
%     for infotmp = rawData.textdata'
%         tmp = regexp(infotmp{1}, ' ', 'split');
%         info.(tmp{1}) = str2double(tmp{2});
%     end
%     
%     lonBTH = (info.xllcorner+info.cellsize/2):info.cellsize:(info.xllcorner+(info.ncols-0.5)*info.cellsize);
%     latBTH = (info.yllcorner+info.cellsize/2):info.cellsize:(info.yllcorner+(info.nrows-0.5)*info.cellsize);   
%     [xBTH, yBTH] = meshgrid(lonBTH, latBTH);
%     
%     areasBTH = (cos(yBTH/180*pi)*r.*info.cellsize/180*pi)*(info.cellsize/180*pi.*r);
%     
%     % convert from ton/year to microg/m2/s
%     dataTmp = dataTmp./areasBTH*10^3*10^3*10^6/(365*24*3600);
%     
%     data = zeros(size(xgrid,1), size(xgrid,2),12);
%     monthlyFactors = meic.(sectorName{1}).monthlyFactors;
%     for monthId = 1:12
%         factors = squeeze(monthlyFactors(:,:,monthId));
%         data(:,:,monthId) = interp2(xBTH, yBTH, dataTmp, xgrid, ygrid).*factors;
%     end
%     bth.(sectorName{1}).emission = data;
% end

%% generate the temporal emission data
for sector = sectors
    % first dimension is x, second y, and third time
    sectorName = convertSectorName(sector{1},'EDGAR');
    emissionData = zeros(size(xgrid,2),size(xgrid,1), totalHours);
    for hourId = 0:(totalHours-1)
        currentTime = initialTime + hourId/24 + timeZones;
        disp(currentTime(1))
        
        % get hourly profiles
        
        temporalFactors = getTemporalProfile(currentTime, sectorName, countryId);
        temporalFactors = temporalFactors';
        
        % get monthId
        monthId = month(currentTime(1));
   
        % combine the emissions: Mix -> MEIC -> BTH
        emissionData(:,:, hourId+1) = fillInData(emissionData(:,:, hourId+1),mix.(sector{1}).emission(:,:,monthId)');
        emissionData(:,:, hourId+1) = fillInData(emissionData(:,:, hourId+1),meic.(sector{1}).emission(:,:,monthId)');
%         emissionData(:,:, hourId+1) = fillInData(emissionData(:,:, hourId+1),bth.(sector{1}).emission(:,:,monthId)');
        
        emissionData(:,:, hourId+1) = emissionData(:,:, hourId+1).*temporalFactors;
    end
    %%
    filenameF = fullfile(output,[sectorName '.bin']);
    fileID = fopen(filenameF, 'w');
    fwrite(fileID, emissionData(:),'float');
    fclose(fileID);
end

