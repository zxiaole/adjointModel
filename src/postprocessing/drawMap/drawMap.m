%% Draw the conc map
% Nov. 22, 2021 Xiaole Zhang
%%
clear

load colormap_edited.mat
addpath('../../utils')
filenames = {
    'TRO.bin';'IND.bin'; 'RCO.bin'};
fontsize = 18;

%% get basic configurations
basicConfigs

%% model configs
[configs] = readConfigs(configFile);

x_min =  configs.domain.x_min;
y_min = configs.domain.y_min;
Delta_x = configs.domain.Delta_x;
Delta_y = configs.domain.Delta_y;
sizes = [configs.domain.Nx configs.domain.Ny configs.domain.Nz configs.domain.Nt/(3600/configs.domain.Delta_t)];

Nx = sizes(1);
Ny = sizes(2);
Nt = sizes(4)-2;
dataN = Nx*Ny*Nt;

xcoord = x_min:Delta_x:x_min+(Nx-1)*Delta_x;
ycoord = y_min:Delta_y:y_min+(Ny-1)*Delta_y;

inputFormat = 'yyyy-MM-dd_HH-mm-SS';
forwardFirstTime = datetime(configs.domain.Date_min,'InputFormat',inputFormat,'TimeZone','UTC');
forwardLastTime =  forwardFirstTime + (configs.domain.Nt*configs.domain.Delta_t/3600 - 1)/24;

[sourceConfigs] = readConfigs(emissionConfigFile);

%% get measurements
siteName = 'Xian';
measSubFolder = 'rawData/measurements';
measFile = fullfile(forwardFolders, measSubFolder, [siteName '.csv']);
measTable = readtable(measFile,'Format','%{MM/dd/uuuu HH:mm}D %f %f %f %f %f %f');

measTable.Date = datetime(measTable.Date, 'TimeZone','Asia/Shanghai');

% get the period for the measurements
measTimeEnd = measTable.Date(end);
hourEnd = hours(measTimeEnd - forwardFirstTime);
measTimeStart = measTable.Date(1);
hourStart = hours(measTimeStart-forwardFirstTime);

%% read conc
dataConc = 0;
for si=1:length(filenames )
    close all
    figure('Position', [100 100 800 800])
    
    filename = fullfile(forwardFolders, 'results', filenames{si});
    
    
    fileID = fopen(filename);
    A = fread(fileID,dataN,'float');
    fclose(fileID);
    surfConc = reshape(A, Nx, Ny, Nt);
    
    surfConcMean = squeeze(mean(surfConc(:,:,hourStart:hourEnd),3));
    
    imagesc(xcoord, ycoord, surfConcMean');
    colormap(cmap)
    axis xy
    hold on
    
    %% draw city
    % detailed regions
    % 610000: Shanxi;
    %
    cityShapeFile = fullfile(forwardFolders, '/rawData/map/市WGS84.shp');
    consolidateProvince = [61];
    
    roi = shaperead(cityShapeFile);
    
    provinceCodes = {roi.provinceCo};
    provinceNames = {roi.provinceNa};
    cityCodes = {roi.cityCode};
    cityNames = {roi.cityName};
    
    namesCH = {roi.NAME};
    for i=1:length(namesCH)
        tmp = regexp(namesCH{i}, '\|','split');
        namesCH{i}=tmp{1};
    end
    
    cityCodeMap = zeros( Ny,Nx);
    
    codeList = [];
    nameList = [];
    
    for regionId = 1:length(roi)
        provCode = provinceCodes{regionId};
        provName = provinceNames{regionId};
        cityCode = cityCodes{regionId};
        cityName = cityNames{regionId};
        
        if(sum(consolidateProvince==(provCode/10000)))
            rx = roi(regionId).X;
            ry = roi(regionId).Y;
            
            blockNum = sum(isnan(rx));
            blockId = find(isnan(ry));
            blockId = [0 blockId];
            for i =1:blockNum
                % convert to image coordinates
                ix = rx((blockId(i)+1):(blockId(i+1)-1));
                iy = ry((blockId(i)+1):(blockId(i+1)-1));
                plot(ix, iy, 'k--', 'color', ones(1,3)*0.5, 'linewidth', 0.5)
            end
        end
    end
    
    %% map configs
    cityShapeFile = fullfile(forwardFolders, '/rawData/map/省WGS84.shp');
    consolidateProvince = [61];
    roi = shaperead(cityShapeFile);
    provinceCodes = {roi.provinceCo};
    for regionId = 1:length(roi)
        provCode = provinceCodes{regionId};
        rx = roi(regionId).X;
        ry = roi(regionId).Y;
        if(sum(consolidateProvince~=(provCode/10000)))
            blockNum = sum(isnan(rx));
            blockId = find(isnan(ry));
            blockId = [0 blockId];
            for i =1:blockNum
                % convert to image coordinates
                ix = rx((blockId(i)+1):(blockId(i+1)-1));
                iy = ry((blockId(i)+1):(blockId(i+1)-1));
                plot(ix, iy,'k-')
            end
        end
        
    end
    
    for regionId = 1:length(roi)
        provCode = provinceCodes{regionId};
        rx = roi(regionId).X;
        ry = roi(regionId).Y;
        if(sum(consolidateProvince==(provCode/10000)))
            linestyle =  'm-';
            blockNum = sum(isnan(rx));
            blockId = find(isnan(ry));
            blockId = [0 blockId];
            for i =1:blockNum
                % convert to image coordinates
                ix = rx((blockId(i)+1):(blockId(i+1)-1));
                iy = ry((blockId(i)+1):(blockId(i+1)-1));
                plot(ix, iy,linestyle, 'linewidth', 2)
            end
        end
    end
    axis equal
    axis tight
    
    h1 = scatter(sourceConfigs.source.Abscissa, sourceConfigs.source.Ordinate, 'ko');
    set(gca, 'xlim', [xcoord(1) xcoord(end)], 'ylim', [ycoord(1) ycoord(end)]);
    cb=cbarrow2();
    
    cb.FontSize = fontsize;
    cb.FontName = 'Arial';
    
    cb.Label.String = 'Concentraion (μg\cdotm^{-3})';
    %%
    
    xlabel('Longitude', 'fontsize', fontsize);
    ylabel('Latitude', 'fontsize', fontsize);
    set(gca, 'fontsize', fontsize, 'fontname', 'Arial');
    set(gcf,'PaperPositionMode','auto')
    pause(0.1)
    
    print([siteName '_' filenames{si}(1:3)],'-dpng','-r300')
    
end