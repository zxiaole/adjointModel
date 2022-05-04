%% Generate the city code map
% Nov. 22, 2021 Xiaole Zhang
%%
clear
addpath('../')
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


%% map configs
cityShapeFile = fullfile(forwardFolders, '/rawData/map/县WGS84.shp');

% detailed regions
% 610000: Shanxi;
% 
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
        cityFlag = 1;
    else
        cityFlag = 0;
    end

    rx = roi(regionId).X;
    ry = roi(regionId).Y;
    
    
    blockNum = sum(isnan(rx));
    blockId = find(isnan(ry));
    blockId = [0 blockId];
    for i =1:blockNum
        % convert to image coordinates
        ix = (rx((blockId(i)+1):(blockId(i+1)-1)) - x_min)/Delta_x + 1;
        iy = (ry((blockId(i)+1):(blockId(i+1)-1)) - y_min)/Delta_y + 1;
        maskTmp = poly2mask(ix,iy,Ny,Nx);
        
        if(cityFlag)
            cityCodeMap(maskTmp) = cityCode;
            currentCode = cityCode;
            currentName = cityName;
        else
            cityCodeMap(maskTmp) = provCode;
            currentCode = provCode;
            currentName = provName;
        end
    end
    
    if(isempty(find(codeList == currentCode, 1)))
        codeList = [codeList;currentCode];
        nameList = [nameList; {currentName}];
    end
end

codeList = [codeList; 0];
nameList = [nameList; {'Others'}];

cityCodeMap = cityCodeMap';
saveFileName = '../regionMask.mat';
save(saveFileName, 'cityCodeMap', 'nameList', 'codeList')

