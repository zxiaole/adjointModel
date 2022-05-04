%% Reverse the meteo data in time for the adjoint model
% revers the direction of wind vector
% Nov. 17, Xiaole Zhang
clear
addpath('../../utils')

%% configuration
forwardFolders = '../../../forwardRun';
adjointFolders = '../../../adjointRun';

genMeteoFlag = 0;

dataSubFolder = 'data';
configSubFolder = 'config';
resultSubFolder = 'results';

subfolders = {'meteo'; 'dep'; 'meteo/Kz_TM'};


filenameSave = fullfile(adjointFolders, dataSubFolder);
if ~exist(filenameSave, 'dir')
    mkdir(filenameSave)
end
if ~exist(fullfile(adjointFolders, resultSubFolder), 'dir')
    mkdir(fullfile(adjointFolders, resultSubFolder))
end

copyfile(fullfile(forwardFolders, dataSubFolder, 'ground'), fullfile(filenameSave, 'ground'))

% configurations for Polair3d
configFile = fullfile(forwardFolders, configSubFolder, 'caseConfig.cfg');
emissionConfigFile = fullfile(forwardFolders, configSubFolder, 'source.dat');

%% get configurations
[configs] = readConfigs(configFile);

sizes = [configs.domain.Nx configs.domain.Ny configs.domain.Nz configs.domain.Nt/(3600/configs.domain.Delta_t)];
typeSize = 4; % float type is 4 bytes
x_min = configs.domain.x_min;
Delta_x = configs.domain.Delta_x;
Nx = sizes(1);
y_min = configs.domain.y_min;
Delta_y = configs.domain.Delta_y;
Ny = sizes(2);

xConc = linspace(x_min, x_min+(Nx-1)*Delta_x, Nx);
yConc = linspace(y_min, y_min+(Ny-1)*Delta_y, Ny);

[emissions] = readConfigs(emissionConfigFile);
load constants
r = earthR;
gridArea = (cos(emissions.source.Ordinate/180*pi)*r.*Delta_x/180*pi)*(Delta_y/180*pi.*r);

%
inputFormat = 'yyyy-MM-dd_HH-mm-SS';

forwardFirstTime = datetime(configs.domain.Date_min,'InputFormat',inputFormat,'TimeZone','UTC');
forwardLastTime =  forwardFirstTime + (configs.domain.Nt*configs.domain.Delta_t/3600 - 1)/24;

%% generate the meteo data for adjoint run
if(genMeteoFlag)
    for i=1:length(subfolders)
        subfolder = subfolders{i};
        fileFolder = fullfile(forwardFolders,dataSubFolder, subfolder);
        filenames = dir([fileFolder, '/*.bin']);
        
        for j=1:length(filenames)
            filename = filenames(j).name;
            disp(filename)
            
            filesize = filenames(j).bytes;
            sizeTmp = sizes;
            
            layerN = floor(filesize/sizes(1)/sizes(2)/sizes(4)/typeSize);
            sizeTmp(3) = layerN;
            
            residual = filesize/sizes(1)/sizes(2)/sizes(4)/typeSize/layerN;
            residual = residual - floor(residual);
            
            
            if(residual < eps)
                factorReverse = 1;
            elseif(abs(residual-1/sizes(1))<eps)
                factorReverse = -1; % for wind vector
                sizeTmp(1) = sizeTmp(1) + 1;
            elseif(abs(residual-1/sizes(2))<eps)
                factorReverse = -1; % for wind vector
                sizeTmp(2) = sizeTmp(2) + 1;
            else
                error('size error');
            end
            dataN = prod(sizeTmp);
            
            fileID = fopen(fullfile(fileFolder, filename));
            A = fread(fileID,dataN,'float');
            fclose(fileID);
            datatmp = reshape(A, sizeTmp(1), sizeTmp(2),sizeTmp(3),sizeTmp(4));
            datatmp = datatmp(:,:,:,end:-1:1);
            
            if(factorReverse==-1)
                datatmp = datatmp*factorReverse;
            end
            
            filenameSaveFull = fullfile(filenameSave, subfolder);
            if ~exist(filenameSaveFull, 'dir')
                mkdir(filenameSaveFull)
            end
            
            fileSaveId = fopen(fullfile(filenameSaveFull, filename), 'w');
            fwrite(fileID, datatmp(:),'float');
            fclose(fileSaveId);
        end
    end
end
%% new configs
species = 'XA_PM25';
startDate = configs.domain.Date_min;%'2021-02-01_00-00-00';

newConfigs.domain.Date_min = startDate;
newConfigs.domain.Nt = 2016; % 7 days
newConfigs.options.With_point_emission = 'yes';
newConfigs.options.With_surface_emission = 'no';
newConfigs.surface_emission.Fields = species;
newConfigs.deposition.Fields = species;
newConfigs.scavenging.Fields = species;

newConfigs.source.Species = species;
newConfigs.source.Date_beg = startDate;
newConfigs.source.Date_end = '2021-02-01_03-00-00';
newConfigs.source.Rate = gridArea;

newConfigs.save.Species = species;
newConfigs.species = species;

newConfigs.meteo.Data_dir = '../data';
newConfigs.meteo.Data_path = '../data/meteo';
newConfigs.meteo.Data_dir_emission = '../data';

configFilesToBeCopied = dir(fullfile(forwardFolders, configSubFolder, '/*.cfg'));
dataFilesToBeCopied = dir(fullfile(forwardFolders, configSubFolder, '/*.dat'));
fileList = [configFilesToBeCopied;dataFilesToBeCopied]';

%
inputFormat = 'yyyy-MM-dd_HH-mm-SS';
outputFormat = 'yyyy-mm-dd_HH-MM-SS';

firstTime = datetime(startDate,'InputFormat',inputFormat) ;

%% get measurements
siteName = 'Xian';
measSubFolder = 'rawData/measurements';
measFile = fullfile(forwardFolders, measSubFolder, [siteName '.csv']);
measTable = readtable(measFile,'Format','%{MM/dd/uuuu HH:mm}D %f %f %f %f %f %f');

measTable.Date = datetime(measTable.Date, 'TimeZone','Asia/Shanghai');

%% get the period for the measurements
measTimeEnd = measTable.Date(end);
hourStart = hours(forwardLastTime-measTimeEnd);
measTimeStart = measTable.Date(1);
hourEnd = hours(forwardLastTime-measTimeStart);

parfor hourId=hourStart:hourEnd
    runAdjointModel(firstTime,hourId,outputFormat,newConfigs, fileList)
end

