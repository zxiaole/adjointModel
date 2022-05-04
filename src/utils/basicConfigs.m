
%% basic configurations
forwardFolders = '../../../forwardRun';
adjointFolders = '../../../adjointRun';

dataSubFolder = 'data';
configSubFolder = 'config';
resultSubFolder = 'results';
measSubFolder = 'rawData/measurements';

subfolders = {'meteo'; 'dep'; 'meteo/Kz_TM'};

filenameSave = fullfile(adjointFolders, dataSubFolder);
% configurations for Polair3d
configFile = fullfile(forwardFolders, configSubFolder, 'caseConfig.cfg');
emissionConfigFile = fullfile(forwardFolders, configSubFolder, 'source.dat');
forwardConfigFile = fullfile(forwardFolders, configSubFolder, 'case-saver.cfg');