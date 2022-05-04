function runAdjointModel(firstTime,hourId,outputFormat,newConfigs, fileList)
% Xiaole Zhang, Nov.17, 2021
%   run the polair3d model
currentTime = firstTime+hourId/24;
timeStr = datestr(currentTime, outputFormat);

newConfigs.domain.Date_min = timeStr;
newConfigs.source.Date_beg = timeStr;

timeEndStr = datestr(currentTime+1/24, outputFormat);
newConfigs.source.Date_end = timeEndStr;

newConfigs.save.Averaged = 'no';

for fileTmp = fileList
    genAdjointConfig(fileTmp.name, newConfigs,timeStr)
end

%% run the program
disp(fullfile('../../../adjointRun',timeStr))
cd(fullfile('../../../adjointRun',timeStr))
tic
% system('source ~/.bashrc')
system('export PATH="../../bin:$PATH";  polair3dAdjoint config/caseConfig.cfg >& log')
tc = toc;
disp(tc)

% return to the working directory
cd('../../src/preprocessing/adjoint')
end

