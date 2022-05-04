function genAdjointConfig(fileName, newConfigs, folderName)
%Generate the config files for the adjoint runs
%   Xiaole Zhang, Nov. 17, 2021
%%
forwardFolders = '../../../forwardRun';
adjointFolders = fullfile('../../../adjointRun',folderName);

newConfigPath = fullfile(adjointFolders,  'config');
if ~exist(newConfigPath, 'dir')
    mkdir(newConfigPath)
end

if ~exist(fullfile(adjointFolders, 'results'), 'dir')
    mkdir(fullfile(adjointFolders, 'results'))
end

fod = fopen(fullfile(adjointFolders, 'config', fileName), 'w');
fid = fopen(fullfile(forwardFolders,'config', fileName) , 'r');
% disp(fullfile(forwardFolders,'config', fileName))

domainsToBeModified = strtrim(fields(newConfigs));
domainName = 'default';
domainFlag = 0;

while ~feof(fid)
    tline = fgetl(fid);
    tmp = regexp(tline, '.*\[(?<domain>\w*)\].*', 'names');
    if(~isempty(tmp))
        domainName = tmp.domain;
        domainFlag = 1;
    else
        domainFlag = 0;
    end
    
    % exclude comments
    tlineSplit = regexp(tline, '#', 'split');
    if(isempty(strtrim(tlineSplit{1}))||~ismember(domainName, domainsToBeModified))
        fprintf(fod, [tline '\n']);
        continue
    end
  
    tlineParams = strtrim(tlineSplit{1});
    
    [items, posStart, posEnd, posComb, delimiters]  = regexp(tlineParams, '(\s)*\t(\s)*|(\s)*=(\s)*|(\s)*:(\s)*', 'split');
    if(length(items)>2)
        [items, posStart, posEnd, posComb, delimiters] = regexp(tlineParams, '(\s)*\t(\s)*|(\s)*=(\s)*|(\s)*:(\s)*|\s+', 'split');
    end
    
    newLine = [];
    if(length(items)>1)
        paramsToBeModified = strtrim(fields(newConfigs.(domainName)));
        
        if(mod(length(items),2)~=0)
            error(['Unpaired parameters: ' tlineParams])
        end
        %         disp(items)
        for i=1:2:length(items)
            newLine = [newLine items{i} delimiters{i}];
            
            paramNow = strtrim(items{i});
            flag = ismember(paramNow, paramsToBeModified);
            if(flag)
                newLine = [newLine num2str(newConfigs.(domainName).(paramNow))];
            else
                newLine = [newLine items{i+1}];
            end
            
            if(i+1~=length(items))
                newLine = [newLine delimiters{i+1}];
            end
        end
        
        for commentId = 2:length(tlineSplit)
            newLine = [newLine ' #' tlineSplit{commentId}];
        end
        
        fprintf(fod, [newLine '\n']);
    else
        if(strcmp(domainName, 'species')&&~domainFlag)
            fprintf(fod, [newConfigs.(domainName) '\n']);
        else
            fprintf(fod, [tline '\n']);
        end
    end
end

fclose(fid);
fclose(fod);
end

