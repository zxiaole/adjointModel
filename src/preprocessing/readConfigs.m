function [configs] = readConfigs(config)
% 2021.11.14, Xiaole Zhang
%   read configurations of the calculation
fid = fopen(config, 'r');

while ~feof(fid)
    tline = strtrim(fgetl(fid));
    
    % exclude comments
    tline = regexp(tline, '#', 'split');
    if(isempty(tline{1}))
        continue
    end
    tline = tline{1};
    
    items = regexp(tline, '(\s)*\t(\s)*|(\s)*=(\s)*|(\s)*:(\s)*|\s+', 'split');
    
    if(length(items)>1)
        if(mod(length(items),2)~=0)
            error(['Unpaired parameters: ' tline])
        end
%         disp(items)
        for i=1:2:length(items)
            val = str2double(items{i+1});
            if(isnan(val))
                configs.(items{i})=items{i+1};
            else
                configs.(items{i})=val;
            end
        end
    end
end

fclose(fid);

end

