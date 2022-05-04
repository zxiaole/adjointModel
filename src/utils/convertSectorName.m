function outputName = convertSectorName(sector,inventoryName)
% convert the name of sector to the corresponding name in the inventory
%  Xiaole ZHang, Nov. 14, 2021
if(strcmpi(inventoryName, 'EDGAR'))
    switch sector
        case {'industry'}
            outputName = 'IND';
        case {'agriculture'}
            outputName = 'AGS';
        case {'transportation'}
            outputName = 'TRO';
        case {'power'}
            outputName = 'ENE';
        case {'residential'}
            outputName = 'RCO';
    end
elseif(strcmpi(inventoryName, 'MIX'))
    switch sector
        case {'industry'}
            outputName = 'INDUSTRY';
        case {'agriculture'}
            outputName = '';
        case {'transportation'}
            outputName = 'TRANSPORT';
        case {'power'}
            outputName = 'POWER';
        case {'residential'}
            outputName = 'RESIDENTIAL';
    end
end
end

