%% Process weekly profile
% Dec. 13, 2020
%%
function processWeekly_Hourly()
data = importdata('weekly_profiles.csv');

countryCode = data.textdata(2:end,1);
sectorCode = data.textdata(2:end,2);
dataWeek = data.data;

uniqueCountryCode = unique(countryCode);
uniqueSectorCode = unique(sectorCode);

weeklyFactors = ones(length(uniqueSectorCode), ...
    length(uniqueCountryCode), 7)/7;

for i = 1:length(uniqueCountryCode)
    countryId = uniqueCountryCode{i};
    for j = 1:length(uniqueSectorCode)
        sectorId = uniqueSectorCode{j};
        id = strcmp(countryId, countryCode) & strcmp(sectorId, sectorCode);
        if(sum(id)>0)
            dataTmp = dataWeek(id,:);
            for dayId =1:7
                idTmp = dataTmp(:,1) == dayId;
                if(sum(idTmp)~=0)
                    weeklyFactors(j, i, dayId) = mean(dataTmp(idTmp,2));
                else
                    disp(countryId)
                    disp(sectorId)
                end
            end
        end
        
        if(abs(sum(weeklyFactors(j, i, :))-1)>0.0001)
            disp(sum(weeklyFactors(j, i, :)))
        end
    end
end
weeklyFactors = weeklyFactors*7;% /(1/7)

if(sum(isnan(weeklyFactors(:)))>0)
    error('nan')
end

%% hourly
data = importdata('hourly_profiles.csv');

countryCode = data.textdata(2:end,1);
sectorCode = data.textdata(2:end,2);
dataHour = data.data;


hourlyFactors =  ones(length(uniqueSectorCode), ...
    length(uniqueCountryCode), 12, 3, 24)/24;
for i = 1:length(uniqueCountryCode)
    disp(i)
    countryId = uniqueCountryCode{i};
    for j = 1:length(uniqueSectorCode)
        sectorId = uniqueSectorCode{j};
        for k = 1:12
            for h = 1:3
                id = strcmp(countryId, countryCode) ...
                    & strcmp(sectorId, sectorCode) ...
                    & dataHour(:,1) == k ...
                    & dataHour(:,2) == h;
                if(sum(id)>0)
                    dataTmp = mean(dataHour(id,3:end),1);
                    hourlyFactors(j, i, k, h, :) = dataTmp/sum(dataTmp);
                end
                
                if(abs(sum(hourlyFactors(j, i, k, h, :) )-1)>0.0001)
                    disp(sum(hourlyFactors(j, i, k, h, :) ))
                end
            end
        end        
    end
end

hourlyFactors=hourlyFactors*24;
if(sum(isnan(hourlyFactors(:)))>0)
    error('nan')
end
save weekdayProfiles.mat uniqueCountryCode uniqueSectorCode weeklyFactors hourlyFactors

end