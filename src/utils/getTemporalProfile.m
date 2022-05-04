function temporalFactors = getTemporalProfile(currentDate, sectorName, mask)
% Nov. 14, 2021, Xiaole Zhang
%   Generate the hourly profile

load weekHourProfiles.mat

subSectorId = find(strcmp(sectorName, uniqueSectorCode) == 1);

% month id
monthId = month(currentDate);

% day id
weekDayIndex = weekday(currentDate);
weekDayIndex = weekDayIndex-1;

weekDayIndex(weekDayIndex==0)=7;

dayTypeId = ones(size(weekDayIndex));

dayTypeId(weekDayIndex==6)=2;
dayTypeId(weekDayIndex==7)=3;

hourId = hour(currentDate);
hourId(hourId==0) = 24;

%% correct by weekday
weekdayFactorTmp = ones(size(currentDate));
idValidWeek = sub2ind(size(weeklyFactors), ones(length(weekdayFactorTmp(:)),1)*subSectorId, ...
    mask(:), weekDayIndex(:));
weekdayFactorTmp = weeklyFactors(idValidWeek);

%% correct by hour
hourFactorTmp = ones(size(currentDate));
idValidHour = sub2ind(size(hourlyFactors),...
    ones(length(hourFactorTmp(:)),1)*subSectorId, ...
    mask(:), ...
    monthId(:), ...
    dayTypeId(:), ...
    hourId(:));
hourFactorTmp = hourlyFactors(idValidHour);

%%
temporalFactors = reshape(weekdayFactorTmp.*hourFactorTmp, size(currentDate,1),  size(currentDate,2));
end

