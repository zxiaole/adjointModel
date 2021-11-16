
%%
clear
regionName = 'Global';%'Austrilia';%'SouthAmerica';%'Asia'; %'Europe';

filenameSave = ['./' regionName '/'] ;
mkdir(filenameSave)

setDomain
load([regionName '_Mask.mat'])
load('weekHourProfiles.mat')

x_min = domain.(regionName).xMin;
Delta_x = domain.(regionName).xDelta;
Nx = domain.(regionName).Nx;

y_min = domain.(regionName).yMin;
Delta_y = domain.(regionName).yDelta;
Ny = domain.(regionName).Ny;

[idx, idy] = meshgrid(1:Nx, 1:Ny);

xgrid = x_min:Delta_x:(x_min+(Nx-1)*Delta_x);
ygrid = y_min:Delta_y:(y_min+(Ny-1)*Delta_y);

timeZones = round(xgrid/360*24);
timeZones = repmat(timeZones, length(ygrid),1)/24; % in days

negIndex = xgrid<0;
xgrid(negIndex) = xgrid(negIndex)+360;
xgrid= max(xgrid, 0.05+0.01);

[xgrid, ygrid] = meshgrid(xgrid, ygrid);

%%
inputFormat = 'yyyyMMdd HH:mm';
firstDay = '20150101 00:00';
initialTime = datetime(firstDay,'InputFormat',inputFormat);

%%
sectorName = 'RCO';
hourlyProfiles = importdata("hourly_profiles.csv");
hourlyProfilesData = hourlyProfiles.data;
hourlyProfilesCode = hourlyProfiles.textdata(2:end, 1:2);

factorMap = ones(900, 1800);

%%

f = filesep;

dataFolder = '../data/edgar_v50';

edgarCat.PM25.IND = {'IND','CHE', 'FOO_PAP','IRO', 'NFE', ...
    'NMM', 'PRU_SOL', 'REF_TRF', 'PRO'};

edgarCat.PM25.RCO = {'RCO', 'SWD_INC', 'SWD_LDF'};
edgarCat.PM25.TRO = {'TRO_noRES', 'TRO_RES', 'TNR_Other','TNR_Ship','TNR_Aviation_CDS',...
    'TNR_Aviation_CRS','TNR_Aviation_LTO'};


edgarCat.PM25.ENE = {'ENE'};
edgarCat.PM25.AGS = {'AGS', 'AWB','MNM'};
years = 2012:2015;

edgarCat.BC.RCO = {'RCO', 'SWD_INC', 'SWD_LDF'};
edgarCat.BC.TRO = {'TRO_noRES', 'TRO_RES', 'TNR_Other','TNR_Ship','TNR_Aviation_CDS',...
    'TNR_Aviation_CRS','TNR_Aviation_LTO'};

edgarCat.BC.IND = {'IND','CHE', 'FOO_PAP','IRO', 'NFE', ...
    'NMM', 'PRU_SOL', 'REF_TRF', 'PRO'};
edgarCat.BC.ENE = {'ENE'};
edgarCat.BC.AGS = {'AGS', 'AWB','MNM'};

for pollutantsName = fieldnames(edgarCat)'
    disp(pollutantsName)
    for sectorName = fieldnames(edgarCat.(pollutantsName{1}))'
        sectorEmssions = zeros( Nx, Ny, 365*24+1);
        for subSectorName = edgarCat.(pollutantsName{1}).(sectorName{1})
            disp(subSectorName)
            subSectorId = find(strcmp(subSectorName{1}(1:3), uniqueSectorCode) == 1);
            
            if(strcmp(pollutantsName, 'PM25'))
                fileTemp = 'PM2.5';
            else
                fileTemp = pollutantsName{1};
            end
            
            year=years(1);
            filename = [dataFolder f pollutantsName{1} f subSectorName{1} f 'v50_' fileTemp '_' num2str(year) '_' subSectorName{1} '.0.1x0.1.zip'];
            if(~isfile(filename))
                continue;
            end
            
            emiAvg = 0;
            for year = years
                filename = [dataFolder f pollutantsName{1} f subSectorName{1} f 'v50_' fileTemp '_' num2str(year) '_' subSectorName{1} '.0.1x0.1.zip'];
                unzip(filename,[dataFolder f pollutantsName{1} f subSectorName{1}])
                filename = [filename(1:end-3) 'nc'];
                
                emiAvg = emiAvg + ncread(filename,['emi_' lower(fileTemp)]);
                
                if(year == 2015)
                    emi = ncread(filename,['emi_' lower(fileTemp)]);
                    emi = emi(:,:)';
                    lat = ncread(filename,'lat');
                    lon = ncread(filename,'lon');
                    [lonGrid, latGrid] = meshgrid(lon(:), lat(:));
                end
                
                delete(filename)
            end
            
            % calculate the average emissions between 2012 and 2015
            emiAvg = emiAvg(:,:)'/length(years);
            
            avgFactors = ones(size(emiAvg));
            
            idTmp = emiAvg>0 & emi>0;
            avgFactors(idTmp) = emiAvg(idTmp)./emi(idTmp);
            
            %% load monthly emission for 2015 and rescaled to the average emissions
            filenameMonthly = [dataFolder f pollutantsName{1} f subSectorName{1} f 'v50_' fileTemp '_' num2str(year) '_monthly_' subSectorName{1} '_nc.zip'];
            targetFolder = [dataFolder f pollutantsName{1} f subSectorName{1}];
            unzip(filenameMonthly,targetFolder);
            
            dataMonthly = zeros(12, Ny, Nx);
            for monthId = 1:12
                filenameMonthly = [targetFolder f subSectorName{1} f 'v50_' fileTemp '_2015_' num2str(monthId) '_' subSectorName{1} '.0.1x0.1.nc'];
                monthDataTmp1 = ncread(filenameMonthly,['emi_' lower(fileTemp)])'.*avgFactors;
                
                monthDataTmp = interp2(lonGrid, latGrid, monthDataTmp1, xgrid, ygrid); % interpolation to the location of the farm
                
                % convert from kg /m2/s to microg/m2/s
                monthDataTmp = monthDataTmp*10^9;
                
                dataMonthly(monthId,:,:) = monthDataTmp; 
            end
            
 
            %% generate hourly emission
            for i = 0:365*24
                disp(initialTime + i/24);
                currentTime = initialTime + i/24 + timeZones;
                
                hourId = hour(currentTime);
                hourId(hourId==0) = 24;
                weekDayIndex = weekday(currentTime);
                weekDayIndex = weekDayIndex-1;
                
                weekDayIndex(weekDayIndex==0)=7;
                
                
                % month id
                monthId = month(currentTime);
                
                dayTypeId = ones(size(weekDayIndex));
                
                dayTypeId(weekDayIndex==6)=2;
                dayTypeId(weekDayIndex==7)=3;
                
                
                %% get month matrix
                monthDataTmp = zeros(Ny, Nx); 
                idValidMonth = sub2ind(size(dataMonthly), monthId(:), idy(:), idx(:));
                monthDataTmp(:) = dataMonthly(idValidMonth);
                
                %% correct by weekday
                idtmp = monthDataTmp(:)>0;
                
                weekdayFactorTmp = ones(Ny, Nx);
                idValidWeek = sub2ind(size(weeklyFactors), ones(sum(idtmp),1)*subSectorId, ...
                    mask(idtmp), weekDayIndex(idtmp));
                weekdayFactorTmp(idtmp) = weeklyFactors(idValidWeek);   
                monthDataTmp = monthDataTmp.*weekdayFactorTmp;
                
                %% correct by hour
                hourFactorTmp = ones(Ny, Nx);
                idValidHour = sub2ind(size(hourlyFactors),...
                    ones(sum(idtmp),1)*subSectorId, ...
                    mask(idtmp), ...
                    monthId(idtmp), ...
                    dayTypeId(idtmp), ...
                    hourId(idtmp));
                hourFactorTmp(idtmp) = hourlyFactors(idValidHour);
                monthDataTmp = monthDataTmp.*hourFactorTmp;
                
                sectorEmssions(:, :, i+1) = sectorEmssions(:, :, i+1) + monthDataTmp';
            end
        end
        
        
        filenameF = [filenameSave pollutantsName{1} '_' sectorName{1} '.bin'];
        fileID = fopen(filenameF, 'w');
        fwrite(fileID, sectorEmssions(:),'float');
        fclose(fileID);
    end
    
    if(strcmp(pollutantsName{1}, 'BC'))
        tmp = 0;
        for sectorName = fieldnames(edgarCat.(pollutantsName{1}))'
            filenameF = [filenameSave pollutantsName{1} '_' sectorName{1} '.bin'];
            fileID = fopen(filenameF, 'r');
            tmp = tmp + fread(fileID,'float');
            fclose(fileID);
        end
        filenameF = [filenameSave pollutantsName{1} '_total.bin'];
        fileID = fopen(filenameF, 'w');
        fwrite(fileID, tmp,'float');
        fclose(fileID);
    end
end

%%

