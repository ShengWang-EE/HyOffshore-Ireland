function [windMatrix,numDays] = loadWindDataFromNC(fileNameList)
% this function can read and organize the wind speed and other climate data from .nc file
%% get the file name according to date (year-month-day)
deleteLine = [];
for i = 1:length(fileNameList)
    fileName = fileNameList(i).name;
    if ~contains(fileName,'nc')
        deleteLine = [deleteLine,i];
    end
end
fileNameList(deleteLine) = [];
%% get data
for i = 1:length(fileNameList)
    fileName = fileNameList(i).name;
    fileInfo = ncinfo(fileName);
    % get date
    fileNameLength = length(fileName);
    for j = 1:fileNameLength-1
        if fileName(j) == '2' && fileName(j+1) == '0'
            year1 = str2num(fileName(j:j+3));
            month1 = str2num(fileName(j+4:j+5));
            day1 = str2num(fileName(j+6:j+7));
            break
        end
    end
    % start from 2012.1.1
    iday = daysact(datetime(2012,1,1), datetime(year1,month1,day1));
    numDays = daysact(datetime(2012,1,1), datetime(2022,12,31));
    islot.begin = iday*24+1; islot.end = (iday+1)*24;
    lon = ncread(fileName, 'lon'); lat = ncread(fileName, 'lat');
    ULML = ncread(fileName, 'ULML'); VLML = ncread(fileName, 'VLML');
    if i == 1
        windVelocity = cell(length(lon), length(lat));
        for ilon = 1:length(lon)
            for ilat = 1:length(lat)
                windVelocity{ilon,ilat}.ULML = zeros(numDays*24,1);
                windVelocity{ilon,ilat}.VLML = zeros(numDays*24,1);
            end
        end
    end
    counter = 0;
    for ilon = 1:length(lon)
        for ilat = 1:length(lat)
            counter = counter + 1;
%             windVelocity{ilon,ilat}.lon = lon(ilon);
%             windVelocity{ilon,ilat}.lat = lat(ilat);
%             windVelocity{ilon,ilat}.ULML(islot.begin:islot.end,1) = ULML(ilon,ilat,:);
%             windVelocity{ilon,ilat}.VLML(islot.begin:islot.end,1) = VLML(ilon,ilat,:);
            % form into a matrix
            windMatrix.lon(1,counter) = lon(ilon);
            windMatrix.lat(1,counter) = lat(ilat);
            windMatrix.ULML(islot.begin:islot.end,counter) = ULML(ilon,ilat,:);
            windMatrix.VLML(islot.begin:islot.end,counter) = VLML(ilon,ilat,:);
        end
    end
end
% add default value if some are missing
while find(windMatrix.ULML==0)
    for i = 1:numDays * 24
        if windMatrix.ULML(i,1) == 0
            windMatrix.ULML(i,:) = windMatrix.ULML(i-24,:);
            windMatrix.VLML(i,:) = windMatrix.VLML(i-24,:);
        end
    end
end
end

