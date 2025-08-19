function [selectedWindData] = windDataSample(windMatrix,day1,numDays,dayBand)
% this function is to get the sample of data
% get the near day data
iday1.begin = day1 - dayBand; iday1.end = day1 + dayBand;
if iday1.begin <= 0
    selectedWindData.ULML = [windMatrix.ULML((numDays+iday1.begin)*24:numDays*24,:); windMatrix.ULML(1:iday1.end*24,:)];
    selectedWindData.VLML = [windMatrix.VLML((numDays+iday1.begin)*24:numDays*24,:); windMatrix.VLML(1:iday1.end*24,:)];
elseif iday1.end > numDays
    selectedWindData.ULML = [windMatrix.ULML((iday1.begin)*24:numDays*24,:); windMatrix.ULML(1:(iday1.end-numDays)*24,:)];   
    selectedWindData.VLML = [windMatrix.VLML((iday1.begin)*24:numDays*24,:); windMatrix.VLML(1:(iday1.end-numDays)*24,:)];  
else
    selectedWindData.ULML = [windMatrix.ULML(iday1.begin*24:iday1.end*24,:)];   
    selectedWindData.VLML = [windMatrix.VLML(iday1.begin*24:iday1.end*24,:)];   
end
selectedWindData.speed = sqrt(selectedWindData.ULML.^2+selectedWindData.VLML.^2);
end

