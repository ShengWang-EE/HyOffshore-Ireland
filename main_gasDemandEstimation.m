clear
clc

totalGasDemand = 24.7884; % 13th Dec 2023, for the hour 17:00 to 18:00, Mm3/day
coordinateFile = projectPath('data','inputs','geospatial','list of coodinates.xlsx');
irishPopulationFile = projectPath('data','inputs','demographics','irish population density 2020.csv');
ukPopulationFile = projectPath('data','inputs','demographics','UK population density 2020.csv');

coordinatesTable = readtable(projectPath('data','inputs','geospatial','list of coodinates.xlsx'));
coordinates = table2array(coordinatesTable(:,3:4));
irishPopulationTable = readtable(irishPopulationFile);
UKPopulationTable = readtable(ukPopulationFile);
irishPopulationDensity = table2array(irishPopulationTable);
UKPopulationDensity = table2array(UKPopulationTable);
populationDensity = [irishPopulationDensity; UKPopulationDensity];
% plot3(irishPopulationTable.X,irishPopulationTable.Y,irishPopulationTable.Z);

% 先算一下各点之间的distance
distance1 = zeros(105,105);
for i = 1:105
    for j = i:105
        distance1(i,j) = sqrt((coordinates(i,1) - coordinates(j,1))^2 ...
            + (coordinates(i,2) - coordinates(j,2))^2);
    end
end
distance1(distance1==0) = nan;
minDistance = min(min(distance1));
searchingRadius = minDistance;
for i = 1:105 % number of gas bus
    distanceToBusi = sqrt(sum((populationDensity(:,1:2) - [coordinates(i,1),coordinates(i,2)]).^2,2));
    demandFactor(i) = (distanceToBusi.^(-2) ./ sum(distanceToBusi.^(-2)))' * populationDensity(:,3).^1.75;
end
demandFactor = demandFactor./sum(demandFactor);
plot(demandFactor)
gasDemand = demandFactor' * totalGasDemand;

%%
clear
clc

coordinatesTable = readtable(coordinateFile);
coordinates = table2array(coordinatesTable(:,3:4));

connectionTable = readtable(projectPath('data','inputs','system','Irish 144 gas network data.xlsx'),'sheet',6);
connection = table2array(connectionTable(:,1:2));
curveFactor = 1.0;
for i = 1:144 % number of pipeline
    frombus(i) = connection(i,1); tobus(i) = connection(i,2); 
    pipelineLength(i,1) = deg2km( distance(coordinates(frombus(i),1), coordinates(frombus(i),2), ...
        coordinates(tobus(i),1), coordinates(tobus(i),2)) ) * curveFactor;
    
end
totalLength = sum(pipelineLength(4:end));
