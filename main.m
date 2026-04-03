clear
clc
yalmip('clear')
systemDataFile = projectPath('data','inputs','system','Irish energy system data.xlsx');
windDataDir = projectPath('wind data');
approvedAreaFile = projectPath('data','inputs','geospatial','offshore area (approved).shp');
foreshoreAreaFile = projectPath('data','inputs','geospatial','offshore areas (early planning).shp');
checkpointDir = projectPath('data','checkpoints');
if ~exist(checkpointDir, 'dir')
    mkdir(checkpointDir);
end

mpc.bus = table2array(readtable(systemDataFile,'sheet','bus','range','A2:M1197'));
mpc.gen = table2array(readtable(systemDataFile,'sheet','gen','range','A2:U140'));
mpc.branch = table2array(readtable(systemDataFile,'sheet','branch','range','A2:M1369'));
mpc.gencost = table2array(readtable(systemDataFile,'sheet','gencost','range','A2:G140'));
mpc.Gbus = table2array(readtable(systemDataFile,'sheet','Gbus','range','A2:F145'));
mpc.Gline = table2array(readtable(systemDataFile,'sheet','Gline','range','A2:I145'));
mpc.Gsou = table2array(readtable(systemDataFile,'sheet','Gsou','range','A2:D3'));
mpc.Gcost = table2array(readtable(systemDataFile,'sheet','Gcost','range','A2:A3'));
mpc.GEcon = table2array(readtable(systemDataFile,'sheet','GEcon','range','A2:D14'));
mpc.ptg = table2array(readtable(systemDataFile,'sheet','ptg','range','A2:E8'));
mpc.gentype = table2array(readtable(systemDataFile,'sheet','gentype','range','A2:A140'));
mpc.version = '2';
mpc.baseMVA = 100;

mpc0 = mpc;
% mpc.branch(:,6:8) = mpc.branch(:,6:8) * 10;
% mpc.gen(131:139,9) = mpc.gen(131:139,9) * 10;
% mpc.bus(:,3) = mpc.bus(:,3)*1.1;
[electricityResults] = runopf(mpc);
% [GEresults,GEinfo] = runopf_ge(mpc);
% 通过测试发现纯电力系统能够消纳东部海岸约5GW的风电
%% wind power
% load nc winf speed files
fileNameList = dir(windDataDir);
[windMatrix,numDays] = loadWindDataFromNC(fileNameList); % load wind data from nc file
windMatrix.speed = sqrt(windMatrix.ULML.^2+windMatrix.VLML.^2);
windMatrix.speedHist = reshape(windMatrix.speed,[size(windMatrix.speed,1)*size(windMatrix.speed,2),1]);
[windTurbine] = windTurbinePara();
totalWindCapacity1st  = 5000; %MW
nWindTurbineApproved = totalWindCapacity1st / windTurbine.ratedPower;
% get the geographycal locations of OWFs
OWFapprovedArea = shaperead(approvedAreaFile);
OWFforeshoreArea = shaperead(foreshoreAreaFile);
nApproved = size(OWFapprovedArea,1);
nForeshore = size(OWFforeshoreArea,1);
approvedAreaTotal = 0;
approvedGrossAreaTotal = 0;
for i = 1:nApproved
    OWFapprovedArea(i).X(isnan(OWFapprovedArea(i).X)) = [];
    OWFapprovedArea(i).Y(isnan(OWFapprovedArea(i).Y)) = [];
    OWFapprovedArea(i).center = [mean(OWFapprovedArea(i).X),mean(OWFapprovedArea(i).Y)];
    OWFapprovedArea(i).area = polyarea(deg2km(OWFapprovedArea(i).X),deg2km(OWFapprovedArea(i).Y))*1e6; %m2
    OWFapprovedArea(i).length = -deg2km(OWFapprovedArea(i).BoundingBox(1,1) - OWFapprovedArea(i).BoundingBox(2,1)) * 1e3;
    OWFapprovedArea(i).height = -deg2km(OWFapprovedArea(i).BoundingBox(1,2) - OWFapprovedArea(i).BoundingBox(2,2)) * 1e3;
    OWFapprovedArea(i).grossArea =  OWFapprovedArea(i).length * OWFapprovedArea(i).height;
    approvedAreaTotal = approvedAreaTotal + OWFapprovedArea(i).area;
    approvedGrossAreaTotal = approvedGrossAreaTotal + OWFapprovedArea(i).grossArea;
end
for i = 1:nForeshore
    OWFforeshoreArea(i).X(isnan(OWFforeshoreArea(i).X)) = [];
    OWFforeshoreArea(i).Y(isnan(OWFforeshoreArea(i).Y)) = [];
    OWFforeshoreArea(i).center = [mean(OWFforeshoreArea(i).X),mean(OWFforeshoreArea(i).Y)];
end
% find the nearest location for each polygon
for i = 1:nApproved
    distanceSquare = (windMatrix.lon - OWFapprovedArea(i).center(1)).^2 + (windMatrix.lat - OWFapprovedArea(i).center(2)).^2;
    [~, OWFapprovedArea(i).nearestGrid] = min(distanceSquare);
    OWFapprovedArea(i).ULML = windMatrix.ULML(:,OWFapprovedArea(i).nearestGrid);
    OWFapprovedArea(i).VLML = windMatrix.VLML(:,OWFapprovedArea(i).nearestGrid);
end
% calculate the micrositting of WTs
for i = 1:nApproved
    OWFapprovedArea(i).nWindTurbine = nWindTurbineApproved / approvedGrossAreaTotal * OWFapprovedArea(i).grossArea;
    OWFapprovedArea(i).lengthSquare = sqrt(OWFapprovedArea(i).grossArea / OWFapprovedArea(i).nWindTurbine);
    OWFapprovedArea(i).nColumn = abs(round( deg2km(OWFapprovedArea(i).BoundingBox(1,1) - OWFapprovedArea(i).BoundingBox(2,1)) ...
        * 1e3 / OWFapprovedArea(i).lengthSquare ));
    OWFapprovedArea(i).nRow = abs(round( deg2km(OWFapprovedArea(i).BoundingBox(1,2) - OWFapprovedArea(i).BoundingBox(2,2)) ...
        * 1e3 / OWFapprovedArea(i).lengthSquare ));
end
save(projectPath('data','checkpoints','stop1.mat'))
%% wake effect demo
clear
clc
load(projectPath('data','checkpoints','stop1.mat'))
% select a day
day1 = daysact(datetime(2012,1,1), datetime(2022,12,1)); % typical day 2022-12-01
dayBand = 15;
[selectedWindData] = windDataSample(windMatrix,day1,numDays,dayBand);       % 选中日期前后加减dayBand的天数，所以这里是30天
% calculate wake effect for a single wind turbine
for ix = 1:2000
    for iy = 1:200
        delta_x = 10*ix; delta_y = 10*iy-1000;
        sigma = (windTurbine.rotorRadius + 0.56 / (log(windTurbine.hubHeight) - log(windTurbine.roughness)) .* delta_x) / 2;
        windDecreaseFactorForSingle(ix,iy) = (1 - sqrt(1 - windTurbine.C_T / 2 ./ (sigma / windTurbine.rotorRadius)^2) ) * exp(-delta_y.^2 / 2 / sigma.^2);
    end
end

% demonstration for wake effect for a single time point
resolution = 10; %m
[windDecreaseFactorTotalDemo,gridCoordinate_new,dataPlot] = wakeEffectDemo(OWFapprovedArea(1),selectedWindData.ULML(1,1),selectedWindData.VLML(1,2),windTurbine,resolution);
% plot3(dataPlot(:,1),dataPlot(:,2),dataPlot(:,3));
save(projectPath('data','checkpoints','stop2.mat'))
%% wake effect calculation
clear
clc
load(projectPath('data','checkpoints','stop2.mat'))
% calculate wake effect for wind turbines
windDecreaseFactorTotal = cell(size(selectedWindData.ULML,1),nApproved);
for i = 1:nApproved
    windSiteData = [selectedWindData.ULML(:,OWFapprovedArea(i).nearestGrid), selectedWindData.VLML(:,OWFapprovedArea(i).nearestGrid)];
    xWT = 0:OWFapprovedArea(i).lengthSquare:OWFapprovedArea(i).length;
    yWT = 0:OWFapprovedArea(i).lengthSquare:OWFapprovedArea(i).height;
    % shift coordinates
    for k = 1:size(windSiteData,1)
        WTcoordinate_new = cell(OWFapprovedArea(i).nRow, OWFapprovedArea(i).nColumn);
        % ULML: surface eastward wind (originates in the east and blows in a westward direction.); VLML: surface northward wind
        theta = atan(windSiteData(k,2) / windSiteData(k,1)); % ↙
        for iRow = 1:OWFapprovedArea(i).nRow
            for iCol = 1:OWFapprovedArea(i).nColumn
                WTcoordinate_new{iRow,iCol} = [cos(theta),-sin(theta);sin(theta),cos(theta)] * [xWT(iCol); yWT(iRow)];
            end
        end
        % wake effect from each TW
        for iRow = 1:OWFapprovedArea(i).nRow
            for iCol = 1:OWFapprovedArea(i).nColumn
                windDecreaseFactor = zeros(OWFapprovedArea(i).nRow,OWFapprovedArea(i).nColumn);
                for iRow_ref = 1:OWFapprovedArea(i).nRow % source of the wake effect
                    for iCol_ref = 1:OWFapprovedArea(i).nColumn
                        delta_x = WTcoordinate_new{iRow_ref,iCol_ref}(1) - WTcoordinate_new{iRow,iCol}(1);
                        delta_y = abs(WTcoordinate_new{iRow_ref,iCol_ref}(2) - WTcoordinate_new{iRow,iCol}(2));
                        if delta_x < 0 || delta_y == 0
                            windDecreaseFactor(iRow_ref,iCol_ref) = 0;
                        else
                            sigma = (windTurbine.rotorRadius + 0.56 / (log(windTurbine.hubHeight) - log(windTurbine.roughness)) .* delta_x) / 2;
                            windDecreaseFactor(iRow_ref,iCol_ref) = (1 - sqrt(1 - windTurbine.C_T / 2 ./ (sigma / windTurbine.rotorRadius)^2) ) * exp(-delta_y.^2 / 2 / sigma.^2);
                        end
                    end
                end
                windDecreaseFactorTotal{k,i}(iRow,iCol) = sum(sum(windDecreaseFactor.^2));
            end
        end
    end
end
% calculate OWF generations with wake effect
OWFgeneration0 = zeros(size(windSiteData,1),nApproved);
for i = 1:nApproved
    windSiteData = [selectedWindData.ULML(:,OWFapprovedArea(i).nearestGrid), selectedWindData.VLML(:,OWFapprovedArea(i).nearestGrid)];
    for k = 1:size(windSiteData,1)
        windSpeed_new{i}(:,:,k) = sqrt(windSiteData(k,1).^2+windSiteData(k,2).^2) * (1-windDecreaseFactorTotal{k,i});
        WTgeneration{i}(:,:,k) = windTurbineGeneration(windSpeed_new{i}(:,:,k),windTurbine);
    end
    OWFgeneration0(:,i) = reshape(sum(sum(WTgeneration{i})),[size(windSiteData,1),1]);
end

save(projectPath('data','checkpoints','stop3.mat'))
%% select samples, calculate prediction errors
clear
clc
load(projectPath('data','checkpoints','stop3.mat'))

% get statistics of OWF generation 
OWFgeneration96 = zeros(size(OWFgeneration0,1)*4,size(OWFgeneration0,2));   % 8 wind farm, two of which are close and use the same wind speed data
for i = 1:size(OWFgeneration0,2)
    OWFgeneration96(:,i) = interp(OWFgeneration0(:,i),4); % extend to 96 points
    OWFgeneration24(:,i) = interp(OWFgeneration0(:,i),1);
end
% 96 （不用96个点，用24个点就够了）
% OWFgeneration96(OWFgeneration96<0) = 0;
% nSelectedDays = floor(size(OWFgeneration96,1) / 96);
% OWFgeneration_sum = sum(OWFgeneration96,2);
% OWFgenStatistic96.curve = zeros(96,nApproved);
% OWFgenStatistic96.curve_sum = zeros(96,1);
% for iday = 1:nSelectedDays
%     OWFgenStatistic96.curve = OWFgenStatistic96.curve + OWFgeneration96((iday-1)*96+1:iday*96,:);
%     OWFgenStatistic96.curve_sum = OWFgenStatistic96.curve_sum + OWFgeneration_sum((iday-1)*96+1:iday*96,:);
% end
% OWFgenStatistic96.curve = OWFgenStatistic96.curve / nSelectedDays;
% OWFgenStatistic96.curve_sum = OWFgenStatistic96.curve_sum / nSelectedDays;
% generationError = zeros(96,nApproved,nSelectedDays);
% generationError_sum = zeros(96,nSelectedDays);
% for iday = 1:nSelectedDays
%     generationError(:,:,iday) = OWFgeneration96((iday-1)*96+1:iday*96,:) - OWFgenStatistic96.curve;
%     generationError_sum(:,iday) = OWFgeneration_sum((iday-1)*96+1:iday*96,:) - OWFgenStatistic96.curve_sum;
% end
% OWFgenStatistic96.errorMean = mean(generationError,3);
% OWFgenStatistic96.errorStd = std(generationError,[],3);
% OWFgenStatistic96.errorRelativeStd = OWFgenStatistic96.errorStd ./ OWFgenStatistic96.curve;
% OWFgenStatistic96.errorMean_sum = mean(generationError_sum,2);
% OWFgenStatistic96.errorStd_sum = std(generationError_sum,[],2);
% OWFgenStatistic96.errorRelativeStd_sum = OWFgenStatistic96.errorStd_sum ./ OWFgenStatistic96.curve_sum;
% 24
OWFgeneration24(OWFgeneration24<0) = 0;
nSelectedDays = floor(size(OWFgeneration24,1) / 24);
OWFgeneration_sum = sum(OWFgeneration24,2);
OWFgenStatistic24.curve = zeros(24,nApproved);
OWFgenStatistic24.curve_sum = zeros(24,1);
for iday = 1:nSelectedDays
    OWFgenStatistic24.curve = OWFgenStatistic24.curve + OWFgeneration24((iday-1)*24+1:iday*24,:);
    OWFgenStatistic24.curve_sum = OWFgenStatistic24.curve_sum + OWFgeneration_sum((iday-1)*24+1:iday*24,:);
end
OWFgenStatistic24.curve = OWFgenStatistic24.curve / nSelectedDays;
OWFgenStatistic24.curve_sum = OWFgenStatistic24.curve_sum / nSelectedDays;
generationError = zeros(24,nApproved,nSelectedDays);
generationError_sum = zeros(24,nSelectedDays);
for iday = 1:nSelectedDays
    generationError(:,:,iday) = OWFgeneration24((iday-1)*24+1:iday*24,:) - OWFgenStatistic24.curve;
    generationError_sum(:,iday) = OWFgeneration_sum((iday-1)*24+1:iday*24,:) - OWFgenStatistic24.curve_sum;
end
OWFgenStatistic24.errorMean = mean(generationError,3);
OWFgenStatistic24.errorStd = std(generationError,[],3);
OWFgenStatistic24.errorRelativeStd = OWFgenStatistic24.errorStd ./ OWFgenStatistic24.curve;
OWFgenStatistic24.errorMean_sum = mean(generationError_sum,2);
OWFgenStatistic24.errorStd_sum = std(generationError_sum,[],2);
OWFgenStatistic24.errorRelativeStd_sum = OWFgenStatistic24.errorStd_sum ./ OWFgenStatistic24.curve_sum;
save(projectPath('data','checkpoints','stop4.mat'))
%% case 1: comparasion of solution methods --------------------------------
clear
clc
yalmip('clear')
load(projectPath('data','checkpoints','stop4.mat'))

gppIndexSet = mpc.GEcon(:,3);                                               % 这个mpc的表格就4050MW的OWF，应该是第一轮auction的结果
nOWF = size(find(mpc.gentype==1),1);                                       
OWFindexSet = find(mpc.gentype==1);

mpc0 = mpc;
mpc0.gencost(gppIndexSet,[2,6,7]) = 0;                                      % set gpp fuel cost to 0
mpc0.gencost(OWFindexSet,[6:7]) = repmat(mpc0.gencost(5,[6:7]) / 10,[nOWF,1]); % set the operation cost of OWF to low values (but not 0)
totalOWFcapacity = 5000;                                                    % 5GW
mpc0.gen(OWFindexSet,[2:5,9]) = mpc0.gen(OWFindexSet,[2:5,9]) ./ sum(mpc0.gen(OWFindexSet,9)) * totalOWFcapacity; % normalise to 5 GW
mpc0.gencost(mpc0.gencost(:,6)==0,6) = 0.1; % set 0 cost to a very low cost

% demand (is this average demand curve?)
electricityDemandCurve = table2array(readtable(projectPath('data','inputs','system','Irish energy system data.xlsx'),'sheet','loadCurve','range','B2:B97'));
gasDemandFactor = table2array(readtable(projectPath('data','inputs','system','Irish energy system data.xlsx'),'sheet','loadCurve','range','C2:C97'));

% OWF generation dataset over the selected period
OWFgeneration = zeros(size(OWFgeneration24,1),nOWF);
for iOWF = 1:nOWF
        OWFindex = OWFindexSet(iOWF);
        approvedAreaIndex = mod(iOWF-1,nApproved)+1;
        OWFgenerationRate = OWFgeneration24(:,approvedAreaIndex) ...
            / ( OWFapprovedArea(approvedAreaIndex).nRow * OWFapprovedArea(approvedAreaIndex).nColumn * windTurbine.ratedPower); %因为容量不一样，所以要折算一下
        OWFgeneration(:,iOWF) = mpc0.gen(OWFindex,2) * OWFgenerationRate;
end
n_sample = 24*7; % sample = 1 week
operatingTimePoint = 361; % oeprating time point
sampledOWFgeneration = OWFgeneration(operatingTimePoint - n_sample-1:operatingTimePoint-1,:);
[forecastErrorSample,mu,Sigma,forecastError_test] = calculate_uncertainty_set_parameters(sampledOWFgeneration,OWFgeneration);


% modify the case according to wind and demand
for k = operatingTimePoint
    timeofday = mod(k-1,96)+1;
    mpc1 = mpc0;
    mpc1.gen(OWFindexSet,[2:5,9]) = repmat(OWFgeneration(k,:)',[1,5]); % update the realtime OWF generation
    % demand
    mpc1.bus(:,3:4) = mpc0.bus(:,3:4) * electricityDemandCurve(timeofday) ./ max(electricityDemandCurve);
    mpc1.Gbus(:,3) = mpc0.Gbus(:,3) * gasDemandFactor(timeofday);
end

% 1) compare computation performance
% for k = operatingTimePoint
%     % test run
%     % [Eresults_base,Einfo] = runopf_e(mpc1);
%     % ob.E = Einfo.problem;
%     % ob.GE = sum(GEinfo.problem);
%     % compare uncertainty handling -------------------------------------------
%     eps = 0.05;
%     % 1. deterministic
%     % [GEsolution,GEinfo] = runopf_ge(mpc1);
%     % 2. robust
%     % [GEsolution,GEinfo] = runopf_ge_rbst(mpc1,forecastErrorSample);
%     % 3. wst drcc
%     [GEsolution, GEinfo] = runopf_ge_wst_drcc(mpc1,forecastErrorSample,eps);
%     % 4. chi drcc
%     % [GEsolution, GEinfo] = runopf_ge_chi_drcc(mpc1,forecastErrorSample,eps);
%     % 5. mmt drcc
%     % [GEsolution, GEinfo] = runopf_ge_mmt_drcc(mpc1,mu,Sigma,eps);
%     % 6. variable moment-based drcc
%     % [GEsolution, GEinfo] = runopf_ge_vmmt_drcc(mpc1,mu,Sigma,forecastErrorSample,eps);
% 
%     stats = evaluate_constraints_violation(GEsolution,mpc1,forecastError_test,eps);
%     % %
%     % [HGEresults_robust{counter}, HGEinfo_robust{counter}] = ...
%     %     runopf_hge_robust(mpc1{counter},GEresults_base{counter},mu,sigma,epsilon);
% end

% 2) compare risk control for M3 and M5
epsSet = [0.01,0.05,0.10,0.15,0.20];
for i = 1:size(epsSet,2)
    eps = epsSet(i);

    [GEsolution, GEinfo] = runopf_ge_wst_drcc(mpc1,forecastErrorSample,eps);
    stats = evaluate_constraints_violation(GEsolution,mpc1,forecastError_test,eps);
    EVJP(i,1) = stats.EVJP;
    cost(i,1) = GEsolution.objfcn;
    % 5. mmt drcc
    [GEsolution, GEinfo] = runopf_ge_mmt_drcc(mpc1,mu,Sigma,eps);
    stats = evaluate_constraints_violation(GEsolution,mpc1,forecastError_test,eps);
    EVJP(i,2) = stats.EVJP;
    cost(i,2) = GEsolution.objfcn;

end

save(projectPath('data','checkpoints','stop5.mat'))
%% case 2: run drcc for different scenarios
% updated 20240414, update gas demand trend in the future
clc
clear
yalmip('clear')
load(projectPath('data','checkpoints','stop4.mat'))

nOWF = size(find(mpc0.gentype==1));                                         
OWFindexSet = find(mpc0.gentype==1);
mpc1 = mpc;
gppIndexSet = mpc1.GEcon(:,3);
mpc0.gencost(OWFindexSet,[6:7]) = repmat(mpc0.gencost(5,[6:7]) / 10,[nOWF,1]);
mpc0.ptg(:,5) = 999;                                                        % increase the capacity of PTGs

% wind and demand
electricityDemandCurve = table2array(readtable(projectPath('data','inputs','system','Irish energy system data.xlsx'),'sheet','loadCurve','range','B2:B97'));
gasDemandFactor = table2array(readtable(projectPath('data','inputs','system','Irish energy system data.xlsx'),'sheet','loadCurve','range','C2:C97'));


% short-term
totalOWFcapacity.shortTerm = 5000; % 5GW
totalOWFcapacity.midTerm = 20000;
totalOWFcapacity.longTerm = 37000;
% draw the curve of wind speed
for h = 1:24
    windPowerHist(h,:) = OWFgenStatistic24.curve_sum(h) + (OWFgenStatistic24.errorStd_sum(h) /5) .* randn(600,1)';
end

% build scenarios
counter = 0;
for k = 1:24
    timeofday = mod(k-1,24)+1;
    counter = counter + 1;
    mpcE1{counter} = mpc0; mpcE2{counter} = mpc0; mpcE3{counter} = mpc0;
    mpcE1{counter}.gen(OWFindexSet,[2:5,9]) = mpc0.gen(OWFindexSet,[2:5,9]) ./ sum(mpc0.gen(OWFindexSet,9)) * totalOWFcapacity.shortTerm; % 5 GW
    mpcE2{counter}.gen(OWFindexSet,[2:5,9]) = mpc0.gen(OWFindexSet,[2:5,9]) ./ sum(mpc0.gen(OWFindexSet,9)) * totalOWFcapacity.midTerm; % 20 GW
    mpcE3{counter}.gen(OWFindexSet,[2:5,9]) = mpc0.gen(OWFindexSet,[2:5,9]) ./ sum(mpc0.gen(OWFindexSet,9)) * totalOWFcapacity.longTerm; % 37 GW
    for iOWF = 1:nOWF
        OWFindex = OWFindexSet(iOWF);
        approvedAreaIndex = mod(iOWF-1,nApproved)+1;
        mpcE1{counter}.gen(OWFindex,[2:5,9]) = mpcE1{counter}.gen(OWFindex,[2:5,9]) / max(OWFgenStatistic24.curve(:,approvedAreaIndex)) * OWFgenStatistic24.curve(k,approvedAreaIndex);
        mpcE2{counter}.gen(OWFindex,[2:5,9]) = mpcE2{counter}.gen(OWFindex,[2:5,9]) / max(OWFgenStatistic24.curve(:,approvedAreaIndex)) * OWFgenStatistic24.curve(k,approvedAreaIndex);
        mpcE3{counter}.gen(OWFindex,[2:5,9]) = mpcE3{counter}.gen(OWFindex,[2:5,9]) / max(OWFgenStatistic24.curve(:,approvedAreaIndex)) * OWFgenStatistic24.curve(k,approvedAreaIndex);
    end
    % demand
    mpcE1{counter}.bus(:,3:4) = mpc0.bus(:,3:4) * electricityDemandCurve(timeofday) ./ max(electricityDemandCurve) * 1.022^(2030-2023); % assume demand growth 1.02
    mpcE2{counter}.bus(:,3:4) = mpc0.bus(:,3:4) * electricityDemandCurve(timeofday) ./ max(electricityDemandCurve) * 1.022^(2040-2023);
    mpcE3{counter}.bus(:,3:4) = mpc0.bus(:,3:4) * electricityDemandCurve(timeofday) ./ max(electricityDemandCurve) * 1.022^(2050-2023);
    % transmission capacity
    mpcE1{counter}.branch(:,6:8) = mpc0.branch(:,6:8) * 1.022^(2030-2023);
    mpcE2{counter}.branch(:,6:8) = mpc0.branch(:,6:8) * 1.022^(2040-2023);
    mpcE3{counter}.branch(:,6:8) = mpc0.branch(:,6:8) * 1.022^(2050-2023);
    % GE
    mpcGE1{counter} = mpcE1{counter}; mpcGE2{counter} = mpcE2{counter}; mpcGE3{counter} = mpcE3{counter};
    % set gpp cost to 0
    mpcGE1{counter}.gencost(gppIndexSet,[2,6,7]) = 0; mpcGE2{counter}.gencost(gppIndexSet,[2,6,7]) = 0; mpcGE3{counter}.gencost(gppIndexSet,[2,6,7]) = 0;
    % gas demand 天然气增长的慢点
    mpcGE1{counter}.Gbus(:,3) = mpc0.Gbus(:,3) * gasDemandFactor(timeofday) * 1.162805205;
    mpcGE2{counter}.Gbus(:,3) = mpc0.Gbus(:,3) * gasDemandFactor(timeofday) * 0.970086967;
    mpcGE3{counter}.Gbus(:,3) = mpc0.Gbus(:,3) * gasDemandFactor(timeofday) * 0.822776129;
    % update transmission factor
    mpcGE1{counter}.Gline(:,3) = mpc0.Gline(:,3) * 1.162805205;
    mpcGE2{counter}.Gline(:,3) = mpc0.Gline(:,3) * 0.970086967;
    mpcGE3{counter}.Gline(:,3) = mpc0.Gline(:,3) * 0.822776129;
    % update the gas source
    mpcGE1{counter}.Gsou(:,4) = mpc0.Gsou(:,4) * 1.162805205;
    mpcGE1{counter}.Gsou(:,4) = mpc0.Gsou(:,4) * 0.970086967;
    mpcGE1{counter}.Gsou(:,4) = mpc0.Gsou(:,4) * 0.822776129;
end
% drcc para
epsilon = 0.05;
mu1 = OWFgenStatistic24.errorMean_sum ./ OWFgenStatistic24.curve_sum * totalOWFcapacity.shortTerm;
sigma1 = OWFgenStatistic24.errorStd_sum ./ OWFgenStatistic24.curve_sum * totalOWFcapacity.shortTerm / 5;
mu2 = OWFgenStatistic24.errorMean_sum ./ OWFgenStatistic24.curve_sum * totalOWFcapacity.midTerm;
sigma2 = OWFgenStatistic24.errorStd_sum ./ OWFgenStatistic24.curve_sum * totalOWFcapacity.midTerm / 5;
mu3 = OWFgenStatistic24.errorMean_sum ./ OWFgenStatistic24.curve_sum * totalOWFcapacity.longTerm;
sigma3 = OWFgenStatistic24.errorStd_sum ./ OWFgenStatistic24.curve_sum * totalOWFcapacity.longTerm / 5;
save(projectPath('data','checkpoints','stop6.mat'))
%% deterministic, power system only
clc
clear
load(projectPath('data','checkpoints','stop6.mat'))
counter = 0;
for k = 1:24
    counter = counter + 1;
    timeofday = mod(k-1,24)+1;

    [Eresults1{counter},Einfo1{counter}] = runopf_drcc(mpcE1{counter},mu1(timeofday),sigma1(timeofday),epsilon);
    [Eresults2{counter},Einfo2{counter}] = runopf_drcc(mpcE2{counter},mu2(timeofday),sigma2(timeofday),epsilon);
    [Eresults3{counter},Einfo3{counter}] = runopf_drcc(mpcE3{counter},mu3(timeofday),sigma3(timeofday),epsilon);
    ob.E(counter,:) = [Einfo1{counter}.problem,Einfo2{counter}.problem,Einfo3{counter}.problem];
end
% organizse results for power system only
counter = 0;
for k = 1:24
    counter = counter + 1;
    Eaccomadation.OWFcapacity1(k,:) = mpcE1{counter}.gen(OWFindexSet,9)';
    Eaccomadation.OWFcapacity2(k,:) = mpcE2{counter}.gen(OWFindexSet,9)';
    Eaccomadation.OWFcapacity3(k,:) = mpcE3{counter}.gen(OWFindexSet,9)';
    Eaccomadation.OWFcapacity1_sum = sum(Eaccomadation.OWFcapacity1,2);
    Eaccomadation.OWFcapacity2_sum = sum(Eaccomadation.OWFcapacity2,2);
    Eaccomadation.OWFcapacity3_sum = sum(Eaccomadation.OWFcapacity3,2);

    Eaccomadation.OWFgeneration1(k,:) = Eresults1{counter}.Pg(OWFindexSet)';
    Eaccomadation.OWFgeneration2(k,:) = Eresults2{counter}.Pg(OWFindexSet)';
    Eaccomadation.OWFgeneration3(k,:) = Eresults3{counter}.Pg(OWFindexSet)';
    Eaccomadation.OWFgeneration1_sum = sum(Eaccomadation.OWFgeneration1,2);
    Eaccomadation.OWFgeneration2_sum = sum(Eaccomadation.OWFgeneration2,2);
    Eaccomadation.OWFgeneration3_sum = sum(Eaccomadation.OWFgeneration3,2);
   
    Eaccomadation.OWFremain1_sum = Eaccomadation.OWFcapacity1_sum - Eaccomadation.OWFgeneration1_sum;
    Eaccomadation.OWFremain2_sum = Eaccomadation.OWFcapacity2_sum - Eaccomadation.OWFgeneration2_sum;
    Eaccomadation.OWFremain3_sum = Eaccomadation.OWFcapacity3_sum - Eaccomadation.OWFgeneration3_sum;

    Eaccomadation.OWFaccomadationRate1 = Eaccomadation.OWFgeneration1./ Eaccomadation.OWFcapacity1;
    Eaccomadation.OWFaccomadationRate2 = Eaccomadation.OWFgeneration2./ Eaccomadation.OWFcapacity2;
    Eaccomadation.OWFaccomadationRate3 = Eaccomadation.OWFgeneration3./ Eaccomadation.OWFcapacity3;
    Eaccomadation.OWFaccomadationRate1_all = Eaccomadation.OWFgeneration1_sum./ Eaccomadation.OWFcapacity1_sum;
    Eaccomadation.OWFaccomadationRate2_all = Eaccomadation.OWFgeneration2_sum./ Eaccomadation.OWFcapacity2_sum;
    Eaccomadation.OWFaccomadationRate3_all = Eaccomadation.OWFgeneration3_sum./ Eaccomadation.OWFcapacity3_sum;
end


%% run HGE drcc
clc
clear
load(projectPath('data','checkpoints','stop6.mat'))
counter = 0;
for k = 1:24
    counter = counter + 1;
    timeofday = mod(k-1,24)+1;
    hymax = [];
    % short term
    [HGEresults_drcc1{counter}, HGEinfo_drcc1{counter}] = runopf_hge_drcc(mpcGE1{counter},mu1(timeofday),sigma1(timeofday),epsilon,hymax);
    % mid term
    [HGEresults_drcc2{counter}, HGEinfo_drcc2{counter}] = runopf_hge_drcc(mpcGE2{counter},mu2(timeofday),sigma2(timeofday),epsilon,hymax);
    % long term
    [HGEresults_drcc3{counter}, HGEinfo_drcc3{counter}] = runopf_hge_drcc(mpcGE3{counter},mu3(timeofday),sigma3(timeofday),epsilon,hymax);
    ob.HGE(k,1:3) = [HGEinfo_drcc1{counter}.problem, HGEinfo_drcc2{counter}.problem, HGEinfo_drcc3{counter}.problem];
end
% organizse results for HGE drcc
counter = 0;
for k = 1:24
    counter = counter + 1;
    HGEdrccAccomadation.OWFcapacity1(k,:) = mpcGE1{counter}.gen(OWFindexSet,9)';
    HGEdrccAccomadation.OWFcapacity2(k,:) = mpcGE2{counter}.gen(OWFindexSet,9)';
    HGEdrccAccomadation.OWFcapacity3(k,:) = mpcGE3{counter}.gen(OWFindexSet,9)';
    HGEdrccAccomadation.OWFcapacity1_sum = sum(HGEdrccAccomadation.OWFcapacity1,2);
    HGEdrccAccomadation.OWFcapacity2_sum = sum(HGEdrccAccomadation.OWFcapacity2,2);
    HGEdrccAccomadation.OWFcapacity3_sum = sum(HGEdrccAccomadation.OWFcapacity3,2);

    HGEdrccAccomadation.OWFgeneration1(k,:) = HGEresults_drcc1{counter}.Pg(OWFindexSet)';
    HGEdrccAccomadation.OWFgeneration2(k,:) = HGEresults_drcc2{counter}.Pg(OWFindexSet)';
    HGEdrccAccomadation.OWFgeneration3(k,:) = HGEresults_drcc3{counter}.Pg(OWFindexSet)';
    HGEdrccAccomadation.OWFgeneration1_sum = sum(HGEdrccAccomadation.OWFgeneration1,2);
    HGEdrccAccomadation.OWFgeneration2_sum = sum(HGEdrccAccomadation.OWFgeneration2,2);
    HGEdrccAccomadation.OWFgeneration3_sum = sum(HGEdrccAccomadation.OWFgeneration3,2);
   
    HGEdrccAccomadation.OWFremain1_sum = HGEdrccAccomadation.OWFcapacity1_sum - HGEdrccAccomadation.OWFgeneration1_sum;
    HGEdrccAccomadation.OWFremain2_sum = HGEdrccAccomadation.OWFcapacity2_sum - HGEdrccAccomadation.OWFgeneration2_sum;
    HGEdrccAccomadation.OWFremain3_sum = HGEdrccAccomadation.OWFcapacity3_sum - HGEdrccAccomadation.OWFgeneration3_sum;

    HGEdrccAccomadation.OWFaccomadationRate1 = HGEdrccAccomadation.OWFgeneration1./ HGEdrccAccomadation.OWFcapacity1;
    HGEdrccAccomadation.OWFaccomadationRate2 = HGEdrccAccomadation.OWFgeneration2./ HGEdrccAccomadation.OWFcapacity2;
    HGEdrccAccomadation.OWFaccomadationRate3 = HGEdrccAccomadation.OWFgeneration3./ HGEdrccAccomadation.OWFcapacity3;
    HGEdrccAccomadation.OWFaccomadationRate1_all = HGEdrccAccomadation.OWFgeneration1_sum./ HGEdrccAccomadation.OWFcapacity1_sum;
    HGEdrccAccomadation.OWFaccomadationRate2_all = HGEdrccAccomadation.OWFgeneration2_sum./ HGEdrccAccomadation.OWFcapacity2_sum;
    HGEdrccAccomadation.OWFaccomadationRate3_all = HGEdrccAccomadation.OWFgeneration3_sum./ HGEdrccAccomadation.OWFcapacity3_sum;
end
save(projectPath('data','checkpoints','stop5.mat'))
%%
clc
clear
yalmip('clear')
load(projectPath('data','checkpoints','stop6.mat'))
% mean gas composition (得重新计算，放点ptg去西边和南边）
hymax = [];
counter = 0;
for k = 1:24
    counter = counter + 1;
    mpcHGE1{counter} = mpcGE1{counter}; mpcHGE2{counter} = mpcGE2{counter}; mpcHGE3{counter} = mpcGE3{counter};
    addPTG = [
        82	6	    0.95	0	999 %
        138	6	    0.95	0	999 %
        126	816	    0.95	0	999 %
        55	334	    0.95	0	999 %
        55	352	    0.95	0	999
        55	1111	0.95	0	999
        55	582	    0.95	0	999
    ];
    mpcHGE1{counter}.ptg = [mpcHGE1{counter}.ptg; addPTG];
    mpcHGE2{counter}.ptg = [mpcHGE2{counter}.ptg; addPTG];
    mpcHGE3{counter}.ptg = [mpcHGE3{counter}.ptg; addPTG];
end

counter = 0;
for k = 1:24
    counter = counter + 1;
    timeofday = mod(k-1,24)+1;
    % short term
    [HGEresults_drcc1_new{counter}, HGEinfo_drcc1{counter}] = runopf_hge_drcc(mpcHGE1{counter},mu1(timeofday),sigma1(timeofday),epsilon,hymax);
    % mid term
    [HGEresults_drcc2_new{counter}, HGEinfo_drcc2{counter}] = runopf_hge_drcc(mpcHGE2{counter},mu2(timeofday),sigma2(timeofday),epsilon,hymax);
    % long term
    [HGEresults_drcc3_new{counter}, HGEinfo_drcc3{counter}] = runopf_hge_drcc(mpcHGE3{counter},mu3(timeofday),sigma3(timeofday),epsilon,hymax);
    ob.HGE(k,1:3) = [HGEinfo_drcc1{counter}.problem, HGEinfo_drcc2{counter}.problem, HGEinfo_drcc3{counter}.problem];
end

nGb = size(mpc0.Gbus,1);
counter = 0;
HGEdrccAccomadation_new.meanHyComposition1 = zeros(nGb,1); HGEdrccAccomadation_new.meanHyComposition2 = zeros(nGb,1); HGEdrccAccomadation_new.meanHyComposition3 = zeros(nGb,1);
HGEdrccAccomadation_new.meanQptg_all = zeros(3,2);
for k = 1:24
    counter = counter + 1;
    % offshore wind accomadation
    HGEdrccAccomadation_new.OWFcapacity1(k,:) = mpcHGE1{counter}.gen(OWFindexSet,9)';
    HGEdrccAccomadation_new.OWFcapacity2(k,:) = mpcHGE2{counter}.gen(OWFindexSet,9)';
    HGEdrccAccomadation_new.OWFcapacity3(k,:) = mpcHGE3{counter}.gen(OWFindexSet,9)';
    HGEdrccAccomadation_new.OWFcapacity1_sum = sum(HGEdrccAccomadation_new.OWFcapacity1,2);
    HGEdrccAccomadation_new.OWFcapacity2_sum = sum(HGEdrccAccomadation_new.OWFcapacity2,2);
    HGEdrccAccomadation_new.OWFcapacity3_sum = sum(HGEdrccAccomadation_new.OWFcapacity3,2);

    HGEdrccAccomadation_new.OWFgeneration1(k,:) = HGEresults_drcc1_new{counter}.Pg(OWFindexSet)';
    HGEdrccAccomadation_new.OWFgeneration2(k,:) = HGEresults_drcc2_new{counter}.Pg(OWFindexSet)';
    HGEdrccAccomadation_new.OWFgeneration3(k,:) = HGEresults_drcc3_new{counter}.Pg(OWFindexSet)';
    HGEdrccAccomadation_new.OWFgeneration1_sum = sum(HGEdrccAccomadation_new.OWFgeneration1,2);
    HGEdrccAccomadation_new.OWFgeneration2_sum = sum(HGEdrccAccomadation_new.OWFgeneration2,2);
    HGEdrccAccomadation_new.OWFgeneration3_sum = sum(HGEdrccAccomadation_new.OWFgeneration3,2);
   
    HGEdrccAccomadation_new.OWFremain1_sum = HGEdrccAccomadation_new.OWFcapacity1_sum - HGEdrccAccomadation_new.OWFgeneration1_sum;
    HGEdrccAccomadation_new.OWFremain2_sum = HGEdrccAccomadation_new.OWFcapacity2_sum - HGEdrccAccomadation_new.OWFgeneration2_sum;
    HGEdrccAccomadation_new.OWFremain3_sum = HGEdrccAccomadation_new.OWFcapacity3_sum - HGEdrccAccomadation_new.OWFgeneration3_sum;

    HGEdrccAccomadation_new.OWFaccomadationRate1 = HGEdrccAccomadation_new.OWFgeneration1./ HGEdrccAccomadation_new.OWFcapacity1;
    HGEdrccAccomadation_new.OWFaccomadationRate2 = HGEdrccAccomadation_new.OWFgeneration2./ HGEdrccAccomadation_new.OWFcapacity2;
    HGEdrccAccomadation_new.OWFaccomadationRate3 = HGEdrccAccomadation_new.OWFgeneration3./ HGEdrccAccomadation_new.OWFcapacity3;
    HGEdrccAccomadation_new.OWFaccomadationRate1_all = HGEdrccAccomadation_new.OWFgeneration1_sum./ HGEdrccAccomadation_new.OWFcapacity1_sum;
    HGEdrccAccomadation_new.OWFaccomadationRate2_all = HGEdrccAccomadation_new.OWFgeneration2_sum./ HGEdrccAccomadation_new.OWFcapacity2_sum;
    HGEdrccAccomadation_new.OWFaccomadationRate3_all = HGEdrccAccomadation_new.OWFgeneration3_sum./ HGEdrccAccomadation_new.OWFcapacity3_sum;
    % gas composition
    HGEdrccAccomadation_new.meanHyComposition1 = HGEdrccAccomadation_new.meanHyComposition1 + HGEresults_drcc1_new{counter}.gasComposition(:,2);
    HGEdrccAccomadation_new.meanHyComposition2 = HGEdrccAccomadation_new.meanHyComposition2 + HGEresults_drcc1_new{counter}.gasComposition(:,2);
    HGEdrccAccomadation_new.meanHyComposition3 = HGEdrccAccomadation_new.meanHyComposition3 + HGEresults_drcc1_new{counter}.gasComposition(:,2);
    % Qptg
    HGEdrccAccomadation_new.meanQptg_all(1,:) = HGEdrccAccomadation_new.meanQptg_all(1,:) + sum(HGEresults_drcc1_new{counter}.Qptg);
    HGEdrccAccomadation_new.meanQptg_all(2,:) = HGEdrccAccomadation_new.meanQptg_all(2,:) + sum(HGEresults_drcc2_new{counter}.Qptg);
    HGEdrccAccomadation_new.meanQptg_all(3,:) = HGEdrccAccomadation_new.meanQptg_all(3,:) + sum(HGEresults_drcc3_new{counter}.Qptg);
end
HGEdrccAccomadation_new.meanHyComposition1 = HGEdrccAccomadation_new.meanHyComposition1 / 24;
HGEdrccAccomadation_new.meanHyComposition2 = HGEdrccAccomadation_new.meanHyComposition2 / 24;
HGEdrccAccomadation_new.meanHyComposition3 = HGEdrccAccomadation_new.meanHyComposition3 / 24;
HGEdrccAccomadation_new.meanQptg_all = HGEdrccAccomadation_new.meanQptg_all / 24;


save
%% sensitivity analysis
clc
clear
yalmip('clear')
load(projectPath('data','checkpoints','stop6.mat'))
% mean gas composition (得重新计算，放点ptg去西边和南边）
counter = 0;
for k = 1:24
    counter = counter + 1;
    mpcHGE1{counter} = mpcGE1{counter}; mpcHGE2{counter} = mpcGE2{counter}; mpcHGE3{counter} = mpcGE3{counter};
    addPTG = [
        82	6	    0.95	0	999 %
        138	6	    0.95	0	999 %
        126	816	    0.95	0	999 %
        55	334	    0.95	0	999 %
        55	352	    0.95	0	999
        55	1111	0.95	0	999
        55	582	    0.95	0	999
    ];
    mpcHGE1{counter}.ptg = [mpcHGE1{counter}.ptg; addPTG];
    mpcHGE2{counter}.ptg = [mpcHGE2{counter}.ptg; addPTG];
    mpcHGE3{counter}.ptg = [mpcHGE3{counter}.ptg; addPTG];
end
% % pathway 1
% for i = 1:10
%     epsilon = 0.05;
%     hymax = [];
%     mu = mu3(1);
%     sigma(i) = sigma3(1) * i;
%     [HGEresults{i,1}, HGEinfo{i,1}] = runopf_hge_drcc(mpcHGE3{1},mu,sigma(i),epsilon,hymax);
% end
% % pathway 2
% for i = 1:10
%     epsilon = 0.05;
%     hymax = i*0.1;
%     mu = mu3(1);
%     sigma = sigma3(1);
%     [HGEresults{i,2}, HGEinfo{i,2}] = runopf_hge_drcc(mpcHGE3{1},mu,sigma,epsilon,hymax);
% end
% results
% HGEsensitivity.accomadationRate = zeros(2,10);
% for i = 1:2
%     for j = 1:10
%         OWFgeneration_total = sum(HGEresults{j,i}.Pg(OWFindexSet));
%         OWFcapacity_total = sum(mpcHGE3{1}.gen(OWFindexSet,9));
%         HGEsensitivity.accomadationRate(i,j) = OWFgeneration_total / OWFcapacity_total;
%         HGEsensitivity.Pptg_sum(i,j) = sum(HGEresults{j,i}.Pptg);
%         HGEsensitivity.totalCost(i,j) = sum(HGEresults{j,i}.totalCost);
%         HGEsensitivity.genAndLCeCost(i,j) = sum(HGEresults{j,i}.genAndLCeCost);
%         HGEsensitivity.gasPurchasingCost(i,j) = sum(HGEresults{j,i}.gasPurchasingCost);
%     end
% end

% pathway 1 and 2
for i = 1:10
    for j = 1:10
        epsilon = 0.05;
        hymax = i*0.1;
        mu = mu3(1);
        sigma = sigma3(1) * j;
        [HGEresults{i,j}, HGEinfo{i,j}] = runopf_hge_drcc(mpcHGE3{1},mu,sigma,epsilon,hymax);
    end
end

HGEsensitivity.accomadationRate = zeros(2,10);
for i = 1:10
    for j = 1:10
        OWFgeneration_total = sum(HGEresults{i,j}.Pg(OWFindexSet));
        OWFcapacity_total = sum(mpcHGE3{1}.gen(OWFindexSet,9));
        HGEsensitivity.accomadationRate(i,j) = OWFgeneration_total / OWFcapacity_total;
        HGEsensitivity.Pptg_sum(i,j) = sum(HGEresults{i,j}.Pptg);
        HGEsensitivity.totalCost(i,j) = sum(HGEresults{i,j}.totalCost);
        HGEsensitivity.genAndLCeCost(i,j) = sum(HGEresults{i,j}.genAndLCeCost);
        HGEsensitivity.gasPurchasingCost(i,j) = sum(HGEresults{i,j}.gasPurchasingCost);
        HGEsensitivity.carbonEmission(i,j) = (sum(HGEresults{i,j}.PGs) + sum(HGEresults{i,j}.Qptg(:,1))) / 16 * 44;
    end
end


save
