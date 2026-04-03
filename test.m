clc
clear
% this script will read-in and plot NetCDF data in MATLAB
% Define data file. Add directory path if necessary.
% file = 'MERRA2_100.instM_2d_asm_Nx.199001.nc4';
weatherFile = projectPath('data','inputs','weather_samples','MERRA2_400.tavg1_2d_flx_Nx.20220101.SUB.nc');
psseRawFile = projectPath('data','inputs','network','PSSERaw_from_Dutrain.raw');
file = weatherFile;
ncinfo(file);

% Uncomment to display metadata information
%ncdisp(file);

% read-in variables
% var1 = ncread(file, 'Var_V50M');
var1 = ncread(file, 'ULML');
% MATLAB data is oriented by (Y,X), but the data was written as (X,Y).
% The rot90 and fliplr function correctly orient the data.
var1 = rot90(fliplr(var1));
lats = ncread(file, 'lat');
lons = ncread(file, 'lon');


% ========== Plot Data ========================
pcolor(lons,lats,var1(:,:,1));
shading flat
c = colorbar;
ylabel(c,'Dobsons')
load coast
hold on
plot(long,lat,'black')
title('MERRA-2 Total Column Ozone')
xlabel('degrees longitude')
ylabel('degrees latitude')

%%
runopf(PSSERaw_from_Dutrain);
%%
mpc1 = PSSERaw_from_Dutrain;
[mpc, warnings] = psse2mpc(psseRawFile);
find(mpc.gen(:,1) == 776)
%%
k = 60;
mpc1{k}.branch(:,6:8) = mpc0.branch(:,6:8) * 100;
[results,info] = runopf(mpc1{k});
%% 通过测试发现，风电多的时候，ptg还是有产气的
mpc.gen(131:139,9) = mpc0.gen(131:139,9) * 1;
[result1,info1] = runopf_ge(mpc);
[HGEresults, HGEinfo] = runopf_hge(mpc,result1);
