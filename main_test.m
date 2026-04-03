% this file is to validate the benchmark operating states of irish energy
% systems. Please include matpower 7.1, yalmip, and Gurobi v11 before
% using. 
% - Sheng Wang, 2023 12 15
clc
clear
yalmip('clear')
%% load data
systemDataFile = projectPath('data','inputs','system','Irish energy system data.xlsx');
mpc.bus = table2array(readtable(systemDataFile,'sheet','bus','range','A2:M1197'));
mpc.gen = table2array(readtable(systemDataFile,'sheet','gen','range','A2:U140'));
mpc.branch = table2array(readtable(systemDataFile,'sheet','branch','range','A2:M1369'));
mpc.gencost = table2array(readtable(systemDataFile,'sheet','gencost','range','A2:G140'));
mpc.Gbus = table2array(readtable(systemDataFile,'sheet','Gbus','range','A2:F145'));
mpc.Gline = table2array(readtable(systemDataFile,'sheet','Gline','range','A2:H145'));
mpc.Gsou = table2array(readtable(systemDataFile,'sheet','Gsou','range','A2:D3'));
mpc.Gcost = table2array(readtable(systemDataFile,'sheet','Gcost','range','A2:A3'));
mpc.GEcon = table2array(readtable(systemDataFile,'sheet','GEcon','range','A2:D14'));
mpc.ptg = table2array(readtable(systemDataFile,'sheet','ptg','range','A2:E8'));
mpc.version = '2';
mpc.baseMVA = 100;
%%
% % case 1: electricity system optimal power flow
[information, electricityResults] = runopf(mpc);
% % case 2: gas system optimal power flow
[information, gasResults] = runopf_gas(mpc);
% % case 3: coordinated optimal power flow of electricity and gas systems
mpc.gencost(mpc.GEcon(:,3),:) = 0; % set the generation cost of gas-fired power plant to zero (because the fuel cost has been counted in the gas system cost)
[information, electricityGasResults] = runopf_ge(mpc);
