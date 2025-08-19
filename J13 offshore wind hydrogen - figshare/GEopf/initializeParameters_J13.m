function [GCV, M, fs, a, R, T_stp, Prs_stp, Z, T_gas, eta, CDF,rho_stp] = initializeParameters_J13()
composition_ng = [91.66 	3.88 	0.46 	0.13 	0.00 	1.54 	2.33] / 100;
% CH4, C2H6, C3H8, C4H10, H2, N2, CO2
GCV.CH4 = 3.85 * 1e7;
GCV.C2H6 = 5.74 * 1e7;
GCV.C3H8 = 5.59 * 1e7;
GCV.C4H10 = 12.43 * 1e7;
GCV.hy = 12.75 * 1e6;      % J/m3
GCV.N2 = 0;     % J/m3
GCV.CO2 = 0;
GCV.all = [GCV.CH4, GCV.C2H6, GCV.C3H8, GCV.C4H10, GCV.hy, GCV.N2, GCV.CO2];
% GCV.ng_ref = 41.04 * 1e6;     % J/m3
GCV.ng = GCV.all * composition_ng';
% 
M.CH4 = 16 * 1e-3;
M.C2H6 = 30 * 1e-3;
M.C3H8 = 44 * 1e-3;
M.C4H10 = 58 * 1e-3;
M.hy = 2 * 1e-3;           % kg/mol
M.N2 = 28 * 1e-3;
M.CO2 = 44 * 1e-3;
M.all = [M.CH4, M.C2H6, M.C3H8, M.C4H10, M.hy, M.N2, M.CO2];
M.air = 29 * 1e-3;         % kg/mol
% M.ng_ref = 17.478 * 1e-3;     % kg/mol
M.ng = M.all * composition_ng';
%  ﬂame speed 
fs.CH4 = 148;
fs.C2H6 = 301;
fs.C3H8 = 398;
fs.C4H10 = 513;
fs.hy = 339;           % kg/mol
fs.N2 = 0;
fs.CO2 = 0;
fs.All = [fs.CH4, fs.C2H6, fs.C3H8, fs.C4H10, fs.hy, fs.N2, fs.CO2];
fs.ng = fs.All * composition_ng';
% a
a_CH4 = 0.3;
a_C2H6 = 0.75;
a_C3H8 = 0.9759;
a_C4H10 = 1.0928;
a_hy = 0;         
a_N2 = 0.699;
a_CO2 = 0.9759;
a.All = [a_CH4, a_C2H6, a_C3H8, a_C4H10, a_hy, a_N2, a_CO2];
%
Rgas = 8.31446261815324; % J/(mol*K) 这个数对所与气体是不变的，但是有时候R要化成kg的单位，所以与相对分子质量有关
R.air = Rgas / M.air; % J/(kg*K) 约287
R.all = Rgas ./ M.all;
[R.CH4, R.C2H6, R.C3H8, R.C4H10, R.hy, R.N2, R.CO2] = deal(...
    Rgas./M.CH4, Rgas./M.C2H6, Rgas./M.C3H8, Rgas./M.C4H10, Rgas./M.hy, Rgas./M.N2, Rgas./M.CO2);
R.ng = Rgas / M.ng;
T_stp = 288;               % K
Prs_stp = 101325;          % Pa
Z = 1;                     % dimenssionless
T_gas = 281.15;            % K
eta.electrolysis = 0.7;    % from the energy perspective, the effciency is about 80%
eta.methanation = 0.8;
eta.GFU = 0.4211;           % from the energy perspective, 从1/200换算而来
rho_stp = Prs_stp / (Z * R.ng * T_stp); % kg/m3, not accurate value, Prs_stp / (Z * R_ng * T_stp);
%% CDF
CDF.electricity = 1e4; % MW/hour, 大概数值，从jia文章中拿的
CDF.gas = CDF.electricity * GCV.ng / 3600 / 24;
end