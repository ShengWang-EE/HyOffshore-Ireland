function [HGEsolution, HGEinfo] = runopf_hge_drcc(mpc,GEsolution,mu,sigma,epsilon)
%RUNOPF_HGE Summary of this function goes here
%   Detailed explanation goes here
%
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
baseMVA = 100;
il = find(mpc.branch(:, RATE_A) ~= 0 & mpc.branch(:, RATE_A) < 1e10);
%
nb   = size(mpc.bus, 1);    %% number of buses
nGb  = size(mpc.Gbus,1); % number of gas bus
nGl = size(mpc.Gline,1);
ng = size(mpc.gen,1);
nl = size(mpc.branch,1);
nGpp = size(mpc.GEcon,1);
nGs = size(mpc.Gsou,1);
nGd = size(find(mpc.Gbus(:,3)~=0),1);
nPTG = size(mpc.ptg,1);
iGd = find(mpc.Gbus(:,3)~=0);
iGppGb = mpc.GEcon(:,1); % gas bus index of the GPPs
mpc.gasCompositionForGasSource = repmat([1,0],[nGs,1]);
OWFindexSet = find(mpc.gentype==1);
nOWF = size(OWFindexSet,1);
% natural gas, hydrogen
nGasType = 2;
[GCV, M, fs, a, R, T_stp, Prs_stp, Z_ref, T_gas, eta, CDF,rho_stp] = initializeParameters_J13();
%%
% secuiity limit
FSlimit = 0.2; 
% WImin = 46.5; WImax = 52.85; WIlimit = (WImax-WImin) / (WImax+WImin); % contigent
WImin = 47.2; WImax = 51.41; WIlimit = (WImax-WImin) / (WImax+WImin); % normal
% WImin = 45.8; WImax = 54.29; WIlimit = (WImax-WImin) / (WImax+WImin); % further relax
% paras for sequential programming
alpha_PHI = 1e-1;alpha_x = 1e-1; alpha_Qd = 1e-1;
lambda = 10; % coefficient for gas flow
multiplier = 2; alpha_PHI_max = 1e4; alpha_x_max = 1e4; alpha_Qd_max = 1e4;
iterationMax = 1;
gap = 1e-3; % optimization convergence criterion

%% initial value
PGs_ref = GEsolution.PGs;
GCV_ng = GCV.ng;
S_ng = M.ng/M.air;
R_ng = R.ng;
WI_ng = GCV_ng / sqrt(S_ng);
FS_ng = fs.ng;
GCVall = [GCV.ng, GCV.hy];
FSall = [fs.ng, fs.hy];
Mall = [M.ng, M.hy];
Rall = [R.ng,R.hy];
Sall = Mall / M.air;
% 需要提前计算并迭代的变量：
S_pipeline_ref = S_ng; %
R_pipeline_ref = R_ng;
gasComposition_ref = repmat([1,0],[nGb,1]); % initial gas composition is [1,0]
Qd_ref = repmat(mpc.Gbus(iGd,3),[1,nGasType]) .* gasComposition_ref(iGd,:);
Qgpp_ref = repmat(GEsolution.Qgpp,[1,nGasType]) .* gasComposition_ref(iGppGb,:);
gasFlow_sum_ref = GEsolution.Gf;

PTGbus = mpc.ptg(:,1) ; 
CgsPTG = sparse(PTGbus, (1:nPTG)', 1, nGb, nPTG); % connection matrix
QptgInbus_ref = CgsPTG * GEsolution.Qptg; 

signGf = GEsolution.gamma;
W_ref = nodalGasInjection_hge(PGs_ref,QptgInbus_ref,signGf, gasFlow_sum_ref,mpc,nGs,nGb,nGl);
gasLoad_ng = mpc.Gbus(iGd,3);
%%
solverTime = zeros(iterationMax,1);
for v = 1:iterationMax
    yalmip('clear')
%% state variables
Prs = sdpvar(nGb,1); % bar^2
PGs = sdpvar(nGs,1); % Mm3/day
Qd = sdpvar(nGd,nGasType);% Mm3/day
Qptg = sdpvar(nPTG,2); % [ methane; hydrogen ] % Mm3/day
Pptg = sdpvar(nPTG,1); % electricity consumption, 1/100 MW
Pg = sdpvar(ng,1); % include TPP, GPP and renewable generators, 1/100 MW
Qgpp = sdpvar(nGpp, nGasType); % Mm3/day
Va = sdpvar(nb,1);
gasComposition = sdpvar(nGb,nGasType); 
gasFlow = sdpvar(nGl,nGasType);% Mm3/day
sigma_PHI = sdpvar(nGl,1); % error limit for gas flow
sigma_x = sdpvar(nGb, nGasType); % error for gas composition
sigma_Qd = sdpvar(nGd, nGasType);
sigma_Qgpp = sdpvar(nGpp,nGasType);
sigma_Gf = sdpvar(nGl,nGasType);
gamma = mpc.Gline(:,9);
%
LCg = 0;
% drcc vars
u_PGs = sdpvar(nGs,1); v_PGs = sdpvar(nGs,1); 
phi_PGs = sdpvar(nGs,1); % affine mapping factor, 假设我们用所有风机误差的总和来作为调节依据，这就简化了，其他也有文章这么做的 
phi_Qd = sdpvar(nGd,nGasType);
phi_gasComposition = sdpvar(nGb,nGasType);
phi_Qgpp = sdpvar(nGpp, nGasType);
phi_Qptg = sdpvar(nPTG,2);
phi_Pptg = sdpvar(nPTG,1);
u_Pptg = sdpvar(nPTG,1); v_Pptg = sdpvar(nPTG,1); 
phi_gasFlow = sdpvar(nGl,nGasType);
phi_Pg = sdpvar(ng,1); u_Pg = sdpvar(ng,1); v_Pg = sdpvar(ng,1);
phi_Va = sdpvar(nb,1); u_Va = sdpvar(nl,1); v_Va = sdpvar(nl,1); % 这里的u,v的维度应该不是变量的维度，而是参与的约束的维度，其他的只是刚好变量和约束的维度一样
phi_Prs = sdpvar(nGb,1); u_Prs = sdpvar(nGb,1); v_Prs = sdpvar(nGb,1);
phi_windCapacity = sdpvar(nOWF,1);
%% bounds
Prsmin = mpc.Gbus(:,5); Prsmax = mpc.Gbus(:,6); % bar
PGsmin = mpc.Gsou(:,3); PGsmax = mpc.Gsou(:,4); % Mm3/day
Qptgmin = mpc.ptg(:,4); Qptgmax = mpc.ptg(:,5);
QptgMax_hydrogen = Qptgmax;
Pgmin = mpc.gen(:, PMIN) / baseMVA *0; %Pgmin is set to zero
Pgmax = mpc.gen(:, PMAX) / baseMVA;
WImin = WI_ng/1e6 * (1-WIlimit); WImax = WI_ng/1e6 * (1+WIlimit);
FSmin = FS_ng * (1-FSlimit); FSmax = FS_ng * (1+FSlimit);

refs = find(mpc.bus(:, BUS_TYPE) == REF);
Vau = Inf(nb, 1);       %% voltage angle limits
Val = -Vau;
Vau(refs) = 1;   %% voltage angle reference constraints
Val(refs) = 1;
gasFlowMax = mpc.Gline(:,5);
%% pre-calculation
% calculate the gas flow of each pipeline
newC = mpc.Gline(:,3) .* sqrt(R_ng) ./ sqrt(R_pipeline_ref);
% newC = mpc.Gline(:,3);
%% constraints
sigma = sigma/baseMVA;
mu = mu/baseMVA;
% epsilon = 0.2;
% 1 gas source cons DRCC
gasSourceCons_drcc = [ ...
    PGsmin <= PGs <= PGsmax;
%     u_PGs.^2 + phi_PGs.^2 * sigma <= ( (PGsmax-PGsmin)/2 - v_PGs).^2 * epsilon;
%     -(u_PGs+v_PGs) <= PGs - (PGsmax+PGsmin)/2 <= u_PGs+v_PGs;
%     u_PGs >= 0;
%     0 <= v_PGs <= (PGsmax-PGsmin)/2;
%
%     PGs + sqrt((1-epsilon)/epsilon) * sigma * phi_PGs .^2 <= PGsmax;
%     PGs >= PGsmin;
    ];

% gas demand cons
energyDemand = mpc.Gbus(iGd,3) * GCV_ng; % energy need of these gas bus
gasDemandCons_drcc = [
    Qd * GCVall'/1e9 == energyDemand/1e9;
    phi_Qd * GCVall'/1e9 == 0

    gasComposition(iGd,:) .* repmat(sum(Qd,2),[1,nGasType]) - Qd == 0; 
    gasComposition(iGd,:) .* repmat(sum(phi_Qd,2),[1,nGasType]) + phi_gasComposition(iGd,:) .* repmat(sum(Qd,2),[1,nGasType]) == phi_Qd;
    phi_gasComposition(iGd,:) .* repmat(sum(phi_Qd,2),[1,nGasType]) == 0;

    Qd >= 0;
    sigma_Qd >= 0;
    ]:'gasDemandCons'; % nonconvex ---------------------------------------

% nodal gas flow balance cons
nodalGasFlowBalanceCons_drcc = [
%     consfcn_nodalGasFlowBalance_hge(PGs,Qd,Qgpp,Qptg, gasFlow,mpc,nGasType,nGpp,nGd,iGd) == 0;
    consfcn_nodalGasFlowBalance_hge_drcc(PGs,Qd,Qgpp,Qptg, gasFlow,phi_PGs, phi_Qd, phi_Qgpp, phi_Qptg, phi_gasFlow,mpc,nGasType,nGpp,nGd,iGd) == 0;
    ]:'nodalGasFlowBalanceCons';

% ptg
% QptgMax_hydrogen = QptgMax_hydrogen * 1000;
PTGcons_drcc = [
    ( Qptg(:,1) * 1e6/24/3600 *GCV.CH4 / eta.methanation + Qptg(:,2) * 1e6/24/3600 * GCV.hy ...
        ) /1e6 == Pptg * baseMVA / eta.electrolysis; % w
%     ( phi_Qptg(:,1) * 1e6/24/3600 *GCV.CH4 / eta.methanation + phi_Qptg(:,2) * 1e6/24/3600 * GCV.hy ...
%         ) /1e6 == phi_Pptg / eta.electrolysis; % w

%     0 <= Pptg * baseMVA <= QptgMax_hydrogen /24/3600 * GCV.hy * eta.electrolysis; % 如果全用来制氢，
%     u_Pptg.^2 + sigma .* phi_Pptg.^2 <= (QptgMax_hydrogen /24/3600 * GCV.hy * eta.electrolysis/baseMVA/2 - v_Pptg).^2 * epsilon;
%     - (u_Pptg + v_Pptg) <= Pptg - QptgMax_hydrogen /24/3600 * GCV.hy * eta.electrolysis/baseMVA/2 <= u_Pptg + v_Pptg;
%     u_Pptg >= 0;
%     0 <= v_Pptg <= QptgMax_hydrogen /24/3600 * GCV.hy * eta.electrolysis/baseMVA/2;
%
    Pptg + sqrt((1-epsilon)/epsilon) * sigma * phi_Pptg .^2 <= QptgMax_hydrogen /24/3600 * GCV.hy * eta.electrolysis / baseMVA;
    Pptg >= 0;
    
    0 <= Qptg;
    ];

% gpp
Pgpp = Pg(mpc.GEcon(:,3));% MW
phi_Pgpp = phi_Pg(mpc.GEcon(:,3));
GPPcons_drcc = [
    gasComposition(iGppGb,:) .* repmat(sum(Qgpp,2),[1,nGasType]) - Qgpp == 0;
%     gasComposition(iGppGb,:) .* repmat(sum(phi_Qgpp,2),[1,nGasType]) + phi_gasComposition(iGppGb,:) .* repmat(sum(Qgpp,2),[1,nGasType]) == phi_Qgpp;
%     phi_gasComposition(iGppGb,:) .* repmat(sum(phi_Qgpp,2),[1,nGasType]) == 0;
    sigma_Qgpp >= 0;
    Pgpp * baseMVA == Qgpp/24/3600 * GCVall'*eta.GFU; %Mm3对应MW
%     phi_Pgpp == phi_Qgpp/24/3600 * GCVall'*eta.GFU; 
    Qgpp >= 0;
    ]:'GPPcons';

% electricity flow
% 是不是DRCC里面风电就不能给上下限了？
Pg_nonwind = Pg; Pg_nonwind(OWFindexSet) = [];
Pgmax_nonwind = Pgmax; Pgmax_nonwind(OWFindexSet) = [];
Pgmin_nonwind = Pgmin; Pgmin_nonwind(OWFindexSet) = [];
u_Pg_nonwind = u_Pg; u_Pg_nonwind(OWFindexSet) = [];
v_Pg_nonwind = v_Pg; v_Pg_nonwind(OWFindexSet) = [];
phi_Pg_nonwind = phi_Pg; phi_Pg_nonwind(OWFindexSet) = [];
Pg_wind = Pg(OWFindexSet);
Pgmax_wind = Pgmax(OWFindexSet); Pgmin_wind = Pgmin(OWFindexSet);
u_Pg_wind = u_Pg(OWFindexSet); v_Pg_wind = v_Pg(OWFindexSet); phi_Pg_wind = phi_Pg(OWFindexSet);
[B, Bf, Pbusinj, Pfinj] = makeBdc(mpc.baseMVA, mpc.bus, mpc.branch);
upf = mpc.branch(il, RATE_A) / mpc.baseMVA - Pfinj(il);
upt = mpc.branch(il, RATE_A) / mpc.baseMVA + Pfinj(il);
electricityCons_drcc = [
    - upf <= Bf(il,:)*Va <= upf;
%     u_Va.^2 + sigma .* (Bf(il,:)*phi_Va).^2 <= (upt-v_Va).^2 * epsilon;
%     -(u_Va+v_Va) <= Bf(il,:)*Va <= u_Va+v_Va;
%     u_Va >= 0;
%     0 <= v_Va <= upt;

%     Pgmin <= Pg <= Pgmax;
%     u_Pg.^2 + phi_Pg.^2 * sigma <= ( (Pgmax-Pgmin)/2 - v_Pg).^2 * epsilon;
%     -(u_Pg+v_Pg) <= Pg - (Pgmax+Pgmin)/2 <= u_Pg+v_Pg;
%     u_Pg >= 0;
%     0 <= v_Pg <= (Pgmax-Pgmin)/2;
%     sum(phi_Pg(OWFindexSet))  == 1;
    % correction
%     Pgmin_nonwind <= Pg_nonwind <= Pgmax_nonwind;
%     u_Pg_nonwind.^2 + phi_Pg_nonwind.^2 * sigma <= ( (Pgmax_nonwind-Pgmin_nonwind)/2 - v_Pg_nonwind).^2 * epsilon;
%     -(u_Pg_nonwind+v_Pg_nonwind) <= Pg_nonwind - (Pgmax_nonwind+Pgmin_nonwind)/2 <= u_Pg_nonwind+v_Pg_nonwind;
%     u_Pg_nonwind >= 0;
%     0 <= v_Pg_nonwind <= (Pgmax_nonwind-Pgmin_nonwind)/2;
    %
    Pg_nonwind + sqrt((1-epsilon)/epsilon) * sigma * phi_Pg_nonwind .^2 <= Pgmax_nonwind;
    Pg_nonwind >= 0;

%     Pgmin_wind <= Pg_wind <= Pgmax_wind;
    %
%     u_Pg_wind.^2 + phi_Pg_wind.^2 * sigma <= ( (Pgmax_wind-Pgmin_wind)/2 - v_Pg_wind).^2 * epsilon;
%     -(u_Pg_wind+v_Pg_wind) <= Pg_wind - (Pgmax_wind+Pgmin_wind)/2 <= u_Pg_wind+v_Pg_wind;
%     u_Pg_wind >= 0;
%     0 <= v_Pg_wind <= (Pgmax_wind-Pgmin_wind)/2;
    sum(phi_windCapacity) == 1;
    % 
    Pg_wind + sqrt((1-epsilon)/epsilon) * sigma * phi_windCapacity .^2 <= Pgmax_wind;
    Pg_wind >= 0;

    Va(refs) == 0; % 除了slackbus外，其他相角都没约束。但是一些solver在处理inf的上下限的时候有问题
    ]:'electricityCons';

electricityBalanceCons_drcc = [...
%     consfcn_electricPowerBalance_hge(Va,Pg,Pptg,mpc) == 0;
    consfcn_electricPowerBalance_hge_drcc(Va,Pg,Pptg,phi_Va, phi_Pg, phi_Pptg,mpc) == 0;
    phi_Pg_wind == 0;
    ]:'electricityBalanceCons';

% security index: wobbe index, FS, SG
GCV_nodal = gasComposition * GCVall';
S_nodal = gasComposition * Mall' / M.air;
sqrtS = 0.5 * (S_nodal/sqrt(S_ng) + sqrt(S_ng));
WI_nodal = GCV_nodal ./ sqrtS; % 如果不能自动转化，那就手动化一下
FSnodal = gasComposition * FSall';
gasSecurityCons_drcc = [
    WImin * sqrtS <= GCV_nodal/1e6 <= WImax * sqrtS;
%     (-2*sqrt(S_ng) * phi_gasComposition * GCVall'/1e6 + WImin * phi_gasComposition * Sall') * mu ...
%         + (-2*sqrt(S_ng) * gasComposition * GCVall'/1e6 + WImin * gasComposition * Sall') ...
%         + sqrt((1-epsilon/2)/(epsilon/2)) * sqrt(sigma) * (-2*sqrt(S_ng) * phi_gasComposition * GCVall'/1e6 + WImin * phi_gasComposition * Sall') ...
%         <= -WImin * S_ng;

        FSmin <= FSnodal <= FSmax;
%     phi_gasComposition * FSall' * mu + gasComposition * FSall' + sqrt((1-epsilon/2)/(epsilon/2)) * sqrt(sigma) ...
%         * phi_gasComposition * FSall' <= FSmax;
        S_nodal <= 0.7;
    ]:'gasSecurityCons';
% other security cons
otherSecurityCons_drcc = [
    0 <= gasComposition <= 1;
    ]:'otherSecurityCons';

% SOC reformulation for gas flow
FB = mpc.Gline(:,1); TB = mpc.Gline(:,2);
gasFlow_sum = sum(gasFlow,2);
phi_gasFlowSum = sum(phi_gasFlow,2);
% gamma = {1,-1}
PHI = gamma .* (Prs(FB).^2-Prs(TB).^2);
gasFlowCons_drcc = [
    (gamma-1) .* gasFlowMax / 2 <= gasFlow_sum <= (gamma+1) .* gasFlowMax / 2; %不改了试试
    repmat((gamma-1) .* gasFlowMax / 2,[1,nGasType]) <= gasFlow <= repmat((gamma+1) .* gasFlowMax / 2, [1,nGasType]);

    Prsmin <= Prs <= Prsmax;
%     u_Prs.^2 + phi_Prs.^2 * sigma <= ( (Prsmax-Prsmin)/2 - v_Prs).^2 * epsilon;
%     -(u_Prs+v_Prs) <= Prs - (Prsmax+Prsmin)/2 <= u_Prs+v_Prs;
%     u_Prs >= 0;
%     0 <= v_Prs <= (Prsmax-Prsmin)/2; %

    PHI == gasFlow_sum.^2 ./ newC.^2 ;
%     gamma.*(Prs(FB) .* phi_Prs(FB) - Prs(TB) .* phi_Prs(TB)) == gasFlow_sum .* phi_gasFlowSum ./ newC;
%     gamma .* (phi_Prs(FB).^2 - phi_Prs(TB).^2) == phi_gasFlowSum.^2 ./ newC;

    sigma_PHI >= 0;
    ]:'gasFlowSOCcons';
% -------------------------------------------
% gas composition Taylor
PTGbus = mpc.ptg(:,1); 
CgsPTG = sparse(PTGbus, (1:nPTG)', 1, nGb, nPTG); % connection matrix
QptgInbusMethane = CgsPTG * Qptg(:,1) ; QptgInbusHydrogen = CgsPTG * Qptg(:,2);
phi_QptgInbusMethane = CgsPTG * phi_Qptg(:,1) ; phi_QptgInbusHydrogen = CgsPTG * phi_Qptg(:,2);
QptgInbusForAllGasComposition = [QptgInbusMethane,QptgInbusHydrogen];
QptgInbusForAllGasComposition_phi = [phi_QptgInbusMethane,phi_QptgInbusHydrogen];

for r = 1:nGasType
    nodalGasInjectionForEachComp(:,r) = nodalGasInjection_hge(PGs.*mpc.gasCompositionForGasSource(:,r),QptgInbusForAllGasComposition(:,r),...
            gamma, gasFlow(:,r),mpc,nGs,nGb,nGl);
%     phi_nodalGasInjectionForEachComp(:,r) = nodalGasInjection_hge(phi_PGs.*mpc.gasCompositionForGasSource(:,r),QptgInbusForAllGasComposition_phi(:,r),...
%         gamma, phi_gasFlow(:,r),mpc,nGs,nGb,nGl);
end

gasCompositionCons_drcc = [
    repmat(sum(nodalGasInjectionForEachComp,2),[1,nGasType]) .* gasComposition - nodalGasInjectionForEachComp == 0;
%     repmat(sum(phi_nodalGasInjectionForEachComp,2),[1,nGasType]) .* gasComposition ...
%         + phi_gasComposition .* repmat(sum(nodalGasInjectionForEachComp,2),[1,nGasType]) - phi_nodalGasInjectionForEachComp == 0;
%         phi_gasComposition .* repmat(sum(phi_nodalGasInjectionForEachComp,2),[1,nGasType]) == 0; %?
    (repmat((1+gamma),[1,nGasType]).*gasComposition(mpc.Gline(:,1),:) + repmat((1-gamma),[1,nGasType]).*gasComposition(mpc.Gline(:,2),:))/2 .* repmat(gasFlow_sum,[1,nGasType]) - gasFlow == 0;
% xxx 这里还缺一条
    sigma_x >= 0;
    sum(gasComposition,2) == 1;
%     sum(phi_gasComposition,2) == 0;
    ]:'gasCompositionCons';

% summarize all the cons
constraints = [
    nodalGasFlowBalanceCons_drcc;
    PTGcons_drcc;    
    electricityCons_drcc;
    electricityBalanceCons_drcc;
    gasSecurityCons_drcc;
    otherSecurityCons_drcc;
    gasSourceCons_drcc;
    gasDemandCons_drcc;
    GPPcons_drcc;
    gasFlowCons_drcc;
    gasCompositionCons_drcc;
    sigma_x == 0;
    sigma_Qd == 0;
    sigma_Qgpp == 0;
    sigma_Gf == 0;
    sigma_PHI == 0;
    ];
%% solve the problem
objfcn = obj_operatingCost(Pg,PGs,LCg,Qptg, mpc,CDF) ...
    +  1*alpha_PHI * sum(sum(sigma_PHI)) + 10000* alpha_x * sum(sum(sigma_x)) ...
    + alpha_Qd * 100 * (sum(sum(sigma_Qd)) + sum(sum(sigma_Qgpp)) + sum(sum(sigma_Gf)));
% 要get dual，必须先设置好，而且把约束用cone函数写
options = sdpsettings('verbose',2,'solver','gurobi', 'debug',1);
% options.ipopt.tol = 1e-4;
options.ipopt.max_iter = 30000;
information = optimize(constraints, objfcn, options);
output{v} = information;
solverTime(v) = output{v}.solvertime;
%% results
Prs = value(Prs);
PGs = value(PGs); % Mm3/day
Qd = value(Qd); % Mm3/day
Qptg = value(Qptg);
Pptg = value(Pptg); % MW
Pg = value(Pg);  % MW
Pg_wind = value(Pg_wind);
Pgpp = value(Pgpp);  % MW
Qgpp = value(Qgpp);  % Mm3/day
Va = value(Va);
LCg = value(LCg);
gamma = value(gamma);
gasComposition = value(gasComposition); % hydrogen, gas
gasFlow = value(gasFlow);
gasFlow_sum = sum(value(gasFlow_sum),2);
S_nodal = value(S_nodal);
GCV_nodal = value(GCV_nodal)/1e6; % MJ/m3
WInodal = value(WI_nodal)/1e6; % MJ/m3
FSnodal = value(FSnodal);
PHI = value(PHI);
sigma_PHI = value(sigma_PHI);
sigma_x = value(sigma_x);
sigma_Qd = value(sigma_Qd);
sigma_Qgpp = value(sigma_Qgpp);
sigma_Gf = value(sigma_Gf);
% drcc vars
u_PGs = value(u_PGs);               v_PGs = value(v_PGs);               phi_PGs = value(phi_PGs); 
phi_Qd = value(phi_Qd);             phi_gasComposition = value(phi_gasComposition);
phi_Qgpp = value(phi_Qgpp);         phi_Qptg = value(phi_Qptg);         
u_Pptg = value(u_Pptg);             v_Pptg = value(v_Pptg);             phi_Pptg = value(phi_Pptg);
phi_gasFlow = value(phi_gasFlow);
phi_Pg = value(phi_Pg);             u_Pg = value(u_Pg);                 v_Pg = value(v_Pg);
phi_Va = value(phi_Va);             u_Va = value(u_Va);                 v_Va = value(v_Va); 
phi_Prs = value(phi_Prs);           u_Prs = value(u_Prs);               v_Prs = value(v_Prs);
Pg_nonwind = value(Pg_nonwind);
u_Pg_nonwind = value(u_Pg_nonwind);
v_Pg_nonwind = value(v_Pg_nonwind);
phi_Pg_nonwind = value(phi_Pg_nonwind);
phi_Pg_wind = value(phi_Pg_wind);
phi_windCapacity = value(phi_windCapacity);

[objfcn,totalCost,genAndLCeCost,gasPurchasingCost,gasCurtailmentCost,PTGsubsidy] = ...
    obj_operatingCost(Pg,PGs,LCg,Qptg, mpc, CDF);

[   sol{v}.totalCost,   sol{v}.genAndLCeCost,   sol{v}.gasPurchasingCost,   sol{v}.gasCurtailmentCost, ...
    sol{v}.PTGsubsidy,  sol{v}.objfcn,          sol{v}.Prs_square,          sol{v}.Prs,...
    sol{v}.PGs,         sol{v}.Qd,              sol{v}.Qptg,                sol{v}.Pptg,...
    sol{v}.Pg,          sol{v}.Pgpp,            sol{v}.Qgpp,                sol{v}.Va,...
    sol{v}.LCg,         sol{v}.gamma,           sol{v}.gasComposition,      sol{v}.gasFlow,...
    sol{v}.gasFlow_sum, sol{v}.S_nodal,         sol{v}.GCV_nodal,           sol{v}.WI,...
    sol{v}.PHI,         sol{v}.sigma_PHI,       sol{v}.sigma_x,             sol{v}.sigma_Qd,...
    sol{v}.sigma_Qgpp,  sol{v}.FSnodal,         ] = ...
deal(...
    totalCost,          genAndLCeCost,          gasPurchasingCost,          gasCurtailmentCost,...
    PTGsubsidy,         objfcn,                 Prs_square,                 Prs,...
    PGs,                Qd,                     Qptg,                       Pptg,...
    Pg,                 Pgpp,                   Qgpp,                       Va,...
    LCg,                gamma,                  gasComposition,             gasFlow,...
    gasFlow_sum,        S_nodal,                GCV_nodal,                  WInodal,...
    PHI,                sigma_PHI,              sigma_x,                    sigma_Qd,...
    sigma_Qgpp,         FSnodal);

%% covergence criterion
if v > 1
    criterion.sigma_PHI(v-1) = sum(sum(sigma_PHI))/1e2;
    criterion.sigma_x(v-1) = sum(sum(sigma_x));
    criterion.sigma_Qd(v-1) = sum(sum(sigma_Qd));
    criterion.sigma_Qgpp(v-1) = sum(sum(sigma_Qgpp));
    criterion.delta_sigma_PHI(v-1) = abs( (sum(sum(sigma_PHI)) - sum(sum(sol{v-1}.sigma_PHI)))) / abs( (sum(sum(sigma_PHI)) + sum(sum(sol{v-1}.sigma_PHI))));
    criterion.delta_sigma_x(v-1) = abs( (sum(sum(sigma_x)) - sum(sum(sol{v-1}.sigma_x)))) / abs( (sum(sum(sigma_x)) + sum(sum(sol{v-1}.sigma_x))));
    criterion.delta_sigma_Qd(v-1) = abs( (sum(sum(sigma_Qd)) - sum(sum(sol{v-1}.sigma_Qd)))) / abs( (sum(sum(sigma_Qd)) + sum(sum(sol{v-1}.sigma_Qd))));
    criterion.delta_sigma_Qgpp(v-1) = abs( (sum(sum(sigma_Qgpp)) - sum(sum(sol{v-1}.sigma_Qgpp)))) / abs( (sum(sum(sigma_Qgpp)) + sum(sum(sol{v-1}.sigma_Qgpp))));
    %cost
    criterion.totalCost(v-1) = abs( (totalCost - sol{v-1}.totalCost)) / sol{v-1}.totalCost;
    criterion.genAndLCeCost(v-1) = abs( (genAndLCeCost - sol{v-1}.genAndLCeCost)/sol{v-1}.genAndLCeCost );
    criterion.gasPurchasingCost(v-1) = abs( (gasPurchasingCost - sol{v-1}.gasPurchasingCost)/sol{v-1}.gasPurchasingCost );
    criterion.PTGsubsidy(v-1) = abs( (PTGsubsidy - sol{v-1}.PTGsubsidy)/sol{v-1}.PTGsubsidy );
    criterion.gasCompositionCH4(v-1) = abs(sum(gasComposition(:,1)) - sum(sol{v-1}.gasComposition(:,1)))*2 ./ (sum(gasComposition(:,1)) + sum(sol{v-1}.gasComposition(:,1)));
    criterion.gasCompositionH2(v-1) = abs(sum(gasComposition(:,2)) - sum(sol{v-1}.gasComposition(:,2)))*2 / (sum(gasComposition(:,2)) + sum(sol{v-1}.gasComposition(:,2)));
    criterion.S_nodal(v-1) = max(abs( (S_nodal-sol{v-1}.S_nodal)./sol{v-1}.S_nodal ));
    criterion.Qd(v-1) = abs( ( sum(sum(sigma_Qd)) - sum(sum(sol{v-1}.sigma_Qd)) ) / sum(sum(sol{v-1}.sigma_Qd)) );

    if  ((criterion.sigma_PHI(v-1)<=1e-3) && (criterion.sigma_x(v-1)<=1e-2) && (criterion.sigma_Qd(v-1)<=1e-2) && (criterion.sigma_Qgpp(v-1)<=1e-2) ) && ...
                ( (criterion.sigma_x(v-1)<=gap) || (criterion.delta_sigma_x(v-1)<=gap) ) ...
                && ( (criterion.sigma_Qd(v-1)<=gap) ||  (criterion.delta_sigma_Qd(v-1)<gap) ) ...
                && ( (criterion.sigma_Qgpp(v-1)<=gap) || (criterion.delta_sigma_Qgpp(v-1)<gap) )
            break
        else
        alpha_PHI = min([multiplier * alpha_PHI,alpha_PHI_max]);
        alpha_x = min([multiplier * alpha_x,alpha_x_max]);
        alpha_Qd = min([multiplier * alpha_Qd,alpha_Qd_max]);
    end
end
% ob
% ob.Pg = [ob.Pg, Pg];
%% update S and Z
R_nodal = gasComposition * Rall';
R_pipeline_ref = ( (1+gamma).*R_nodal(FB) + (1-gamma).*R_nodal(TB) ) / 2;
QptgInbus = CgsPTG * sum(Qptg,2);
W_ref = nodalGasInjection_hge(PGs,QptgInbus,gamma,gasFlow_sum,mpc,nGs,nGb,nGl);
gasComposition_ref = gasComposition;
gasFlow_sum_ref = gasFlow_sum;
if output{v}.problem ~= 0
%     error('optimization failed');
%     break
end
end


HGEinfo = output{v}; HGEsolution = sol{v};
end


