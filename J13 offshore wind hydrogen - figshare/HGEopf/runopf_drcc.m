function [HGEsolution, HGEinfo] = runopf_drcc(mpc,mu,sigma,epsilon)
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

%%
solverTime = zeros(iterationMax,1);
for v = 1:iterationMax
    yalmip('clear')
%% state variables
Prs = sdpvar(nGb,1); % bar^2
PGs = sdpvar(nGs,1); % Mm3/day
Qd = sdpvar(nGd,nGasType);% Mm3/day
Qptg = sdpvar(nPTG,2); % [ methane; hydrogen ] % Mm3/day
Pptg = zeros(nPTG,1); % electricity consumption, 1/100 MW
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
phi_PGs = sdpvar(nGs,1); % affine mapping factor
phi_Qd = sdpvar(nGd,nGasType);
phi_gasComposition = sdpvar(nGb,nGasType);
phi_Qgpp = sdpvar(nGpp, nGasType);
phi_Qptg = sdpvar(nPTG,2);
phi_Pptg = sdpvar(nPTG,1);
u_Pptg = sdpvar(nPTG,1); v_Pptg = sdpvar(nPTG,1); 
phi_gasFlow = sdpvar(nGl,nGasType);
phi_Pg = sdpvar(ng,1); u_Pg = sdpvar(ng,1); v_Pg = sdpvar(ng,1);
phi_Va = sdpvar(nb,1); u_Va = sdpvar(nl,1); v_Va = sdpvar(nl,1); 
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
newC = mpc.Gline(:,3);
% newC = mpc.Gline(:,3);
%% constraints
sigma = sigma/baseMVA;
mu = mu/baseMVA;
% electricity flow
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

    Pg_nonwind + sqrt((1-epsilon)/epsilon) * sigma * phi_Pg_nonwind .^2 <= Pgmax_nonwind;
    Pg_nonwind >= 0;

    sum(phi_windCapacity) == 1;
    % 
    Pg_wind + sqrt((1-epsilon)/epsilon) * sigma * phi_windCapacity .^2 <= Pgmax_wind;
    Pg_wind >= 0;

    Va(refs) == 0; 
    ]:'electricityCons';

electricityBalanceCons = [...
    consfcn_electricPowerBalance_hge(Va,Pg,Pptg,mpc) == 0;
    ]:'electricityBalanceCons';

% summarize all the cons
constraints = [
    electricityCons_drcc;
    electricityBalanceCons;
    ];
%% solve the problem
objfcn = obj_operatingCost(Pg,0,0,0, mpc,CDF);

options = sdpsettings('verbose',2,'solver','mosek', 'debug',1);
% options.ipopt.tol = 1e-4;
options.ipopt.max_iter = 30000;
% options.gurobi.IterationLimit = 100000;
% options.gurobi.TuneTimeLimit = 0;
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
Qgpp = value(Qgpp);  % Mm3/day
Va = value(Va);
LCg = value(LCg);
gamma = value(gamma);
gasComposition = value(gasComposition); % hydrogen, gas
gasFlow = value(gasFlow);


[objfcn,totalCost,genAndLCeCost,gasPurchasingCost,gasCurtailmentCost,PTGsubsidy] = ...
    obj_operatingCost(Pg,0,0,0, mpc, CDF);

[   sol{v}.totalCost,   sol{v}.genAndLCeCost,   sol{v}.gasPurchasingCost,   sol{v}.gasCurtailmentCost, ...
    sol{v}.PTGsubsidy,  sol{v}.objfcn,          ...
    sol{v}.Pg,                                      sol{v}.Va...
            ] = ...
deal(...
    totalCost,          genAndLCeCost,          gasPurchasingCost,          gasCurtailmentCost,...
    PTGsubsidy,         objfcn,               ... 
    Pg * baseMVA,                                       Va...
    );




HGEinfo = output{v}; HGEsolution = sol{v};
end


