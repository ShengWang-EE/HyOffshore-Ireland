function [solution, information] = runopf_ge_rbst(mpc,xi_hat)
% 20220810: 纯电力天然气的稳态opf，考虑燃气机组和ptg，考虑切负荷
%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
% create (read-only) copies of individual fields for convenience\
baseMVA=100;
il = find(mpc.branch(:, RATE_A) ~= 0 & mpc.branch(:, RATE_A) < 1e10);
% robust optimisation paras
robust_wind_reduction = max(abs(xi_hat));
%% initialization
nGb = size(mpc.Gbus,1);
nGs = size(mpc.Gsou,1);
nGl = size(mpc.Gline,1);
iGd = find(mpc.Gbus(:,3)~=0);
nGd = size(iGd,1);
nb = size(mpc.bus,1);
ng = size(mpc.gen,1);
id = find(mpc.bus(:,3)~=0);
nd = size(id,1);
nPTG = size(mpc.ptg,1);
i_owf = find(mpc.gentype == 1);
i_nonowf = find(mpc.gentype == 0);

refbus = find(mpc.bus(:, BUS_TYPE) == REF);

[GCV, M, fs, a, R, T_stp, Prs_stp, Z, T_gas, eta, CDF] = initializeParameters_J13();
[B, Bf, Pbusinj, Pfinj] = makeBdc(mpc.baseMVA, mpc.bus, mpc.branch);
upf = mpc.branch(il, RATE_A) - Pfinj(il);
%% state variable
Prs_square = sdpvar(nGb,1); % nodal gas pressure
PGs = sdpvar(nGs,1); % gas production of gas source
Gf = sdpvar(nGl,1); % gas flow in the pipeline
Qptg = sdpvar(nPTG,1); % gas (hydrogen) production of PTG

Va = sdpvar(nb,1); % voltage phase angle
Pg = sdpvar(ng,1); % electricity generation (including gas fired units)

% see if has pre-direction
if size(mpc.Gline,2) == 9 % has direction information
    gamma = mpc.Gline(:,9);
else
    alpha = binvar(nGl,1);          % direction of gas flow, 0-1
    gamma = (alpha-0.5)*2;          % 1,-1
    w = sdpvar(nGl,1);            % auxiliary variable for gas flow
end
%% upper and lower bounds
Prs_square_min = mpc.Gbus(:,5).^2; Prs_square_max = mpc.Gbus(:,6).^2;
PGsmin = mpc.Gsou(:,3); PGsmax = mpc.Gsou(:,4);
Gfmin = -mpc.Gline(:,5); Gfmax = mpc.Gline(:,5); 
Qptgmin = 0; Qptgmax = mpc.ptg(:,5);
Vamin = -Inf(nb,1); Vamax = Inf(nb,1); 
Vamin(refbus) = 1; Vamax(refbus) = 1;
Pgmin = mpc.gen(:, PMIN) *0; %Pgmin is set to zero
Pgmax = mpc.gen(:, PMAX);
%% contraints
Pptg = Qptg/24/3600 * GCV.ng * eta.electrolysis;
boxCons = [
    Prs_square_min <= Prs_square <= Prs_square_max;
    PGsmin <= PGs <= PGsmax;
    Qptgmin <= Qptg <= Qptgmax;
    Pgmin(i_nonowf) <= Pg(i_nonowf) <= Pgmax(i_nonowf);
    Pgmin(i_owf) <= Pg(i_owf) <= Pgmax(i_owf) - robust_wind_reduction';
    ];
electricityNodalBalanceConsDC = [consfcn_electricityNodalBalance(Va,Pg,Pptg,0,mpc,id) == 0;];
electricityBranchFlowConsDC = [- upf <= Bf(il,:)*Va <= upf;];
Pgpp = Pg(mpc.GEcon(:,3));
Qgpp = Pgpp / GCV.ng * 3600 * 24 * eta.GFU;
gasNodalBalanceCons = [consfcn_gasNodalBalance(Gf,Pg,PGs,Qgpp,0,Qptg,mpc,iGd,nGb,nGl,nGs,nGd) == 0;];

gasPipelineFlowCons = [
    (gamma-1) .* Gfmax / 2 <= Gf <= (gamma+1) .* Gfmax / 2;
    ];
% gas flow
FB = mpc.Gline(:,1); TB = mpc.Gline(:,2);
if size(mpc.Gline,2) == 9 % has direction information
    gasFlowCons = [
        Prs_square(FB) - Prs_square(TB) == gamma .* Gf.^2 ./ mpc.Gline(:,3).^2;
        ];
else
    gasFlowCons = [
        Prs_square(FB) - Prs_square(TB) == (2*w - Gf.^2) ./ mpc.Gline(:,3).^2;
        0 <= w <= Gf.^2;
        w <= alpha .* Gfmax.^2;
        w >= Gf.^2 - (1-alpha) .* Gfmax.^2;
        ];
end
%% solotion
cons = [
    boxCons;
    electricityNodalBalanceConsDC;
    electricityBranchFlowConsDC;
    gasNodalBalanceCons;
    gasPipelineFlowCons;
    gasFlowCons;
    ];
objfcn = objfcn_IEGSoperatingCost(Pg,PGs,mpc); % closer to mid level
yalmipOptions = sdpsettings('solver','gurobi');
% yalmipOptions.gurobi.MIPgapabs = 1e1;

information = optimize(cons,objfcn,yalmipOptions);
%% results
Prs = sqrt(value(Prs_square)); % nodal gas pressure
PGs = value(PGs); % gas production of gas source
Gf = value(Gf); % gas flow in the pipeline
Qptg = value(Qptg); % gas (hydrogen) production of PTG
Va = value(Va); % voltage phase angle
Pg = value(Pg); % electricity generation (including gas fired units)
gamma = value(gamma);
Qgpp = value(Qgpp);
[operatingCost,electricityGenerationCost,gasPurchasingCost] = objfcn_IEGSoperatingCost(Pg,PGs,mpc);

[solution.Prs, solution.PGs, solution.Gf, solution.Qptg, solution.Qgpp, ...
    solution.Va, solution.Pg, solution.gamma, ...
    solution.objfcn, solution.electricityGenerationCost,  ...
    solution.gasPurchasingCost] = deal(...
    Prs, PGs, Gf, Qptg, Qgpp, Va, Pg, gamma, ...
    operatingCost, electricityGenerationCost, gasPurchasingCost);
end
