function [solution, information] = runopf_ge_mmt_drcc(mpc,mu,sigma_sqrt,eps)
% 20251117 wasserstein based DRCC
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
il = find(mpc.branch(:, RATE_A) ~= 0 & mpc.branch(:, RATE_A) < 1e10);
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
iOWF = find(mpc.gentype == 1);
i_nonowf = find(mpc.gentype == 0);
nOWF = size(iOWF,1);
refbus = find(mpc.bus(:, BUS_TYPE) == REF);

[GCV, M, fs, a, R, T_stp, Prs_stp, Z, T_gas, eta, CDF] = initializeParameters_J13();
%% state variable
Prs = sdpvar(nGb,1); % nodal gas pressure
PGs = sdpvar(nGs,1); % gas production of gas source
Gf = sdpvar(nGl,1); % gas flow in the pipeline
Qptg = sdpvar(nPTG,1); % gas (hydrogen) production of PTG

Va = sdpvar(nb,1); % voltage phase angle
Pg = sdpvar(ng,1); % electricity generation (including gas fired units)
% adjust factor
alpha_Prs = sdpvar(nGb,1);
alpha_PGs = sdpvar(nGs,1);
alpha_Gf = sdpvar(nGl,1); 
alpha_Qptg = sdpvar(nPTG,1);
alpha_Va = sdpvar(nb,1);
alpha_Pg = sdpvar(ng,1);
% see if has pre-direction
if size(mpc.Gline,2) == 9 % has direction information
    gamma = mpc.Gline(:,9);
else
    % alpha = binvar(nGl,1);          % direction of gas flow, 0-1
    gamma = (alpha-0.5)*2;          % 1,-1
    w = sdpvar(nGl,1);            % auxiliary variable for gas flow
end
%% upper and lower bounds
Prsmin = mpc.Gbus(:,5); Prsmax = mpc.Gbus(:,6);
PGsmin = mpc.Gsou(:,3); PGsmax = mpc.Gsou(:,4);
Gfmin = -mpc.Gline(:,5); Gfmax = mpc.Gline(:,5); 
Qptgmin = 0; Qptgmax = mpc.ptg(:,5);
Vamin = -Inf(nb,1); Vamax = Inf(nb,1); 
Vamin(refbus) = 1; Vamax(refbus) = 1;
Pgmin = mpc.gen(:, PMIN) *0; %Pgmin is set to zero
Pgmax = mpc.gen(:, PMAX);
%% contraints
[B, Bf, Pbusinj, Pfinj] = makeBdc(mpc.baseMVA, mpc.bus, mpc.branch);
Pptg = Qptg/24/3600 * GCV.ng * eta.electrolysis;
Pgpp = Pg(mpc.GEcon(:,3));
Qgpp = Pgpp / GCV.ng * 3600 * 24 * eta.GFU;
upf = mpc.branch(il, RATE_A) - Pfinj(il);

alpha_Pptg = alpha_Qptg/24/3600 * GCV.ng * eta.electrolysis;
alpha_Pgpp = alpha_Pg(mpc.GEcon(:,3));
alpha_Qgpp = alpha_Pgpp / GCV.ng * 3600 * 24 * eta.GFU;
% normal constraints
boxCons = [
    Prsmin <= Prs <= Prsmax;
    PGsmin <= PGs <= PGsmax;
    Qptgmin <= Qptg <= Qptgmax;
    Pgmin <= Pg <= Pgmax;
    ];
electricityBranchFlowConsDC = [- upf <= Bf(il,:)*Va <= upf;];

electricityNodalBalanceConsDC = [consfcn_electricityNodalBalance(Va,Pg,Pptg,0,mpc,id) == 0;];

gasNodalBalanceCons = [consfcn_gasNodalBalance(Gf,Pg,PGs,Qgpp,0,Qptg,mpc,iGd,nGb,nGl,nGs,nGd) == 0;];

gasPipelineFlowCons = [
    (gamma-1) .* Gfmax / 2 <= Gf <= (gamma+1) .* Gfmax / 2;
    ];
% gas flow
FB = mpc.Gline(:,1); TB = mpc.Gline(:,2);
if size(mpc.Gline,2) == 9 % has direction information
    gasFlowCons = [
        Prs(FB).^2 - Prs(TB).^2 == gamma .* Gf.^2 ./ mpc.Gline(:,3).^2;
        ];
else
    gasFlowCons = [
        Prs_square(FB) - Prs_square(TB) == (2*w - Gf.^2) ./ mpc.Gline(:,3).^2;
        0 <= w <= Gf.^2;
        w <= alpha .* Gfmax.^2;
        w >= Gf.^2 - (1-alpha) .* Gfmax.^2;
        ];
end

% moment based DRCC constraints
factor1 = norm(sum(sigma_sqrt,2)); % 推导过程见onenote，适用于一般的根据所有误差累加的ajust factor;推导出这个系数的目的是让约束可以用向量方式写
for i = 1:nOWF
    e_i = zeros(nOWF,1); e_i(i) = 1;
    factor2(i,1) = norm(sigma_sqrt * (alpha_Pg(iOWF)-e_i));
end
boxCons_drcc = [
    % Prsmin + sqrt((1-eps)/eps) * factor1*sqrt(alpha_Prs.^2) <= Prs <= Prsmax - sqrt((1-eps)/eps) * factor1*sqrt(alpha_Prs.^2);
    Prsmin <= Prs <= Prsmax;
    PGsmin + sqrt((1-eps)/eps) * factor1*sqrt(alpha_PGs.^2) <= PGs <= PGsmax - sqrt((1-eps)/eps) * factor1*sqrt(alpha_PGs.^2);
    Qptgmin + sqrt((1-eps)/eps) * factor1*sqrt(alpha_Qptg.^2) <= Qptg <= Qptgmax - sqrt((1-eps)/eps) * factor1*sqrt(alpha_Qptg.^2);
    Pgmin(iOWF) + sqrt((1-eps)/eps) * factor2 <= Pg(iOWF) <= Pgmax(iOWF) - sqrt((1-eps)/eps) * factor2; % for owf
    Pgmin(i_nonowf) + sqrt((1-eps)/eps) * factor1*sqrt(alpha_Pg(i_nonowf).^2) <= Pg(i_nonowf) ...
        <= Pgmax(i_nonowf) - sqrt((1-eps)/eps) * factor1*sqrt(alpha_Pg(i_nonowf).^2);
    ];

electricityBranchFlowConsDC_drcc = [- upf <= Bf(il,:)*Va <= upf;];
gasPipelineFlowCons_drcc = [
    (gamma-1) .* Gfmax / 2 <= Gf <= (gamma+1) .* Gfmax / 2;
    ];

electricityNodalBalanceConsDC_drcc = [consfcn_electricPowerBalance_hge_drcc(Va,Pg,Pptg,alpha_Va, alpha_Pg, alpha_Pptg, mpc) == 0;];

gasNodalBalanceCons_drcc = [consfcn_gasNodalBalance_drcc(Gf,PGs,Qgpp,Qptg,alpha_Gf,alpha_PGs,alpha_Qgpp,alpha_Qptg,mpc,iGd,nGb,nGl,nGs,nGd) == 0;];

if size(mpc.Gline,2) == 9 % has direction information
    C = mpc.Gline(:,3)*1;
    gasFlowCons_drcc = [
        Prs(FB).^2 - Prs(TB).^2 == gamma./ C .^2 .* Gf .^2 ;
        alpha_Prs(FB) - alpha_Prs(TB) == gamma./ C .^2 .* alpha_Gf;
        alpha_Prs(FB).^2 - alpha_Prs(TB).^2 == alpha_Gf.^2;
        ];
end

%% solotion
cons = [
    % boxCons;
    % electricityNodalBalanceConsDC;
    % electricityBranchFlowConsDC;
    % gasNodalBalanceCons;
    % gasPipelineFlowCons;
    % gasFlowCons;
        % Pgmin <= Pg <= Pgmax;
    % PGsmin <= PGs <= PGsmax;
    boxCons_drcc;
    electricityNodalBalanceConsDC_drcc;
    electricityBranchFlowConsDC_drcc;
    gasNodalBalanceCons_drcc;
    gasPipelineFlowCons_drcc;
    gasFlowCons_drcc;
    ];
objfcn = objfcn_IEGSoperatingCost(Pg,PGs,mpc); % closer to mid level
yalmipOptions = sdpsettings('solver','gurobi');
% yalmipOptions.gurobi.MIPgap = 0.1;
% yalmipOptions.gurobi.MIPgapabs = 1e1;
% yalmipOptions.ipopt.max_iter = 1e4;
% yalmipOptions.relax = 1;

information = optimize(cons,objfcn,yalmipOptions);
%% results
Prs = value(Prs); % nodal gas pressure
PGs = value(PGs); % gas production of gas source
Gf = value(Gf); % gas flow in the pipeline
Qptg = value(Qptg); % gas (hydrogen) production of PTG
Va = value(Va); % voltage phase angle
Pg = value(Pg); % electricity generation (including gas fired units)
gamma = value(gamma);
Qgpp = value(Qgpp);
alpha_Prs = value(alpha_Prs);
alpha_PGs = value(alpha_PGs);
alpha_Gf = value(alpha_Gf);
alpha_Qptg = value(alpha_Qptg);
alpha_Va = value(alpha_Va);
alpha_Pg = value(alpha_Pg);

[operatingCost,electricityGenerationCost,gasPurchasingCost] = objfcn_IEGSoperatingCost(Pg,PGs,mpc);

[solution.Prs, solution.PGs, solution.Gf, solution.Qptg, solution.Qgpp, ...
    solution.Va, solution.Pg, solution.gamma, ...
    solution.objfcn, solution.electricityGenerationCost,  ...
    solution.gasPurchasingCost] = deal(...
    Prs, PGs, Gf, Qptg, Qgpp, Va, Pg, gamma, ...
    operatingCost, electricityGenerationCost, gasPurchasingCost);
end
