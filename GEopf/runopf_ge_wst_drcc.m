function [solution, information] = runopf_ge_wst_drcc(mpc,xi_hat,eps)
% 20251117 moment based DRCC
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
% create (read-only) copies of individual fields for convenience
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
[n_s,n_xi] = size(xi_hat);
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
% other auxilary variables
tau_Prs = sdpvar(nGb,1);        lambda_Prs = sdpvar(nGb,1);     sigma_Prs = sdpvar(nGb,n_s);
tau_PGs = sdpvar(nGs,1);        lambda_PGs = sdpvar(nGs,1);     sigma_PGs = sdpvar(nGs,n_s);
tau_Gf = sdpvar(nGl,1);         lambda_Gf = sdpvar(nGl,1);      sigma_Gf = sdpvar(nGl,n_s);
tau_Qptg = sdpvar(nPTG,1);        lambda_Qptg = sdpvar(nPTG,1);     sigma_Qptg = sdpvar(nPTG,n_s);
tau_Va = sdpvar(nb,1);        lambda_Va = sdpvar(nb,1);     sigma_Va = sdpvar(nb,n_s);
tau_Pg = sdpvar(ng,1);        lambda_Pg = sdpvar(ng,1);     sigma_Pg = sdpvar(ng,n_s);

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

rho = 1/sqrt(n_s) * 2;
% eps = 0.2;
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

% wst based DRCC constraints
% 1) Prs
% cons_Prs = [];
% for i = 1:nGb % 假设一条一条约束处理
%     a = alpha_Prs(i) * ones(n_xi,1);
%     b = Prsmax(i) - Prs(i);
%     cons1 = [tau_Prs(i) + 1/eps * (lambda_Prs(i)*rho + 1/n_s * sum(sigma_Prs(i,:)) ) <= 0;];
%     cons2 = [];
%     for i_s = 1:n_s % for all samples
%         cons2 = [cons2;
%                  a'*xi_hat(i_s,:)' - b - tau_Prs(i) <= sigma_Prs(i,i_s)];
%     end
%     cons3 = [norm(a,2) <= lambda_Prs(i)];
%     cons_Prs = [cons_Prs;
%                 cons1;
%                 cons2;
%                 cons3;];
% end
% cons_Prs = [
%     cons_Prs;
%     sigma_Prs >= 0;
%     lambda_Prs >= 0;
%     ];

% 1) Box constraints
% 1.1) Prs cons
cons_Prs = [
    tau_Prs + 1/eps * ( lambda_Prs * rho + 1/n_s * sum(sigma_Prs, 2) ) <= 0;
    alpha_Prs * sum(xi_hat,2).' - (Prsmax-Prs) * ones(1,n_s) - tau_Prs*ones(1,n_s) <= sigma_Prs;
    lambda_Prs >= sqrt(n_xi) * alpha_Prs;
    lambda_Prs >= - sqrt(n_xi) * alpha_Prs;
    sigma_Prs >= 0;
    lambda_Prs >= 0;
    Prsmin <= Prs <= Prsmax;
    ]; % 认为天然气系统不随之调节
% 1.2) PGs cons
cons_PGs = [
    % % upper bound
    tau_PGs + 1/eps * ( lambda_PGs * rho + 1/n_s * sum(sigma_PGs, 2) ) <= 0;
    alpha_PGs * sum(xi_hat,2).' - (PGsmax-PGs) * ones(1,n_s) - tau_PGs*ones(1,n_s) <= sigma_PGs;
    lambda_PGs >= sqrt(n_xi) * alpha_PGs;
    lambda_PGs >= - sqrt(n_xi) * alpha_PGs;
    sigma_PGs >= 0;
    lambda_PGs >= 0;
    % % lower bound
    % PGs >= PGsmin;
    PGsmin <= PGs <= PGsmax;
    ];
% 1.3) Qptg cons
cons_Qptg = [
    % % upper bound
    tau_Qptg + 1/eps * ( lambda_Qptg * rho + 1/n_s * sum(sigma_Qptg, 2) ) <= 0;
    alpha_Qptg * sum(xi_hat,2).' - (Qptgmax-Qptg) * ones(1,n_s) - tau_Qptg*ones(1,n_s) <= sigma_Qptg;
    lambda_Qptg >= sqrt(n_xi) * alpha_Qptg;
    lambda_Qptg >= - sqrt(n_xi) * alpha_Qptg;
    sigma_Qptg >= 0;
    lambda_Qptg >= 0;
    % % lower bound
    % Qptg >= Qptgmin;
    Qptgmin <= Qptg <= Qptgmax;
    ];
% 1.4) nonowf cons
tau_nonowf = tau_Pg(i_nonowf,:);        lambda_nonowf = lambda_Pg(i_nonowf,:);
sigma_nonowf = sigma_Pg(i_nonowf,:);    alpha_nonowf = alpha_Pg(i_nonowf,:);
Pg_nonowfmax = Pgmax(i_nonowf,:);       Pg_nonowfmin = Pgmin(i_nonowf,:);
Pg_nonowf = Pg(i_nonowf,:);
cons_nonowf = [
    % upper bound
    tau_nonowf + 1/eps * ( lambda_nonowf * rho + 1/n_s * sum(sigma_nonowf, 2) ) <= 0;
    alpha_nonowf * sum(xi_hat,2).' - (Pg_nonowfmax-Pg_nonowf) * ones(1,n_s) - tau_nonowf*ones(1,n_s) <= sigma_nonowf;
    lambda_nonowf >= sqrt(n_xi) * alpha_nonowf;
    lambda_nonowf >= - sqrt(n_xi) * alpha_nonowf;
    sigma_nonowf >= 0;
    lambda_nonowf >= 0;
    % lower bound
    Pg_nonowf >= Pg_nonowfmin;
    ];

% 1.5) owf cons
cons_owf = [];
for i = 1:nOWF % 假设一条一条约束处理
    ii = iOWF(i);
    e_i = zeros(n_xi,1); e_i(i) = 1;
    a = alpha_Pg(ii) * ones(n_xi,1) - e_i;
    b = Pgmax(ii) - Pg(ii);
    cons1 = [tau_Pg(ii) + 1/eps * (lambda_Pg(ii)*rho + 1/n_s * sum(sigma_Pg(ii,:)) ) <= 0;];
    cons2 = [];
    for i_s = 1:n_s % for all samples
        cons2 = [cons2;
                 a'*xi_hat(i_s,:)' - b - tau_Pg(ii) <= sigma_Pg(ii,i_s)];
    end
    cons3 = [norm(a,2) <= lambda_Pg(ii)];
    cons_owf = [cons_owf;
                cons1;
                cons2;
                cons3;];
end
cons_owf = [
    cons_owf;
    sigma_Pg >= 0;
    lambda_Pg >= 0;
    ];
% test
% cons_owf = [
%     Pgmin(iOWF) <= Pg(iOWF) <= Pgmax(iOWF);
%     ];

boxCons_drcc = [
    % cons_Prs;
    % cons_PGs;
    % cons_Qptg;
    cons_nonowf;
    cons_owf;
    ];
% 2) other cons
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
cons_alpha = [
    0 <= [alpha_Prs;alpha_PGs;alpha_Gf;alpha_Qptg;alpha_Va;alpha_Pg] <= 0.5;
    ];

%% solotion
cons = [
    boxCons;
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
    % cons_alpha;
    ];
objfcn = objfcn_IEGSoperatingCost(Pg,PGs,mpc); % closer to mid level
yalmipOptions = sdpsettings('solver','gurobi');
% yalmipOptions.gurobi.MIPgap = 0.1;
% yalmipOptions.gurobi.MIPgapabs = 1e1;
yalmipOptions.ipopt.max_iter = 1e4;
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
alpha_Pf = 0;

[operatingCost,electricityGenerationCost,gasPurchasingCost] = objfcn_IEGSoperatingCost(Pg,PGs,mpc);

[solution.Prs, solution.PGs, solution.Gf, solution.Qptg, solution.Qgpp, ...
    solution.Va, solution.Pg, solution.gamma, ...
    solution.objfcn, solution.electricityGenerationCost, solution.gasPurchasingCost, ...
    solution.alpha_Prs, solution.alpha_PGs, solution.alpha_Gf, solution.alpha_Qptg, solution.alpha_Va, solution.alpha_Pg, ...   
    solution.alpha_Pf] = deal(...
    Prs, PGs, Gf, Qptg, Qgpp, Va, Pg, gamma, ...
    operatingCost, electricityGenerationCost, gasPurchasingCost, ...
    alpha_Prs, alpha_PGs, alpha_Gf, alpha_Qptg, alpha_Va, alpha_Pg, alpha_Pf);
end
