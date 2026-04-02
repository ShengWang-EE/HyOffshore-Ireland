function [solution, information] = runopf_e(mpc)
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

refbus = find(mpc.bus(:, BUS_TYPE) == REF);
%% state variable
Va = sdpvar(nb,1); % voltage phase angle
Pg = sdpvar(ng,1); % electricity generation (including gas fired units)

% set as zero
Pptg = 0;
%% upper and lower bounds
Vamin = -Inf(nb,1); Vamax = Inf(nb,1); 
Vamin(refbus) = 1; Vamax(refbus) = 1;
Pgmin = mpc.gen(:, PMIN) / baseMVA *0; %Pgmin is set to zero
Pgmax = mpc.gen(:, PMAX) / baseMVA;
%% contraints
boxCons = [
    Pgmin <= Pg <= Pgmax;
    ];
electricityNodalBalanceConsDC = [consfcn_electricityNodalBalance(Va,Pg,Pptg,0,mpc,id) == 0;];
electricityBranchFlowConsDC = [consfcn_electricityBranchFlow(Va, mpc, il) <= 0;];

%% solotion
cons = [
    boxCons;
    electricityNodalBalanceConsDC;
    electricityBranchFlowConsDC;
    ];
objfcn = objfcn_operating_cost_e(Pg,mpc); % closer to mid level
yalmipOptions = sdpsettings('verbose',2,'solver','gurobi','debug',1);

information = optimize(cons,objfcn,yalmipOptions);
%% results
Va = value(Va); % voltage phase angle
Pg = value(Pg); % electricity generation (including gas fired units)
objfcn = value(objfcn);

[solution.Va, solution.Pg, solution.objfcn] = deal(Va, Pg, objfcn);
end
