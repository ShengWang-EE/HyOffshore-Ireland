function [information, results] = runopf_gas(mpc)
%% grab dimenssions and parameters
nGb = size(mpc.Gbus,1);
nGs = size(mpc.Gsou,1);
nGl = size(mpc.Gline,1);
iGd = find(mpc.Gbus(:,3)~=0);
nGd = size(iGd,1);
%% state variable
Prs_square = sdpvar(nGb,1);     % nodal gas pressure
PGs = sdpvar(nGs,1);            % gas production of gas source
Gf = sdpvar(nGl,1);             % gas flow in the pipeline
alpha = binvar(nGl,1);          % direction of gas flow, 0-1
gamma = (alpha-0.5)*2;          % 1,-1
w = sdpvar(nGl,1);            % auxiliary variable for gas flow
%% upper and lower bounds
Prs_square_min = mpc.Gbus(:,5).^2; Prs_square_max = mpc.Gbus(:,6).^2;
PGsmin = mpc.Gsou(:,3); PGsmax = mpc.Gsou(:,4);
Gfmin = -mpc.Gline(:,5); Gfmax = mpc.Gline(:,5); 
%% contraints
boxCons = [
    Prs_square_min <= Prs_square <= Prs_square_max;
    PGsmin <= PGs <= PGsmax;
    ];
Pg = 0; Qgpp = 0; LCg = 0; Qptg = 0;
gasNodalBalanceCons = [consfcn_gasNodalBalance_gasOnly(Gf,Pg,PGs,Qgpp,LCg,Qptg,mpc,iGd,nGb,nGl,nGs,nGd) == 0;];

gasPipelineFlowCons = [
    (gamma-1) .* Gfmax / 2 <= Gf <= (gamma+1) .* Gfmax / 2;
    ];
% gas flow
FB = mpc.Gline(:,1); TB = mpc.Gline(:,2);
mpc.Gline(:,3) = mpc.Gline(:,3)/1.475;
gasFlowCons = [
%     Prs_square(FB) - Prs_square(TB) == gamma .* Gf.^2 ./ mpc.Gline(:,3).^2;
    Prs_square(FB) - Prs_square(TB) == (2*w - Gf.^2) ./ mpc.Gline(:,3).^2;
    0 <= w <= Gf.^2;
    w <= alpha .* Gfmax.^2;
    w >= Gf.^2 - (1-alpha) .* Gfmax.^2;
    ];

%% solotion
cons = [
    boxCons;
    gasNodalBalanceCons;
    gasPipelineFlowCons;
    gasFlowCons;
    ];
objfcn = sum(PGs' * mpc.Gcost);
yalmipOptions = sdpsettings('verbose',2,'solver','gurobi','usex0',0,'debug',1);

information = optimize(cons,objfcn,yalmipOptions);
%% results
Prs = sqrt(value(Prs_square)); % nodal gas pressure
PGs = value(PGs); % gas production of gas source
Gf = value(Gf); % gas flow in the pipeline
gamma = value(gamma);
objfcn = value(objfcn);

results = mpc;
results.Gline = [results.Gline, Gf, gamma];
results.Gbus = [results.Gbus, Prs];
results.Gsou = [results.Gsou, PGs];
results.objfcn = objfcn;


end

function [g] = consfcn_gasNodalBalance_gasOnly(Gf,Pg,PGs,Qgpp,LCg,Qptg,mpc,iGd,nGb,nGl,nGs,nGd)
PGsbus = mpc.Gsou(:,1) ; 
Cgs_PGs = sparse(PGsbus, (1:nGs)', 1, nGb, nGs); % connection matrix
Qdbus = iGd; 
Cgs_Qd = sparse(Qdbus, (1:nGd)', 1, nGb, nGd); % connection matrix
g = Cgs_PGs*PGs - Cgs_Qd * mpc.Gbus(iGd,3); % supply-demand, Mm3/day

for  m = 1:nGl
    fb = mpc.Gline(m,1); tb = mpc.Gline(m,2);
    g(fb) = g(fb) - Gf(m);
    g(tb) = g(tb) + Gf(m);
end
% LCg
g(iGd) = g(iGd) + LCg;

end


