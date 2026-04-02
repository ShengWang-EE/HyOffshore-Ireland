function [g] = consfcn_gasNodalBalance(Gf,Pg,PGs,Qgpp,LCg,Qptg,mpc,iGd,nGb,nGl,nGs,nGd)
%% initialization
nGPP = size(mpc.GEcon(:,3),1);
nPTG = size(mpc.ptg,1);
%%
% 1 supply-demand
PGsbus = mpc.Gsou(:,1) ; 
Cgs_PGs = sparse(PGsbus, (1:nGs)', 1, nGb, nGs); % connection matrix
Qdbus = iGd; 
Cgs_Qd = sparse(Qdbus, (1:nGd)', 1, nGb, nGd); % connection matrix
g = Cgs_PGs*PGs - Cgs_Qd * mpc.Gbus(iGd,3); % supply-demand, Mm3/day
%2 ptg (因为GB里面有重复的，所以得一个一个累加，不能用向量累加）
for i = 1:nPTG
    GB = mpc.ptg(i,1);
    g(GB) = g(GB) + Qptg(i);
end
% 3 gpp
for i = 1:nGPP
    GB = mpc.GEcon(i,1);
    g(GB) = g(GB)-Qgpp(i);
end
% gas flow
for  m = 1:nGl
    fb = mpc.Gline(m,1); tb = mpc.Gline(m,2);
    g(fb) = g(fb) - Gf(m);
    g(tb) = g(tb) + Gf(m);
end
% LCg
g(iGd) = g(iGd) + LCg;

end

