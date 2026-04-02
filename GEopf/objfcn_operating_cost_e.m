function [objfcn] = objfcn_operating_cost_e(Pg,mpc)
Pg = mpc.baseMVA*Pg;
%% 
C = mpc.gencost(:,6);
electricityGenerationCost = sum(C.*Pg);

objfcn = electricityGenerationCost;


end