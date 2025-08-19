function [objfcn,electricityGenerationCost,gasPurchasingCost] = objfcn_IEGSoperatingCost(Pg,PGs,mpc)
Pg = mpc.baseMVA*Pg;
%% 
electricityGenerationCost = sum( Pg' * mpc.gencost(:,6));
gasPurchasingCost = sum(PGs' * mpc.Gcost);

objfcn = electricityGenerationCost + gasPurchasingCost;


end