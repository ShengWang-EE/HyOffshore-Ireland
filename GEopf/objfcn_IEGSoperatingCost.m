function [objfcn,electricityGenerationCost,gasPurchasingCost] = objfcn_IEGSoperatingCost(Pg,PGs,mpc)
%% 
electricityGenerationCost = sum( Pg' * mpc.gencost(:,6));
gasPurchasingCost = sum(PGs' * mpc.Gcost);

objfcn = electricityGenerationCost + gasPurchasingCost;


end