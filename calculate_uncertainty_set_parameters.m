function [forecastErrorSample,meanError,Sigma,forecastError_test] = calculate_uncertainty_set_parameters(sampledOWFgeneration, allOWFgeneration)
%% trainning set
[nk,nOWF] = size(sampledOWFgeneration);

predictValue = circshift(sampledOWFgeneration, 1);
forecastErrorSample = predictValue(2:end,:) - sampledOWFgeneration(2:end,:);
% 需要手动让预测更准确一点，不然DRCC优化模型会infeasible
error_percent = forecastErrorSample ./ repmat(max(sampledOWFgeneration),[nk-1,1]);

forecastErrorSample = forecastErrorSample/5; %最大误差控制在15%以内


meanError = mean(forecastErrorSample);
Sigma = cov(forecastErrorSample);

% sigma_sqrt = real(sqrtm(covarError));

%% test set
predictValue_all = circshift(allOWFgeneration, 1);
forecastError_all = predictValue_all(2:end,:) - allOWFgeneration(2:end,:);
forecastError_all = forecastError_all/5;
% augmentation to 8760 points

mu_all = mean(forecastError_all); Sigma_all = cov(forecastError_all);
n_test = 8760;

forecastError_test = my_mvnrnd(mu_all, Sigma_all, n_test);   % N_test x n_xi
end