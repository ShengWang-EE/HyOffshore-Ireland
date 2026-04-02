clear
clc
% WI error
GCV.ng = 37.9; GCV.hy = 12.75;
RD.ng = 0.607; RD.hy = 0.069;

count = 1;
for hyComp = 0:0.01:1
    WI_true(count) = (GCV.ng * (1-hyComp) + GCV.hy * hyComp) ...
        / sqrt(RD.ng * (1-hyComp) + RD.hy * hyComp);
    WI_estimate(count) = 2 * (GCV.ng * (1-hyComp) + GCV.hy * hyComp) ...
        / (sqrt(RD.ng) + (RD.ng * (1-hyComp) + RD.hy * hyComp) / sqrt(RD.ng));
    relativeError(count) = abs(WI_true - WI_estimate) / (WI_true + WI_estimate) * 2;
    count = count + 1;
end

% 作图
figure;
yyaxis left
plot(hyComp, WI_true, 'b-', 'LineWidth', 1.5); hold on;
plot(hyComp, WI_estimate, 'g--', 'LineWidth', 1.5);
ylabel('Wobbe Index (MJ/Nm^3)');

yyaxis right
plot(hyComp, relativeError, 'r-', 'LineWidth', 1.5);
ylabel('Relative Error');

xlabel('Hydrogen Fraction');
title('Wobbe Index True vs Estimate and Relative Error');
legend('WI True','WI Estimate','Relative Error','Location','best');
grid on;