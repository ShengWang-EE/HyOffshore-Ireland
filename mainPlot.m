%% accomadation
if ~exist(projectPath('figs'), 'dir')
    mkdir(projectPath('figs'));
end

data2030 = readtable(projectPath('data','inputs','system','Ireland future gas demand.xlsx'),'Sheet','results','Range','A29:H53');
data2040 = readtable(projectPath('data','inputs','system','Ireland future gas demand.xlsx'),'Sheet','results','Range','K29:R53');
data2050 = readtable(projectPath('data','inputs','system','Ireland future gas demand.xlsx'),'Sheet','results','Range','U29:AB53');
data2030.Generation = data2030.Generation/1e3;data2030.RemainningCapacity = data2030.RemainningCapacity/1e3;
data2040.Generation = data2040.Generation/1e3;data2040.RemainningCapacity = data2040.RemainningCapacity/1e3;
data2050.Generation = data2050.Generation/1e3;data2050.RemainningCapacity = data2050.RemainningCapacity/1e3;
data2030.Generation_1 = data2030.Generation_1/1e3;data2030.RemainningCapacity_1 = data2030.RemainningCapacity_1/1e3;
data2040.Generation_1 = data2040.Generation_1/1e3;data2040.RemainningCapacity_1 = data2040.RemainningCapacity_1/1e3;
data2050.Generation_1 = data2050.Generation_1/1e3;data2050.RemainningCapacity_1 = data2050.RemainningCapacity_1/1e3;

fig = figure;
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultAxesFontSize', 12);
colorDarker = colororder('gem12');colorDarker = colorDarker * 0.6;
colorLighter = colororder('glow12');
% fig 1
subplot1 = axes('Position', [0.1, 0.6, 0.8, 0.18]); % [left, bottom, width, height]
yyaxis left;          % 激活左侧纵坐标轴
fig1 = bar(data2030.Time-0.2,[data2030.Generation,data2030.RemainningCapacity],'stacked',...
    'BarWidth',0.3);
fig1(1).FaceColor = colorLighter(1,:); fig1(2).FaceColor = 'white';
fig1(1).FaceAlpha = 0.5; 
hold on;
fig2 = bar(data2030.Time_1+0.2,[data2030.Generation_1,data2030.RemainningCapacity_1],'stacked','BarWidth',0.3);
fig2(1).FaceColor = colorLighter(2,:); fig2(2).FaceColor = 'white';
fig2(1).FaceAlpha = 0.5; 
ax1 = gca; % 获取当前轴的句柄
ax1.YColor = 'k'; % 设置轴线和刻度标签颜色为黑色
ylabel('OWF generation (GW)');          % 设置左轴标签
ylim([0,8]);


yyaxis right;         
fig3 = plot(data2030.Time-0.2,data2030.AccomadationRate,'-o',...
    'Color',colorDarker(1,:),'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerEdgeColor', colorDarker(1,:), 'MarkerFaceColor', 'w');
fig4 = plot(data2030.Time_1-0.2,data2030.AccomadationRate_1,'-s',...
    'Color',colorDarker(2,:),'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerEdgeColor', colorDarker(2,:), 'MarkerFaceColor', 'w');
ax2 = gca; % 获取当前轴的句柄
ax2.YColor = 'k'; % 设置轴线和刻度标签颜色为黑色
ylabel('Accomadation rate'); 
ylim([0,1]);
xlabel('Time (hour)'); 

fig4.DisplayName = 'Accomadation rate in Schme 2';
fig3.DisplayName = 'Accomadation rate in Schme 1';
fig1(1).DisplayName = 'OWF generation in Scheme 1';
fig1(2).DisplayName = 'Remaining capacity of OWFs';
fig2(1).DisplayName = 'OWF generation in Scheme 2';
% fig2(2).DisplayName = 'Remaining capacity of OWFs';
legend([fig1,fig2(1),fig3,fig4]);


lgd = legend('Location', 'northoutside','NumColumns',2);
lgd.Position = [0.417 0.8 0.16 0.06];
text(0, -0.18, '(a)', 'Units', 'normalized', 'VerticalAlignment', 'top','FontWeight','bold');

% fig 2
subplot2 = axes('Position', [0.1, 0.35, 0.8, 0.18]); % [left, bottom, width, height]
yyaxis left;          % 激活左侧纵坐标轴
fig1 = bar(data2040.Time-0.2,[data2040.Generation,data2040.RemainningCapacity],'stacked',...
    'BarWidth',0.3);
fig1(1).FaceColor = colorLighter(1,:); fig1(2).FaceColor = 'white';
fig1(1).FaceAlpha = 0.5; 
hold on;
fig2 = bar(data2040.Time_1+0.2,[data2040.Generation_1,data2040.RemainningCapacity_1],'stacked','BarWidth',0.3);
fig2(1).FaceColor = colorLighter(2,:); fig2(2).FaceColor = 'white';
fig2(1).FaceAlpha = 0.5; 
ax1 = gca; % 获取当前轴的句柄
ax1.YColor = 'k'; % 设置轴线和刻度标签颜色为黑色
ylabel('OWF generation (GW)');          % 设置左轴标签
ylim([0,23]);
fig1(1).DisplayName = 'OWF generation in Scheme 1';
fig2(1).DisplayName = 'OWF generation in Scheme 2';
fig2(2).DisplayName = 'Remaining capacity of OWFs';
% legend('Location', 'northoutside','NumColumns',3);

yyaxis right;         
fig3 = plot(data2040.Time-0.2,data2040.AccomadationRate,'-o',...
    'Color',colorDarker(1,:),'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerEdgeColor', colorDarker(1,:), 'MarkerFaceColor', 'w');
fig4 = plot(data2040.Time_1-0.2,data2040.AccomadationRate_1,'-s',...
    'Color',colorDarker(2,:),'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerEdgeColor', colorDarker(2,:), 'MarkerFaceColor', 'w');
ax2 = gca; % 获取当前轴的句柄
ax2.YColor = 'k'; % 设置轴线和刻度标签颜色为黑色
ylabel('Accomadation rate'); 
ylim([0,1]);
xlabel('Time (hour)'); 
text(0, -0.18, '(b)', 'Units', 'normalized', 'VerticalAlignment', 'top','FontWeight','bold');

% fig 3
subplot1 = axes('Position', [0.1, 0.1, 0.8, 0.18]); % [left, bottom, width, height]
yyaxis left;          % 激活左侧纵坐标轴
fig1 = bar(data2050.Time-0.2,[data2050.Generation,data2050.RemainningCapacity],'stacked',...
    'BarWidth',0.3);
fig1(1).FaceColor = colorLighter(1,:); fig1(2).FaceColor = 'white';
fig1(1).FaceAlpha = 0.5; 
hold on;
fig2 = bar(data2050.Time_1+0.2,[data2050.Generation_1,data2050.RemainningCapacity_1],'stacked','BarWidth',0.3);
fig2(1).FaceColor = colorLighter(2,:); fig2(2).FaceColor = 'white';
fig2(1).FaceAlpha = 0.5; 
ax1 = gca; % 获取当前轴的句柄
ax1.YColor = 'k'; % 设置轴线和刻度标签颜色为黑色
ylabel('OWF generation (GW)');          % 设置左轴标签
ylim([0,40]);
fig1(1).DisplayName = 'OWF generation in Scheme 1';
fig2(1).DisplayName = 'OWF generation in Scheme 2';
fig2(2).DisplayName = 'Remaining capacity of OWFs';
% legend('Location', 'northoutside','NumColumns',3);

yyaxis right;         
fig3 = plot(data2050.Time-0.2,data2050.AccomadationRate,'-o',...
    'Color',colorDarker(1,:),'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerEdgeColor', colorDarker(1,:), 'MarkerFaceColor', 'w');
fig4 = plot(data2050.Time_1-0.2,data2050.AccomadationRate_1,'-s',...
    'Color',colorDarker(2,:),'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerEdgeColor', colorDarker(2,:), 'MarkerFaceColor', 'w');
ax2 = gca; % 获取当前轴的句柄
ax2.YColor = 'k'; % 设置轴线和刻度标签颜色为黑色
ylabel('Accomadation rate'); 
ylim([0,1]);
xlabel('Time (hour)'); 
text(0.00, -0.18, '(c)', 'Units', 'normalized', 'VerticalAlignment', 'top','FontWeight','bold');

set(fig, 'Position', [100, 100, 550, 800]);  % 同样的参数
exportgraphics(gcf, projectPath('figs','fig owf accomadation.pdf'), 'ContentType', 'vector');
%% power demand
clc
clear
peakDemand = [7.02
7.38
7.59
7.8
8.02
8.22
8.35
8.48
8.64
8.82
9.03
];
year = [2022:2032]';
p = polyfit(year, peakDemand, 1);
futureYear = [2033:2050]';
futureDemand = polyval(p, futureYear);
demandAll = [peakDemand;futureDemand];
yearAll = [year;futureYear];
fig = figure;
plot(year,peakDemand,'-');
hold on;
plot(futureYear,futureDemand,'--');
set(fig, 'Position', [100, 100, 600, 200]);  % 同样的参数
xlabel('Year'); xlim([2022,2050]);
ylabel('Peak electricity demand (GW)')
exportgraphics(gcf, projectPath('figs','fig power demand forecast.pdf'), 'ContentType', 'vector');
%%
clear
clc
EdemandCurve = [3987
3920
3849
3825
3773
3762
3716
3687
3643
3623
3595
3568
3541
3547
3540
3528
3513
3523
3516
3516
3533
3589
3637
3680
3818
4012
4126
4315
4516
4839
5006
5163
5240
5276
5269
5301
5371
5410
5418
5435
5423
5404
5425
5423
5450
5445
5445
5430
5440
5478
5494
5491
5500
5458
5410
5376
5367
5400
5360
5364
5384
5444
5484
5542
5585
5699
5851
5984
6043
6102
6100
6091
6013
5921
5863
5782
5740
5720
5645
5578
5542
5475
5377
5279
5211
5145
5013
4896
4793
4657
4546
4426
4351
4333
4241
4150
];
GdemandCurve = [1
0.98564325
0.974523157
0.969008808
0.969923334
0.979108775
0.989382623
0.995857754
0.994495774
0.977663023
0.954995188
0.932144461
0.914881069
0.912315649
0.918711513
0.931030891
0.945154315
0.952691139
0.957865502
0.961541903
0.965205426
0.972472834
0.981014411
0.99020104
0.999213682
1.005318957
1.011362842
1.018264111
1.026734814
1.037633777
1.049555657
1.061231538
1.071161785
1.077024498
1.079323629
1.078113802
1.073913898
1.067269189
1.059801717
1.052771696
1.047179084
1.043152473
1.041400916
1.041598677
1.042854335
1.043755601
1.043026745
1.039403145
1.032042461
1.020038375
1.005141643
0.98900788
0.973658345
0.96039161
0.952091307
0.95013172
0.954786711
0.966506828
0.981114192
0.995355292
1.006290545
1.008266139
1.006430189
1.002913548
0.999606841
1.001131423
1.000941097
0.995525964
0.981718105
0.953466588
0.917877834
0.87927912
0.842539807
0.813957298
0.793261837
0.780559568
0.774523295
0.771786463
0.771488135
0.771786888
0.771378022
0.768697364
0.765887369
0.764680748
0.767053273
0.775242362
0.789707717
0.810664584
0.837625319
0.868528289
0.903554386
0.941553982
0.981128366
1.02030541
1.058147135
1.093331103
];
EdemandFactor = EdemandCurve / max(EdemandCurve);
GdemandFactor = GdemandCurve / max(GdemandCurve);
fig = figure;
hours = [0:0.25:23.75]';
plot(hours,[EdemandFactor,GdemandFactor]);
xlabel('Time (h)');
ylabel('Demand curve');
lgd = legend({'Electricity demand curve','Gas demand curve'});
lgd.Location = 'southeast';
set(fig, 'Position', [100, 100, 600, 200]);  % 同样的参数

exportgraphics(gcf, projectPath('figs','fig demand curve.pdf'), 'ContentType', 'vector');
%% demand
clear
clc
fig = figure;
subplot(2,2,[1 2]);
gasDemand = [12.76551	5.04116	2.58382	0.00903	0.14455	0.58723	7.72435
12.72034	5.58322	2.71933	0.00903	0.16262	0.58723	7.22746
14.68982	8.86268	4.2642	0.01807	0.23489	0.59627	9.54929
14.34652	8.85364	4.12869	0.03614	0.21682	0.55109	9.52218
13.51536	9.35053	4.06545	0.05421	0.22586	0.48785	9.32343
13.56053	8.32965	4.01124	0.08131	0.23489	0.48785	9.30536
13.37081	8.42	3.948	0.10841	0.23489	0.47882	9.81128
13.03654	8.31158	3.99317	0.11745	0.23489	0.45172	9.16985
12.44027	8.33869	3.93897	0.12648	0.23489	0.47882	8.63682
13.01847	8.05862	2.90002	0.12648	0.19876	0.47882	8.77233
13.10377	7.81639	3.34148	0.12142	0.1874	0.45738	8.41421
12.69281	7.41658	3.30707	0.11079	0.17651	0.46421	8.16847
12.63355	7.15003	3.23272	0.09679	0.17567	0.47683	8.00908
12.72369	6.83748	3.20474	0.08059	0.17109	0.47107	7.76903
12.608	6.644	3.13826	0.06371	0.16651	0.46997	7.56519
12.47265	6.40095	3.10101	0.04724	0.16399	0.48037	7.37335
12.4417	6.22178	3.04679	0.03208	0.15914	0.48466	7.17304
12.39477	5.98909	3.00241	0.01881	0.15641	0.48474	6.98473
12.30236	5.80046	2.95556	0.00778	0.15239	0.49112	6.80122
12.22965	5.57756	2.90917	0.00778	0.14917	0.49857	6.62059];

year = [2021:2040]';
fig1 = area(year,gasDemand,'LineStyle', 'none');
legend(fig1, {'Power', 'Industrial & Commercial', 'Residential','Transport','Own use','Isle of Man','Northern Ireland'}, ...
    'Location', 'NorthEast','NumColumns', 3);
xlabel('Year'); xlim([2021,2040]);
ylim([10,60]);
ylabel('Peak gas demand (Mm^3/day)')
text(0, -0.15, '(a)', 'Units', 'normalized', 'VerticalAlignment', 'top','FontWeight','bold');

subplot(2,2,3);
peakDemand = [7.02
7.38
7.59
7.8
8.02
8.22
8.35
8.48
8.64
8.82
9.03
];
year = [2022:2032]';
p = polyfit(year, peakDemand, 1);
futureYear = [2033:2040]';
futureDemand = polyval(p, futureYear);
demandAll = [peakDemand;futureDemand];
yearAll = [year;futureYear];
plot(year,peakDemand,'-');
hold on;
plot(futureYear,futureDemand,'--');
xlabel('Year'); xlim([2022,2040]);
ylabel('Peak electricity demand (GW)')
text(0, -0.15, '(b)', 'Units', 'normalized', 'VerticalAlignment', 'top','FontWeight','bold');

subplot(2,2,4);
EdemandCurve = [3987
3920
3849
3825
3773
3762
3716
3687
3643
3623
3595
3568
3541
3547
3540
3528
3513
3523
3516
3516
3533
3589
3637
3680
3818
4012
4126
4315
4516
4839
5006
5163
5240
5276
5269
5301
5371
5410
5418
5435
5423
5404
5425
5423
5450
5445
5445
5430
5440
5478
5494
5491
5500
5458
5410
5376
5367
5400
5360
5364
5384
5444
5484
5542
5585
5699
5851
5984
6043
6102
6100
6091
6013
5921
5863
5782
5740
5720
5645
5578
5542
5475
5377
5279
5211
5145
5013
4896
4793
4657
4546
4426
4351
4333
4241
4150
];
GdemandCurve = [1
0.98564325
0.974523157
0.969008808
0.969923334
0.979108775
0.989382623
0.995857754
0.994495774
0.977663023
0.954995188
0.932144461
0.914881069
0.912315649
0.918711513
0.931030891
0.945154315
0.952691139
0.957865502
0.961541903
0.965205426
0.972472834
0.981014411
0.99020104
0.999213682
1.005318957
1.011362842
1.018264111
1.026734814
1.037633777
1.049555657
1.061231538
1.071161785
1.077024498
1.079323629
1.078113802
1.073913898
1.067269189
1.059801717
1.052771696
1.047179084
1.043152473
1.041400916
1.041598677
1.042854335
1.043755601
1.043026745
1.039403145
1.032042461
1.020038375
1.005141643
0.98900788
0.973658345
0.96039161
0.952091307
0.95013172
0.954786711
0.966506828
0.981114192
0.995355292
1.006290545
1.008266139
1.006430189
1.002913548
0.999606841
1.001131423
1.000941097
0.995525964
0.981718105
0.953466588
0.917877834
0.87927912
0.842539807
0.813957298
0.793261837
0.780559568
0.774523295
0.771786463
0.771488135
0.771786888
0.771378022
0.768697364
0.765887369
0.764680748
0.767053273
0.775242362
0.789707717
0.810664584
0.837625319
0.868528289
0.903554386
0.941553982
0.981128366
1.02030541
1.058147135
1.093331103
];
EdemandFactor = EdemandCurve / max(EdemandCurve);
GdemandFactor = GdemandCurve / max(GdemandCurve);
hours = [0:0.25:23.75]';
plot(hours,[EdemandFactor,GdemandFactor]);
xlabel('Time (h)');
ylabel('Demand curve');
lgd = legend({'Electricity demand curve','Gas demand curve'});
lgd.Location = 'southeast';

text(0, -0.15, '(c)', 'Units', 'normalized', 'VerticalAlignment', 'top','FontWeight','bold');
% set(fig, 'Position', [100, 100, 600, 300]);  % 同样的参数
exportgraphics(gcf, projectPath('figs','fig future energy demand.pdf'), 'ContentType', 'vector');

%% accomadation
clear
clc
data2030 = readtable(projectPath('data','inputs','system','Ireland future gas demand.xlsx'),'Sheet','results','Range','A29:H53');
data2040 = readtable(projectPath('data','inputs','system','Ireland future gas demand.xlsx'),'Sheet','results','Range','K29:R53');
data2050 = readtable(projectPath('data','inputs','system','Ireland future gas demand.xlsx'),'Sheet','results','Range','U29:AB53');
data2030.Generation = data2030.Generation/1e3;data2030.RemainningCapacity = data2030.RemainningCapacity/1e3;
data2040.Generation = data2040.Generation/1e3;data2040.RemainningCapacity = data2040.RemainningCapacity/1e3;
data2050.Generation = data2050.Generation/1e3;data2050.RemainningCapacity = data2050.RemainningCapacity/1e3;
data2030.Generation_1 = data2030.Generation_1/1e3;data2030.RemainningCapacity_1 = data2030.RemainningCapacity_1/1e3;
data2040.Generation_1 = data2040.Generation_1/1e3;data2040.RemainningCapacity_1 = data2040.RemainningCapacity_1/1e3;
data2050.Generation_1 = data2050.Generation_1/1e3;data2050.RemainningCapacity_1 = data2050.RemainningCapacity_1/1e3;

fig = figure;
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultAxesFontSize', 12);
colorDarker = colororder('gem12');colorDarker = colorDarker * 0.6;
colorLighter = colororder('glow12');
% fig 1
subplot1 = axes('Position', [0.1, 0.5, 0.8, 0.25]); % [left, bottom, width, height]
yyaxis left;          % 激活左侧纵坐标轴
fig1 = bar(data2030.Time-0.2,[data2030.Generation,data2030.RemainningCapacity],'stacked',...
    'BarWidth',0.3);
fig1(1).FaceColor = colorLighter(1,:); fig1(2).FaceColor = 'white';
fig1(1).FaceAlpha = 0.5; 
hold on;
fig2 = bar(data2030.Time_1+0.2,[data2030.Generation_1,data2030.RemainningCapacity_1],'stacked','BarWidth',0.3);
fig2(1).FaceColor = colorLighter(2,:); fig2(2).FaceColor = 'white';
fig2(1).FaceAlpha = 0.5; 
ax1 = gca; % 获取当前轴的句柄
ax1.YColor = 'k'; % 设置轴线和刻度标签颜色为黑色
ylabel('OWF generation (GW)');          % 设置左轴标签
ylim([0,8]);


yyaxis right;         
fig3 = plot(data2030.Time-0.2,data2030.AccomadationRate,'-o',...
    'Color',colorDarker(1,:),'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerEdgeColor', colorDarker(1,:), 'MarkerFaceColor', 'w');
fig4 = plot(data2030.Time_1-0.2,data2030.AccomadationRate_1,'-s',...
    'Color',colorDarker(2,:),'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerEdgeColor', colorDarker(2,:), 'MarkerFaceColor', 'w');
ax2 = gca; % 获取当前轴的句柄
ax2.YColor = 'k'; % 设置轴线和刻度标签颜色为黑色
ylabel('Accomadation rate'); 
ylim([0,1]);
xlabel('Time (hour)'); 

fig4.DisplayName = 'Accomadation rate in Schme 2';
fig3.DisplayName = 'Accomadation rate in Schme 1';
fig1(1).DisplayName = 'OWF generation in Scheme 1';
fig1(2).DisplayName = 'Remaining capacity of OWFs';
fig2(1).DisplayName = 'OWF generation in Scheme 2';
% fig2(2).DisplayName = 'Remaining capacity of OWFs';
legend([fig1,fig2(1),fig3,fig4]);


lgd = legend('Location', 'northoutside','NumColumns',2);
lgd.Position = [0.417 0.8 0.16 0.06];
text(0, -0.18, '(a)', 'Units', 'normalized', 'VerticalAlignment', 'top','FontWeight','bold');

% fig 2
subplot2 = axes('Position', [0.1, 0.14, 0.8, 0.25]); % [left, bottom, width, height]
yyaxis left;          % 激活左侧纵坐标轴
fig1 = bar(data2040.Time-0.2,[data2040.Generation,data2040.RemainningCapacity],'stacked',...
    'BarWidth',0.3);
fig1(1).FaceColor = colorLighter(1,:); fig1(2).FaceColor = 'white';
fig1(1).FaceAlpha = 0.5; 
hold on;
fig2 = bar(data2040.Time_1+0.2,[data2040.Generation_1,data2040.RemainningCapacity_1],'stacked','BarWidth',0.3);
fig2(1).FaceColor = colorLighter(2,:); fig2(2).FaceColor = 'white';
fig2(1).FaceAlpha = 0.5; 
ax1 = gca; % 获取当前轴的句柄
ax1.YColor = 'k'; % 设置轴线和刻度标签颜色为黑色
ylabel('OWF generation (GW)');          % 设置左轴标签
ylim([0,23]);
fig1(1).DisplayName = 'OWF generation in Scheme 1';
fig2(1).DisplayName = 'OWF generation in Scheme 2';
fig2(2).DisplayName = 'Remaining capacity of OWFs';
% legend('Location', 'northoutside','NumColumns',3);

yyaxis right;         
fig3 = plot(data2040.Time-0.2,data2040.AccomadationRate,'-o',...
    'Color',colorDarker(1,:),'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerEdgeColor', colorDarker(1,:), 'MarkerFaceColor', 'w');
fig4 = plot(data2040.Time_1-0.2,data2040.AccomadationRate_1,'-s',...
    'Color',colorDarker(2,:),'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerEdgeColor', colorDarker(2,:), 'MarkerFaceColor', 'w');
ax2 = gca; % 获取当前轴的句柄
ax2.YColor = 'k'; % 设置轴线和刻度标签颜色为黑色
ylabel('Accomadation rate'); 
ylim([0,1]);
xlabel('Time (hour)'); 
text(0, -0.18, '(b)', 'Units', 'normalized', 'VerticalAlignment', 'top','FontWeight','bold');



set(fig, 'Position', [100, 100, 550, 400]);  % 同样的参数
exportgraphics(gcf, projectPath('figs','fig owf accomadation2.pdf'), 'ContentType', 'vector');

%% ====================== response R2 ================================
%% 
% 数据（两列）
data = [
    0.000913   0
    0.036301   0.000114
    0.10137    0.007877
    0.178196   0.04395
    0.210342   0.105365
];

% 横坐标 eps
eps_values = [0.01 0.05 0.1 0.15 0.2];

% 绘图
figure('Position',[300 300 600 280]); % ← 减少高度，让图更扁
h = bar(eps_values, data, 'grouped');
% 小清新颜色
h(1).FaceColor = [0.55 0.75 0.95];   % 柔和浅蓝
h(2).FaceColor = [0.99 0.75 0.65];   % 温柔浅橙

xlabel('\epsilon');
ylabel('Constraint violation probability');
legend('M3: Wassertein-distance based DRCC', 'M5: Moment based DRCC', 'Location', 'northwest');
grid on;

% ---- 添加数值标注（四位小数） ----
xtips1 = h(1).XEndPoints;
ytips1 = h(1).YEndPoints;
labels1 = string(num2str(h(1).YData','%.2f'));
text(xtips1, ytips1, labels1, 'HorizontalAlignment','center', 'VerticalAlignment','bottom');

xtips2 = h(2).XEndPoints;
ytips2 = h(2).YEndPoints;
labels2 = string(num2str(h(2).YData','%.2f'));
text(xtips2, ytips2, labels2, 'HorizontalAlignment','center', 'VerticalAlignment','bottom');

exportgraphics(gcf, projectPath('figs','fig constraint violation probability.pdf'), 'ContentType', 'vector');

%% cost
data = [
    239929.007	251890.6622
    232251.1392	233528
    231715.7357	232325.7419
    231467.3641	231974.8282
    231326.1899	231750.387
];

% 横坐标 eps
eps_values = [0.01 0.05 0.1 0.15 0.2];

% 绘图
figure('Position',[300 300 600 280]); % ← 减少高度，让图更扁
h = bar(eps_values, data, 'grouped');
% 小清新颜色
h(1).FaceColor = [0.55 0.75 0.95];   % 柔和浅蓝
h(2).FaceColor = [0.99 0.75 0.65];   % 温柔浅橙
ylim([231000 252000]);

xlabel('\epsilon');
ylabel('Operating cost (€)');
legend('M3: Wassertein-distance based DRCC', 'M5: Moment based DRCC', 'Location', 'northeast');
grid on;


exportgraphics(gcf, projectPath('figs','fig operating costs.pdf'), 'ContentType', 'vector');
%%
% ================== 合成大图（上下排列） ==================
figure('Position',[100 100 600 600]);  % 竖排更高一些

% ----------- 子图 (a) Constraint violation probability ----------
subplot(2,1,1);

% 数据
data1 = [
    0.000913   0
    0.036301   0.000114
    0.10137    0.007877
    0.178196   0.04395
    0.210342   0.105365
];
eps_values = [0.01 0.05 0.1 0.15 0.2];

h = bar(eps_values, data1, 'grouped');
h(1).FaceColor = [0.55 0.75 0.95];
h(2).FaceColor = [0.99 0.75 0.65];

ylabel('Constraint violation probability');
legend('M3: Wassertein DRCC', 'M5: Moment DRCC', ...
       'Location', 'northwest');
grid on;

% 数值标注
xtips1 = h(1).XEndPoints;
ytips1 = h(1).YEndPoints;
labels1 = string(num2str(h(1).YData','%.2f'));
text(xtips1, ytips1, labels1, 'HorizontalAlignment','center', 'VerticalAlignment','bottom');

xtips2 = h(2).XEndPoints;
ytips2 = h(2).YEndPoints;
labels2 = string(num2str(h(2).YData','%.2f'));
text(xtips2, ytips2, labels2, 'HorizontalAlignment','center', 'VerticalAlignment','bottom');

% 子图 (a) 标签
text(-0.02, -0.15, '(a)', 'Units','normalized', 'FontSize', 10, 'FontWeight','bold');


% ----------- 子图 (b) Operating cost ----------
subplot(2,1,2);

data2 = [
    239929.007	251890.6622
    232251.1392	233528
    231715.7357	232325.7419
    231467.3641	231974.8282
    231326.1899	231750.387
];

h = bar(eps_values, data2, 'grouped');
h(1).FaceColor = [0.55 0.75 0.95];
h(2).FaceColor = [0.99 0.75 0.65];

ylim([231000 252000]);
xlabel('\epsilon');
ylabel('Operating cost (€)');
legend('M3: Wassertein DRCC', 'M5: Moment DRCC', ...
       'Location', 'northeast');
grid on;

% 子图 (b) 标签
text(-0.02, -0.15, '(b)', 'Units','normalized', 'FontSize', 10, 'FontWeight','bold');

exportgraphics(gcf, projectPath('figs','fig_combined_DRCC.pdf'), ...
    'ContentType','vector', 'BackgroundColor','none');
%% convergence
% Data
iter = 1:6;

TotalCost   = [5.19E-04, 1.17E-09, 1.10E-10, 1.89E-12, 1.01E-13, 2.42E-15];
GasFlow     = [0.00188, 1.84E-04, 1.49E-05, 1.16E-06, 2.64E-07, 9.53E-09];
GasComp     = [1.99E-05, 3.15E-06, 1.87E-07, 7.43E-09, 1.32E-09, 4.24E-10];
GasDemand   = [0.16524, 0.1618, 0.16155, 0.16152, 0.16152, 2.84E-05];
GasConsGPP  = [2.28E-04, 4.04E-05, 3.25E-06, 2.53E-07, 5.73E-08, 4.15E-08];

% Plot
figure('Position',[200 200 600 300]); % Wide figure like your example
semilogy(iter, TotalCost, '-s', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
semilogy(iter, GasFlow, '-o', 'LineWidth', 1.5, 'MarkerSize', 8);
semilogy(iter, GasComp, '-^', 'LineWidth', 1.5, 'MarkerSize', 8);
semilogy(iter, GasDemand, '-v', 'LineWidth', 1.5, 'MarkerSize', 8);
semilogy(iter, GasConsGPP, '-d', 'LineWidth', 1.5, 'MarkerSize', 8);

% Formatting
grid on;
xlabel('Iteration');
ylabel('Relative value of slack variables');

legend({'Total cost', 'Gas flow', 'Gas composition', ...
        'Gas demand', 'Gas consumption of GPP'}, ...
        'NumColumns', 2, 'Location', 'southwest');

set(gca, 'FontSize', 12);

exportgraphics(gcf, projectPath('figs','fig convergence.pdf'), ...
    'ContentType','vector', 'BackgroundColor','none');

%%

% ================== 数据 ==================
% Accomodation rate (10×10)
A = [
0.4679 0.4679 0.4681 0.4678 0.4701 0.4679 0.4678 0.4679 0.4352 0.3675;
0.4665 0.4665 0.4665 0.4665 0.4665 0.4665 0.4665 0.4665 0.4343 0.3674;
0.4649 0.4649 0.4649 0.4649 0.4649 0.4649 0.4649 0.4649 0.4343 0.3674;
0.4630 0.4630 0.4630 0.4630 0.4630 0.4630 0.4630 0.4630 0.4343 0.3678;
0.4607 0.4614 0.4607 0.4607 0.4607 0.4607 0.4607 0.4607 0.4343 0.3677;
0.4580 0.4580 0.4580 0.4580 0.4580 0.4580 0.4580 0.4580 0.4343 0.3677;
0.4545 0.4545 0.4545 0.4545 0.4545 0.4545 0.4545 0.4545 0.4326 0.3635;
0.4501 0.4501 0.4501 0.4501 0.4501 0.4501 0.4501 0.4501 0.4301 0.3608;
0.4442 0.4442 0.4442 0.4442 0.4442 0.4442 0.4442 0.4442 0.4262 0.3727;
0.4179 0.4186 0.4181 0.4186 0.4162 0.4187 0.4179 0.4186 0.4167 0.3691];

% Carbon emission (10×10)
C = [
87.1302258 87.1290426 87.5017506 87.1432308 88.5073074 87.1290426 87.1290426 87.1290426 86.2590540 86.7122094;
85.5727572 85.5855582 85.5727572 85.5715842 85.5713190 85.5906582 85.5716250 85.5714822 85.1731110 85.1731824;
83.7583710 83.7583914 83.7583914 83.7583914 83.7583914 83.7583914 83.7583914 83.7583914 81.8530008 81.8529906;
81.6157080 81.6159222 81.6157794 81.6157692 81.6160242 81.6159426 81.6159120 81.6158202 79.0882296 79.1029074;
79.0478376 79.7974662 79.0494084 79.0480212 79.0465728 79.0479498 79.0493778 79.0466136 75.8026872 79.0625460;
75.9129798 75.9129798 75.9129798 75.9129798 75.9129798 75.9129798 75.9129900 75.9129798 71.8288386 76.0767612;
72.0005454 72.0005352 72.0005352 72.0005352 72.0005352 72.0005352 72.0005352 72.0005352 67.0152750 67.5157482;
66.9800544 66.9800544 66.9800544 66.9800544 66.9800544 66.9800544 66.9800544 66.9800544 60.9648390 63.0096738;
60.3034710 60.3034710 60.3034710 60.3034710 60.3034710 60.3034710 60.3034710 60.3033078 53.8768794 60.4860000;
30.9859578 31.2161004 31.3401120 31.2161004 30.3213972 31.3785354 31.2917742 31.2161004 31.1559102 29.9588688];

% 横纵坐标（0.1~1.0）
x = 0.1:0.1:1.0;   % Max hydrogen composition
y = 0.1:0.1:1.0;   % Proportion to original deviation

% ================== 画图 ==================
figure('Position',[200 100 600 600]);

% ---------- (a) Accomodation rate ----------
subplot(2,1,1);
imagesc(x, y, A);        % 或者 imagesc(A); 自己设 xticks/yticks
set(gca,'YDir','normal');    % y 轴向上递增
colormap(gca, summer);       % 绿色系
colorbar;

% 坐标轴设置
xlabel('Max hydrogen composition');
ylabel('Proportion to original deviation');
title('Accomodation rate');

% 坐标刻度和标签
xticks(x);
xticklabels(arrayfun(@(v)sprintf('%.1f',v), x,'UniformOutput',false));
yticks(y);
yticklabels(arrayfun(@(v)sprintf('%.1f',v), y,'UniformOutput',false));

% 画网格线（白色小格子效果）
set(gca,'XGrid','on','YGrid','on',...
        'GridColor',[1 1 1],'GridAlpha',1,...
        'Layer','top','FontSize',10);
% 子图 (a) 标签
text(-0.02, -0.15, '(a)', 'Units','normalized', 'FontSize', 10, 'FontWeight','bold');

% ---------- (b) Carbon emission ----------
subplot(2,1,2);
imagesc(x, y, C);
set(gca,'YDir','normal');
colormap(gca, summer);   % 可以改成别的，如 parula / hot / etc.
cb = colorbar;
ylabel(cb,'Carbon emission (Mm^3/day)');

xlabel('Max hydrogen composition');
ylabel('Proportion to original deviation');
title('Carbon emission (Mm^3/day)');

xticks(x);
xticklabels(arrayfun(@(v)sprintf('%.1f',v), x,'UniformOutput',false));
yticks(y);
yticklabels(arrayfun(@(v)sprintf('%.1f',v), y,'UniformOutput',false));

set(gca,'XGrid','on','YGrid','on',...
        'GridColor',[1 1 1],'GridAlpha',1,...
        'Layer','top','FontSize',10);

% 子图 (b) 标签
text(-0.02, -0.15, '(b)', 'Units','normalized', 'FontSize', 10, 'FontWeight','bold');

exportgraphics(gcf, projectPath('figs','fig sensitivity.pdf'), ...
    'ContentType','vector', 'BackgroundColor','none');
