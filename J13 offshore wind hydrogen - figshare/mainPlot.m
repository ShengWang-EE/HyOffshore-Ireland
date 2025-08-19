%% accomadation
data2030 = readtable('Ireland future gas demand.xlsx','Sheet','results','Range','A29:H53');
data2040 = readtable('Ireland future gas demand.xlsx','Sheet','results','Range','K29:R53');
data2050 = readtable('Ireland future gas demand.xlsx','Sheet','results','Range','U29:AB53');
data2030.Generation = data2030.Generation/1e3;data2030.RemainningCapacity = data2030.RemainningCapacity/1e3;
data2040.Generation = data2040.Generation/1e3;data2040.RemainningCapacity = data2040.RemainningCapacity/1e3;
data2050.Generation = data2050.Generation/1e3;data2050.RemainningCapacity = data2050.RemainningCapacity/1e3;
data2030.Generation_1 = data2030.Generation_1/1e3;data2030.RemainningCapacity_1 = data2030.RemainningCapacity_1/1e3;
data2040.Generation_1 = data2040.Generation_1/1e3;data2040.RemainningCapacity_1 = data2040.RemainningCapacity_1/1e3;
data2050.Generation_1 = data2050.Generation_1/1e3;data2050.RemainningCapacity_1 = data2050.RemainningCapacity_1/1e3;

fig = figure;
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultAxesFontSize', 12);
colorDarker = colororder('gem12');
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
    'Color',colorDarker(1,:),'LineWidth', 1, 'MarkerSize', 4, 'MarkerEdgeColor', colorDarker(1,:), 'MarkerFaceColor', 'w');
fig4 = plot(data2030.Time_1-0.2,data2030.AccomadationRate_1,'-s',...
    'Color',colorDarker(2,:),'LineWidth', 1, 'MarkerSize', 4, 'MarkerEdgeColor', colorDarker(2,:), 'MarkerFaceColor', 'w');
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
fig2(2).DisplayName = 'Remaining capacity of OWFs';



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
    'Color',colorDarker(1,:),'LineWidth', 1, 'MarkerSize', 4, 'MarkerEdgeColor', colorDarker(1,:), 'MarkerFaceColor', 'w');
fig4 = plot(data2040.Time_1-0.2,data2040.AccomadationRate_1,'-s',...
    'Color',colorDarker(2,:),'LineWidth', 1, 'MarkerSize', 4, 'MarkerEdgeColor', colorDarker(2,:), 'MarkerFaceColor', 'w');
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
    'Color',colorDarker(1,:),'LineWidth', 1, 'MarkerSize', 4, 'MarkerEdgeColor', colorDarker(1,:), 'MarkerFaceColor', 'w');
fig4 = plot(data2050.Time_1-0.2,data2050.AccomadationRate_1,'-s',...
    'Color',colorDarker(2,:),'LineWidth', 1, 'MarkerSize', 4, 'MarkerEdgeColor', colorDarker(2,:), 'MarkerFaceColor', 'w');
ax2 = gca; % 获取当前轴的句柄
ax2.YColor = 'k'; % 设置轴线和刻度标签颜色为黑色
ylabel('Accomadation rate'); 
ylim([0,1]);
xlabel('Time (hour)'); 
text(0.00, -0.18, '(c)', 'Units', 'normalized', 'VerticalAlignment', 'top','FontWeight','bold');

set(fig, 'Position', [100, 100, 550, 550]);  % 同样的参数
exportgraphics(gcf, 'fig owf accomadation.pdf', 'ContentType', 'vector');
%% hydrogen share
