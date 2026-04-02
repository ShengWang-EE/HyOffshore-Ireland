function stats = evaluate_constraints_violation(solution,mpc,xi_test,eps)
%% unpack data and grab dimenssions
if ~isfield(solution,'alpha_Prs')
    [solution.alpha_Prs, solution.alpha_PGs, solution.alpha_Gf, solution.alpha_Qptg, solution.alpha_Va, ...
        solution.alpha_Pg,solution.alpha_Pf] = deal(0);
end

[Prs_base, PGs_base, Gf_base, Qptg_base, Qgpp_base, Va_base, Pg_base, gamma_base, ...
    operatingCost, electricityGenerationCost, gasPurchasingCost, ...
    alpha_Prs, alpha_PGs, alpha_Gf, alpha_Qptg, alpha_Va, alpha_Pg,alpha_Pf ...
    ] = deal(...
    solution.Prs, solution.PGs, solution.Gf, solution.Qptg, solution.Qgpp, ...
    solution.Va, solution.Pg, solution.gamma, ...
    solution.objfcn, solution.electricityGenerationCost, solution.gasPurchasingCost, ...
    solution.alpha_Prs, solution.alpha_PGs, solution.alpha_Gf, solution.alpha_Qptg, solution.alpha_Va, solution.alpha_Pg, ...
    solution.alpha_Pf);

[N_test, n_xi] = size(xi_test);

ng   = size(mpc.gen,1);
nb   = size(mpc.bus,1);
i_owf = find(mpc.gentype == 1);
%% prepare some data
[~, ~, ~, ~, ~, ~, ~, ~, PMAX, PMIN] = idx_gen;
[~, ~, ~, ~, ~, RATE_A]              = idx_brch;
il = find(mpc.branch(:, RATE_A) ~= 0 & mpc.branch(:, RATE_A) < 1e10);
[B, Bf, Pbusinj, Pfinj] = makeBdc(mpc.baseMVA, mpc.bus, mpc.branch);
Pgmin = mpc.gen(:, PMIN) * 0;   % 你原来是 0
Pgmax0 = mpc.gen(:, PMAX);
Prsmin = mpc.Gbus(:,5); Prsmax = mpc.Gbus(:,6);
upf   = mpc.branch(il, RATE_A) - Pfinj(il);   % 线路上界（正负对称）
%% 
joint_violation_flag = false(N_test,1);   % 场景联合违约标记
Pg_violation_cnt     = zeros(ng,1);       % 每台机组违约次数
line_violation_cnt   = zeros(length(il),1); % 每条受限线路违约次数

Prs_violation_cnt     = zeros(nb,1);       % 每个节点气压违约次数

%% 主循环：对每个测试场景做一次"后仿真"
for k = 1:N_test
    xi_k   = xi_test(k,:).';        % n_xi x 1
    sum_xi = sum(xi_k);             % 对应你代码里的 unified affine mapping
    Pgmax = Pgmax0; Pgmax(i_owf) = Pgmax0(i_owf) + xi_k*1; % 根据不确定量修改风机上限

    % ---- 根据 affine 规则生成该场景下的实际状态 ----
    Pg_k = Pg_base + alpha_Pg * sum_xi;    % ng x 1
    Va_k = Va_base + alpha_Va * sum_xi;    % nb x 1
    
    Prs_k = Prs_base + alpha_Prs * sum_xi;
    Gf_k  = Gf_base  + alpha_Gf  * sum_xi;

    % ---- 检查机组出力上下限违反情况 ----
    % viol_Pg_lower = Pg_k < Pgmin - 1e-4;
    viol_Pg_lower = 0;
    viol_Pg_upper = Pg_k > Pgmax + 1e-4;
    viol_Pg       = viol_Pg_lower | viol_Pg_upper;

    if any(viol_Pg)
        Pg_violation_cnt = Pg_violation_cnt + viol_Pg;
    end

    % ---- 检查线路潮流越限（基于 DC）----
    Pf_k = Bf * Va_k + alpha_Pf * sum_xi;            % 所有支路潮流
    Pf_lim_viol = (Pf_k(il) > upf + 1e-4) | (Pf_k(il) < -upf - 1e-4);

    if any(Pf_lim_viol)
        line_violation_cnt = line_violation_cnt + Pf_lim_viol;
    end

    % 你还可以在这里加：
    %   气压越界、气流越界、节点功率平衡大偏差等
    % viol_Prs_lower = Prs_k < Prsmin - 1e-4;
    % viol_Prs_upper = Prs_k > Prsmax + 1e-4;
    % viol_Prs       = viol_Prs_lower | viol_Prs_upper;
    % 
    % if any(viol_Prs)
    %     Prs_violation_cnt = Prs_violation_cnt + viol_Prs;
    % end
    % ---- 联合违约（EVJP 对应的指示函数）----
    if any(viol_Pg) 
            % || any(Pf_lim_viol)
        joint_violation_flag(k) = true;
    end
end

% ==== 统计结果 ====
EVJP = sum(joint_violation_flag) / N_test;     % 联合违约概率

Pg_violation_rate   = Pg_violation_cnt   / N_test;
line_violation_rate = line_violation_cnt / N_test;

% ==== 输出结构体 ====
stats.Ntest        = N_test;
stats.eps_model    = eps;                 % 你设的机会约束水平（仅记录）
stats.EVJP         = EVJP;                % Empirical violation of joint probability
stats.count_joint  = sum(joint_violation_flag);

stats.Pg_violation_rate   = Pg_violation_rate;
stats.line_violation_rate = line_violation_rate;

end

