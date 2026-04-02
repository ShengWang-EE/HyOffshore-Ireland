function X = my_mvnrnd(mu, Sigma, n)
%MY_MVNRND Generate multivariate normal random vectors (no toolboxes)
%   X = MY_MVNRND(mu, Sigma, n) returns an n-by-d matrix X whose rows are
%   drawn from N(mu, Sigma).
%
%   mu    : 1 x d mean vector
%   Sigma : d x d covariance matrix
%   n     : number of samples

    % 保证 mu 是行向量
    mu = mu(:)';               
    d  = length(mu);

    % 数值上强制对称
    Sigma = (Sigma + Sigma')/2;

    % 先尝试 Cholesky
    [~, p] = chol(Sigma);
    if p == 0
        L = chol(Sigma, 'lower');      % 正定情况
    else
        % 不是严格正定，用特征分解修复
        [U, D] = eig(Sigma);
        D = max(D, 0);                 % 把负特征值截断为 0
        L = U * sqrt(D);
    end

    % Z ~ N(0, I)
    Z = randn(n, d);                   % n x d

    % X = mu + L*Z  （注意行列对应）
    X = Z * L' + repmat(mu, n, 1);     % n x d
end
