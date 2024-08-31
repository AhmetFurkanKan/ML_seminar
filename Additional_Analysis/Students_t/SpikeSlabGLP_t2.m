function [store_B,store_z,store_phi,store_q,store_R2,store_gam,store_s2,y,x,u,T,k] = SpikeSlabGLP_t2(y,x,u,abeta,bbeta,Abeta,Bbeta,M,N,nu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear regression with a spike-and-slab prior and Student's t likelihood
%
% Inputs:
%   y: T x 1 vector of observations of the response variable
%   x: T x k matrix of observations of the predictors 
%   u: T x l matrix of observations of the predictors always included
%   abeta, bbeta: parameters of beta prior on q
%   Abeta, Bbeta: parameters of beta prior on R2
%   M: total number of draws
%   N: number of initial draws to be discarded
%   nu: degrees of freedom for Student's t-distribution
%
% Outputs:
%   store_B, store_z, store_phi, store_q, store_R2, store_gam, store_s2, y, x, u, T, k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(u); l = 0; else; l = size(u, 2); end
if isempty(abeta) || isempty(bbeta); abeta = 1; bbeta = 1; end
if isempty(Abeta) || isempty(Bbeta); Abeta = 1; Bbeta = 1; end
if isempty(M); M = 110000; N = 10000; end
if isempty(N) && ~isempty(M); N = round(M/11); end
if isempty(nu); nu = 3; end  % Default value for degrees of freedom

[T,k] = size(x);
varx = var(x(:,1), 1);

% Edges of grids for q and R2
edgesq = [0:.001:.1 .11:.01:.9 .901:.001:1];
edgesR2 = [0:.001:.1 .11:.01:.9 .901:.001:1];

% Storage and preliminary calculations
if l > 0; store_phi = zeros(l, M); end
store_B = zeros(k, M);
store_s2 = zeros(M, 1);
store_gam = zeros(M, 1);
store_z = zeros(k, M);
store_q = zeros(M, 1);
store_R2 = zeros(M, 1);

areaq = edgesq(2:end) - edgesq(1:end-1);
areaR2 = edgesR2(2:end) - edgesR2(1:end-1);
areaqR2 = repmat(areaq, length(areaR2), 1) .* repmat(areaR2', 1, length(areaq));
intq = edgesq(1:end-1) + areaq/2;
intR2 = edgesR2(1:end-1) + areaR2/2;
INTQ = repmat(intq, length(intR2), 1);
INTR2 = repmat(intR2', 1, length(intq));

xx = x' * x;
xy = x' * y;
yx = y' * x;
yy = y' * y;
if l > 0
    xu = x' * u;
    yu = y' * u;
    uu = u' * u;
    invuu = uu \ eye(l);
    cholinvuu = chol(invuu);
    ux = u' * x;
    invuuuy = invuu * u' * y;
    invuuux = invuu * u' * x;
elseif l == 0
    xu = zeros(k, 1);
    yu = zeros(1, 1);
    uu = zeros(1, 1);
    ux = zeros(1, k);
    phi = 0;
end

QR2 = INTQ .* (1 - INTR2) ./ INTR2;
for ttt = 0:k
    prodINT(:, :, ttt+1) = INTQ .^ (ttt + ttt/2 + abeta - 1) .* (1 - INTQ) .^ (k - ttt + bbeta - 1) .* ...
        INTR2 .^ (Abeta - 1 - ttt/2) .* (1 - INTR2) .^ (ttt/2 + Bbeta - 1) .* areaqR2;
end

% Starting values (based on Lasso under the assumption that the x's explain 50% of var(resy))
if l > 0
    phi0 = invuuuy;
    resy = y - u * phi0;
    store_phi(:, 1) = phi0;
else
    resy = y;
end
b0 = lasso(x, resy, 'Lambda', sqrt(8 * k * varx * var(resy, 1) / 2) / (2 * length(resy)), 'Standardize', 0);
store_B(:, 1) = b0;
b = b0;
s2 = sum((resy - x * b).^2) / (T - 2);  % Initial residual variance
store_s2(1) = s2;
z = (b ~= 0);
store_z(:, 1) = z;
tau = sum(z);

% Gibbs sampling
tic
for i = 2:M
    if i == 1000 * floor(.001 * i)
        i
    end
    
    % Draw q and R2
    pr_qR2 = exp(-(.5 * k * varx * b' * b / s2) .* QR2) .* squeeze(prodINT(:, :, tau+1));
    cdf_qR2 = cumsum(pr_qR2(:) / sum(pr_qR2(:)));
    aux = sum(cdf_qR2 < rand(1)) + 1;
    
    q = INTQ(aux);
    store_q(i) = q;
    R2 = INTR2(aux);
    store_R2(i) = R2;
    gam = sqrt((1 / (k * varx * q)) * R2 / (1 - R2));
    store_gam(i) = gam;
    
    % Draw phi
    if l > 0
        phihat = invuuuy - invuuux * store_B(:, i-1);
        Vphi = invuu * s2;
        phi = (randn(1, l) * sqrt(s2) * cholinvuu)' + phihat;
        store_phi(:, i) = phi';
    end
    
    % Draw z
    for j = 1:k
        z0 = z; z0(j) = 0;
        z1 = z; z1(j) = 1;
        tau0 = sum(z0);
        tau1 = sum(z1);
        W0 = xx(logical(z0), logical(z0)) + eye(tau0) / gam^2;
        W1 = xx(logical(z1), logical(z1)) + eye(tau1) / gam^2;
        bhat0 = W0 \ (xy(logical(z0)) - xu(logical(z0), :) * phi);
        bhat1 = W1 \ (xy(logical(z1)) - xu(logical(z1), :) * phi);
        
        % Log-posterior under Student's t likelihood
        log_pr_z0 = tau0 * log(q) + (k - tau0) * log(1 - q) - tau0 * log(gam) - ...
            .5 * 2 * sum(log(diag(chol(W0)))) - (T / 2) * log(1 + (resy - x(:, logical(z0)) * bhat0)' * (resy - x(:, logical(z0)) * bhat0) / (s2 * nu));

        log_pr_z1 = tau1 * log(q) + (k - tau1) * log(1 - q) - tau1 * log(gam) - ...
            .5 * 2 * sum(log(diag(chol(W1)))) - (T / 2) * log(1 + (resy - x(:, logical(z1)) * bhat1)' * (resy - x(:, logical(z1)) * bhat1) / (s2 * nu));
        
        z(j) = (rand(1) <= (1 / (exp(log_pr_z0 - log_pr_z1) + 1)));
    end
    tau = sum(z);
    W = xx(logical(z), logical(z)) + eye(tau) / gam^2;
    bhat = W \ (xy(logical(z)) - xu(logical(z), :) * phi);
    
    % Draw b
    b = zeros(k, 1);
    if tau > 0
        b(logical(z)) = (randn(1, tau) * sqrt(s2) * chol(inv(W)))' + bhat;
    end
    store_B(:, i) = b;
    
    % Draw s2 (variance)
    s2 = 1 / gamrnd(.5 * (T + nu), 2 / ((nu + (resy - x * b)' * (resy - x * b))));
    store_s2(i) = s2;
    store_z(:, i) = z;
end
timeGibbs = toc;
end
