function [xhat, Lambdahat] = myAnalysis_L1(y, Omega, M, Lambdahat, noise_level)

Omega = Omega(Lambdahat,:);
zeroTol = 1e-6;
d = size(Omega, 2);

%% Solving for x from y
cvx_begin
    cvx_quiet(true);
    cvx_precision best    
    variables xhat(d);    
    minimize(norm(Omega*xhat, 1));
    subject to
%     y == M*xhat;
    norm(y - M*xhat) <= noise_level;
    
cvx_end

Lambdahat = find(abs(Omega*xhat) < zeroTol);
xhat = pinv([M; Omega(Lambdahat, :)]) * [y; zeros(length(Lambdahat),1)];
