function [xhat, out] = ICDEL1(y, M, Omega, k, params, xinit)
%% Analysis Sparse Recovery via Iterative Cosupport Detection Estimation (ICDE) Algorithm
% [xhat, out] = ICDEL1(y, M, Omega, beta, params, xinit)
% 
% This aims to find an approximate (sometimes exact) solution of
%    xhat = argmin || Omega * x ||_0   subject to   || y - M * x ||_2 <= epsilon.
% ===============================
% Outputs:
%   xhat :              estimate of the target cosparse vector x0.
%   out :
%       out.Lambdahat:  estimate of the cosupport of x0.
%       out.iter:       the number of iterations
% ===============================
% Inputs:
%   y :             observation/measurement vector of a target cosparse solution x0, given by relation  y = M * x0 + noise.
%   M :             measurement matrix. This should be given as a matrix which implements linear transformation.
%   Omega :         analysis operator. Like M, this should be given as a matrix which implements linear transformation.
%   k:              cosparsity level
%   params :        parameters that govern the behavior of the algorithm (mostly).
%       params.num_iteration:                       the number of iterations.
%       params.stopping_coefficient_size:           when the maximum analysis coefficient is smaller than this, ICDE terminates.
%       params.stopping_relative_solution_change:   when the maximum relative solution_change is smaller than this, ICDE terminates.
%       params.stopping_cosparsity:                 when the minimum cosparsity is smaller than this, ICDE terminates.
%       params.noise_level :                        this corresponds to epsilon above.
%	xinit :                                      	initial estimate of x0. default : xinit = zeros(d, 1);
% ===============================
%
% Examples:
%
% d = 220;
% p = 240;
% m = 90; 
% k = 200;
% operator_type = 'random';	% random tight fram
% [x0, y, Omega, M, Lambda] = Generate_Problem(d, p, m, k, operator_type);
% beta = 0.9;
% [xhat, Lambdahat] = ICDEL1(y, M, Omega, beta);
% SupportDetection(Omega*x0, Omega*xhat)
%
% % Written by Heping Song, Jiangsu University, China, songhp@ujs.edu.cn
% ===============================


%% First, set up CVX package
% cd cvx
% cvx_setup
% cd ..
% %select SeDuMi solver for optimization
% cvx_solver sedumi

%%  
p = size(Omega, 1);
d = size(Omega, 2);

if nargin < 4
    beta = 0.7;
end

if nargin < 5
	params.num_iteration = 50;
	params.stopping_relative_solution_change = 1e-6;
	params.stopping_coefficient_size = 1e-6;
	params.stopping_residual_size = 1e-6;
    params.stopping_cosparsity = p*0.5;
    params.noise_level = 1e-6;
    xinit = zeros(d, 1);
end


iter = 0;
Lambdahat = 1:p; % cosupport


while iter < params.num_iteration
   
	iter = iter + 1; 
    
	xhat = myAnalysis_L1(y, Omega, M, Lambdahat, params.noise_level); 
  
    Omega_x = Omega * xhat;
    abscoef = abs(Omega_x);
    maxcoef = max(abscoef);
    qq = quantile(abscoef, k/p);
    to_be_removed = find(abscoef >= qq);
    Lambdahat = setdiff(1:p,to_be_removed);
            
            
            
    
    if check_stopping_criteria(xhat, xinit, maxcoef, Lambdahat, params)
        break;
    end
    
    xinit = xhat;


end

out.Lambdahat = Lambdahat;
out.iter = iter-1;
return;

function r = check_stopping_criteria(xhat, xinit, maxcoef, Lambdahat, params)

	r = 0;

    if isfield(params, 'stopping_relative_solution_change') && norm(xhat-xinit)/norm(xhat) < params.stopping_relative_solution_change
        r = 1;
        return;
    end

    if isfield(params, 'stopping_coefficient_size') && maxcoef < params.stopping_coefficient_size
        r = 1;
        return;
    end
    
    if isfield(params, 'stopping_cosparsity') && length(Lambdahat) < params.stopping_cosparsity
        r = 1;
        return;
    end
