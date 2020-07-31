function [xhat, out] = ICDEL2(y, M, MH, Omega, OmegaH, xinit, k, params)
%% Analysis Sparse Recovery via Iterative Cosupport Detection Estimation (ICDE) Algorithm
% [xhat, out] = ICDEL2a(y, M, MH, Omega, OmegaH, xinit, beta, params)
% 
% This aims to find an approximate (sometimes exact) solution of
%    xhat = argmin || Omega * x ||_0   subject to   || y - M * x ||_2 <= epsilon.
% ===============================
% Outputs:
%   xhat                estimate of the target cosparse vector x0.
%   out :
%       out.Lambdahat:  estimate of the cosupport of x0.
%       out.iter:       the number of iterations
% ===============================
% Inputs:
%	y :         observation/measurement vector of a target cosparse solution x0, 
%               given by relation  y = M * x0 + noise.
%   M :         measurement matrix. This should be given either as a matrix or as a function
%               handle which implements linear transformation.
%   MH :        conjugate transpose of M. 
%   Omega :     analysis operator. Like M, this should be given either as a matrix or as a function
%               handle which implements linear transformation.
%   OmegaH :	conjugate transpose of OmegaH.
%	xinit :     initial estimate of x0 that ICDE will start with. can be zeros(d, 1).
%   k:          cosparsity level
%
%
%   params :    parameters that govern the behavior of the algorithm (mostly).
%       params.num_iteration :                      maximum number of iterations.
%       params.stopping_relative_solution_change:   when the maximum relative solution_change is smaller than this, ICDE terminates.
%       params.stopping_coefficient_size :          when the maximum analysis coefficient is smaller than this, ICDE terminates.
%       params.stopping_lagrange_multiplier_size:   when the minimum lagrange multiplier is larger than this, ICDE terminates.
%       params.stopping_cosparsity:                 when the minimum cosparsity is smaller than this, ICDE terminates.
%       params.max_inner_iteration:                 maximum number of iterations in conjugate gradient method.
%       params.noise_level :                        this corresponds to epsilon above.
%       params.l2solver :                           legitimate values are 'pseudoinverse' or 'cg'. determines which method
%                                                       is used to compute
%                                                       argmin || Omega_Lambdahat * x ||_2   subject to  || y - M * x ||_2 <= epsilon.
%       params.l2_accuracy :                        when l2solver is 'cg', this determines how accurately the above problem is solved.
% ===============================
% Examples:
%
% d = 200;
% p = 240;
% m = 80; 
% k = 180;
% operator_type = 'random';	% random tight fram
% [x0, y, Omega, M, Lambda] = Generate_Problem(d, p, m, k, operator_type);
% 
% [xhat, Lambdahat] = ICDEL2(y, M, Omega, beta, params, xinit)
% SupportDetection(Omega(x0), Omega(xhat))
%  
%
% % Written by Heping Song, Jiangsu University, China, songhp@ujs.edu.cn
% ===============================


d = length(xinit(:));

if isa(Omega, 'function_handle')
    p = length(Omega(zeros(d,1)));
else
    p = size(Omega, 1);
end


if nargin < 7
    k = 2;
	params.num_iteration = 50;    
    params.stopping_coefficient_size = 1e-6;
    params.stopping_relative_solution_change = 1e-6;
    params.l2solver = 'pseudoinverse'; %  'cg' or  'pseudoinverse'
    params.l2_accuracy = 1e-6;
    params.stopping_cosparsity = p*0.5;
    params.max_inner_iteration = 100;
    params.noise_level = 1e-6;
%     params.stopping_lagrange_multiplier_size = 1e-1;
end

if nargin < 8
    params.num_iteration = 50;    
    params.stopping_coefficient_size = 1e-6;
    params.stopping_relative_solution_change = 1e-6;
    params.l2solver = 'pseudoinverse'; %  'cg' or  'pseudoinverse'
    params.l2_accuracy = 1e-6;
    params.stopping_cosparsity = p*0.5;
    params.max_inner_iteration = 100;
    params.noise_level = 1e-6;
%     params.stopping_lagrange_multiplier_size = 1e-1;
end


iter = 0;
Lambdahat = 1:p;	% cosupport
lagmult = 1e-4;     % initial lagrange multiplier to be used in truncated least squares problem

% initialization has no effect
% xinit = ArgminOperL2Constrained(y, M, MH, Omega, OmegaH, Lambdahat, xinit, lagmult, params);
% xinit = [M; lagmult*Omega(Lambdahat,:)]\[y; zeros(length(Lambdahat), 1)];

while iter < params.num_iteration
    
    iter = iter + 1;
    
    if isa(M, 'function_handle')
            xhat = xinit + MH(y - M(xinit));
            Omega_x = Omega(xhat);
            abscoef = abs(Omega_x);
            maxcoef = max(abscoef);
            qq = quantile(abscoef, k/p);
            to_be_removed = find(abscoef >= qq);
            Lambdahat = setdiff(1:p,to_be_removed);
            
    else
            xhat = xinit + M'*(y - M * xinit);            
            Omega_x = Omega*xhat;
            abscoef = abs(Omega_x);
            maxcoef = max(abscoef);
            qq = quantile(abscoef, k/p);
            to_be_removed = find(abscoef >= qq);
            Lambdahat = setdiff(1:p,to_be_removed);


    end
    
    [xhat, ~, lagmult] = ArgminOperL2Constrained(y, M, MH, Omega, OmegaH, Lambdahat, xinit, lagmult, params);
%     xhat = [M; lagmult*Omega(Lambdahat,:)]\[y; zeros(length(Lambdahat), 1)];
    
    
    if check_stopping_criteria(xhat, xinit, maxcoef, lagmult, Lambdahat, params)
        break;
    end
    
    xinit = xhat;
    
end    

out.Lambdahat = Lambdahat;
out.iter = iter-1;
return;   



function r = check_stopping_criteria(xhat, xinit, maxcoef, lagmult, Lambdahat, params)

    r = 0;
    
    if isfield(params, 'stopping_relative_solution_change') && norm(xhat-xinit)/norm(xhat) < params.stopping_relative_solution_change
        r = 1;
        return;
    end

    if isfield(params, 'stopping_coefficient_size') && maxcoef < params.stopping_coefficient_size
        r = 1;
        return;
    end

    if isfield(params, 'stopping_lagrange_multiplier_size') && lagmult > params.stopping_lagrange_multiplier_size
        r = 1;
        return;
    end

    if isfield(params, 'stopping_cosparsity') && length(Lambdahat) < params.stopping_cosparsity
        r = 1;
        return;
    end
    
    
%%