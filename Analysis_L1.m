%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
% Analysis Greedy Like Algorithms (AGLA)
% Version 1.0 
%     
% Copyright 2012 Raja Giryes, Sangnam Nam, Michael Elad, Remi Gribonval, and Mike. E. Davies
% 
% For all details please refer to README.TXT
%
% This software is a free software distributed under the terms of the GNU 
% Public License version 3 (http://www.gnu.org/licenses/gpl.txt). You can 
% redistribute it and/or modify it under the terms of this licence, for 
% personal and non-commercial use and research purpose.                                        
%                                                                                                                                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xhat, Lambdahat] = Analysis_L1(y, Omega, M, noise_level)

%add CVX path
% addpath ./cvx/builtins
% addpath ./cvx/commands
% addpath ./cvx/functions
% addpath ./cvx/lib
% addpath ./cvx/structures
% addpath ./cvx
% addpath c:\my files\study\phd\cvx
% addpath c:\my files\study\phd\cvx\structures
% addpath c:\my files\study\phd\cvx\lib
% addpath c:\my files\study\phd\cvx\functions
% addpath c:\my files\study\phd\cvx\commands
% addpath c:\my files\study\phd\cvx\builtins

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
