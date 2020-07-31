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

function [x0,y,Omega,M,Lambda]=Generate_Problem(d, p, m, k, type, noiselevel)

% Building an analysis problem, which includes the ingredients:
%   - Omega - the analysis operator of size p*d
%   - M is anunderdetermined measurement matrix of size m*d (m<d)
%   - x0 is a vector of length d that satisfies ||Omega*x0||=p-k
%   - Lambda is the true location of these k zeros in Omega*x0
%   - we are given a measurement vector y=Mx0

if nargin==0 % default values
    d=100;
    p=110;
    m=50;
    k=85;
end;

if nargin < 6
    noiselevel = 0;
end

if nargin < 5
    type = 'random';
end


Omega=Generate_Analysis_Operator(d, p,type);


Lambda = [];
select_size = k;
while(1)
    Lambda=randperm(p);
    Lambda=sort(Lambda(1:select_size));
    cur_rank = rank(Omega(Lambda,:));
    if (cur_rank == k)
        break;
    end
    select_size = select_size + k - cur_rank;
end

if 1 % general M matrix
    M=randn(m,d);
    M_sum = sqrt(sum(M.^2,1));
    M_sum_mat = repmat(M_sum,m,1);
    M = M./M_sum_mat;
else % inpainting
    M=eye(d);
    pos=randperm(d);
    pos=pos(1:m);
    M=M(sort(pos),:);
end

% The signal is drawn at random from the null-space defined by the rows
% of the matreix Omega(Lambda,:)
if strcmp(type,'random')
    [U,D,V]=svd(Omega(Lambda,:));
    NullSpace=V(:,k+1:d);
    x0=NullSpace*randn(d-k,1); 
elseif strcmp(type,'TV') || strcmp(type,'TV_mult') || strcmp(type,'DCT')
    OmegaL=Omega(Lambda,:);
    x0=0;
    while norm(x0)< 1e-2
    x0=randn(d,1);
    x0=x0-pinv(OmegaL)*OmegaL*x0;
    end
    x0=x0/norm(x0); 
else
    disp('error in type');
    disp('should be TV or random');
end
    
y=M*x0;


t_norm = norm(y,2);
n = randn(m, 1);
y = y + noiselevel * t_norm * n / norm(n, 2);

return

