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

function Omega = Generate_Analysis_Operator(d, p,type)
%Generate operator. 
%type = 'random' for random tight fram
%type = 'TV' for total variation operator

if nargin < 3
    type = 'random';
end

if strcmp(type, 'random')
% generate random tight frame with equal column norms

if p == d
    T = randn(d);
    [Omega, discard] = qr(T);
else
    Omega = randn(p, d);
    T = zeros(p, d);
    tol = 1e-8;
    max_j = 200;
    j = 1;
    while sum(sum(abs(T-Omega))) > tol*(p*d) && j < max_j
        j = j + 1;
        T = Omega;
        [U, S, V] = svd(Omega);
        Omega = U * [eye(d); zeros(p-d,d)] * V';
        Omega = diag(1./sqrt(diag(Omega*Omega')))*Omega;
    end
    %disp(j);
end

elseif strcmp(type, 'TV')
    % The Omega matrix
    Omega=zeros(p,d);
    n = floor(sqrt(d));
    n2 = d/n;
    count=1;
    for k=1:1:n,
        for j=1:1:n2-1,
            Image=zeros(n,n2);
            Image(k,j)=1;
            Image(k,j+1)=-1;
            Omega(count,:)=Image(:)';
            count=count+1;
        end;
        Image=zeros(n,n2);
        Image(k,n2)=1;
        Image(k,1)=-1;
        Omega(count,:)=Image(:)';
        count=count+1;
    end;
    for k=1:1:n2,
        for j=1:1:n-1,
            Image=zeros(n,n2);
            Image(j,k)=1;
            Image(j+1,k)=-1;
            Omega(count,:)=Image(:)';
            count=count+1;
        end;
        Image=zeros(n,n2);
        Image(n,k)=1;
        Image(1,k)=-1;
        Omega(count,:)=Image(:)';
        count=count+1;
    end;
    Omega=Omega/sqrt(2);
elseif strcmp(type, 'TV_mult')
    Omega=zeros(p,d);
    n = sqrt(d);
    for count = 1: d-1
        Omega(count,count)=1;
        Omega(count,count+1)=-1;
    end
    Omega(d,d)=1;
    Omega(d,1)=-1;
    for count = d+1: p-1
        Omega(count,2*(count-d)-1)=1;
        Omega(count,2*(count-d))=1;
        Omega(count,2*(count-d)+1)=-1;
        Omega(count,2*(count-d)+2)=-1;        
    end
    Omega(p,d-1)=1;
    Omega(p,d)=1;
    Omega(p,1)=-1;
    Omega(p,1)=-1;        
elseif strcmp(type,'DCT')
    Omega = zeros(d,d);
    Omega(1,:) = 1/sqrt(d);
    for k = 2:d
        v = sqrt(2/d)*cos((1:2:2*d-1)*pi*(k-1)/(2*d))';
        %v = v-mean(v);
        Omega(k,:) = v';%/norm(v);
    end
else
    disp('error in type');
    disp('should be TV or random');
end