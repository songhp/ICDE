%% A simple demo of Iterative Cosupport Detection Estimation (ICDE) on sparse analysis recovery
% demo_ICDEL2



close all; clear; clc;
% d = 200;
% p = 220;
% m = 80;
% k = 190;
% 
% operator_type = 'random';	% random tight fram
% [x0, y, Omega, M, Lambda] = Generate_Problem(d, p, m, k, operator_type);


load demoL2

beta = 0.5;

                   
params.num_iteration = 50;
params.stopping_relative_solution_change = 1e-6;
params.stopping_coefficient_size = 1e-6;
params.stopping_residual_size = 1e-6;
params.stopping_cosparsity = p*0.5;
params.noise_level = 1e-6;
xinit = zeros(d, 1);


iter = 0;
lagmult = 1e-4; % initial lagrange multiplier to be used in truncated least squares problem
Lambdahat = 1:p; % cosupport
stop = 0;


greedy_level = p-k;  % determines how many rows of Omega eliminates at each iteration.


% xinit = Analysis_L1(y, Omega, M, 1e-6);
% xinit = [M; lagmult*Omega(Lambdahat,:)]\[y; zeros(length(Lambdahat), 1)];

while iter < params.num_iteration
    
    iter = iter + 1;
    
	xinit =xinit +  M'*(y - M*xinit); 
%     xinit =xinit + ( M'*(y - M*xinit) )/m; 


%% method 1: eliminate decreasingly
% 	Omega_x = Omega*xinit;
% 	abscoef = abs(Omega_x(Lambdahat));
% 	maxcoef = max(abscoef);
% 	th = maxcoef* beta;
% 	to_be_removed = abscoef >= th;
% 	Lambdahat(to_be_removed) = [];

    
%% method 2: eliminate NON-decreasingly
% 	Omega_x = Omega*xinit;
% 	abscoef = abs(Omega_x);
%     maxcoef = max(abscoef);
%     qq = quantile(abscoef, 1 - greedy_level/p);
%     to_be_removed = find(abscoef >= qq);
%     Lambdahat = setdiff(1:p,to_be_removed);


%% method 3: eliminate NON-decreasingly 
    Omega_x = Omega*xinit;
    abscoef = abs(Omega_x);
    abscoefL = abscoef(Lambdahat); %
    maxcoef = max(abscoefL);
    th = maxcoef * beta;
    [to_be_removed, ~] = find(abscoef >= th);
    Lambdahat = setdiff(1:p,to_be_removed);


    

    xhat = [M; lagmult*Omega(Lambdahat,:)]\[y; zeros(length(Lambdahat), 1)];
    
    SupportDetection(Omega*x0, Omega*xhat)
    hold on
	plot([1 d],[th th],'--g','LineWidth', 2,'DisplayName','threshold');
	plot([1 d],[-th -th],'--g','LineWidth', 2,'DisplayName','threshold');
	disp('Paused. Press any key to continue...');
    grid on
    filename= ['demoICDEL20',num2str(iter)];
%     print(gcf, '-depsc2',[filename, '.eps'])
    print(gcf, '-dpng', [filename, '.png'])

	pause;    
    
	x_RelErr2 = norm(xhat-xinit)/norm(xhat);
    
	disp(['**ICDEL2 iter = ',num2str(iter) ,'  maxcoef= ', num2str(maxcoef),' cosparsity = ',num2str(length(Lambdahat)) ]);
%     disp([ 'x_RelErr2= ', num2str(x_RelErr2), '  r_RelErr2= ', num2str(r_RelErr2) ]);
    	
    if  x_RelErr2 < params.stopping_relative_solution_change...
            || maxcoef < params.stopping_coefficient_size...
            || length(Lambdahat) < 0.5*d
        break; % convergent
    end
    xinit = xhat;

end
