%% A simple demo of Iterative Cosupport Detection Estimation (ICDE) on sparse analysis recovery
% demo_ICDEL1

%% %First, set up CVX package
% cd cvx
% cvx_setup
% cd ..
% % select SeDuMi solver for optimization
% cvx_solver sedumi

close all; clear; clc;
% d = 200;
% p = 220;
% m = 60;
% k = 190;
% 
% operator_type = 'random';	% random tight fram
% [x0, y, Omega, M, Lambda] = Generate_Problem(d, p, m, k, operator_type);


%%

load demoL1

beta = 0.6;


params.num_iteration = 50;
params.stopping_relative_solution_change = 1e-6;
params.stopping_coefficient_size = 1e-6;
params.stopping_residual_size = 1e-6;
params.stopping_cosparsity = p*0.5;
params.noise_level = 1e-6;
xinit = zeros(d, 1);
    

iter = 0;
Lambdahat = 1:p; % cosupport
LC = []; % complement of cosupport in [1,p]
res_ = y;
stop = 0;

while  ~stop
    
	iter = iter + 1; 
    
	xhat = myAnalysis_L1(y, Omega, M, Lambdahat, 1e-6); 
	res = y - M*xhat;    

%% method 1: eliminate decreasingly
%     Omega_x = Omega * xhat;    
% 	abscoef = abs(Omega_x(Lambdahat));
% 	maxcoef = max(abscoef);
% 	th = maxcoef * beta;
% 	to_be_removed = abscoef >= th;
% 	Lambdahat(to_be_removed) = [];



%% method 2: eliminate NON-decreasingly 
	Omega_x = Omega * xhat;
    abscoef = abs(Omega_x);
    abscoefL = abscoef(Lambdahat); %
    maxcoef = max(abscoefL);
    th = maxcoef * beta;
    [to_be_removed, ~] = find(abscoef >= th);
    Lambdahat = setdiff(1:p,to_be_removed);
    
    
    
	SupportDetection(Omega*x0, Omega*xhat)
    hold on
	plot([1 p],[th th],'--g','LineWidth', 2,'DisplayName','threshold');
	plot([1 p],[-th -th],'--g','LineWidth', 2,'DisplayName','threshold');
    set(gca, 'FontSize', 12)
    grid on
    filename= ['demoICDEL10',num2str(iter)];
%     print(gcf, '-depsc',[filename, '.eps'])
    print(gcf, '-dpng', [filename, '.png'])


	disp('Paused. Press any key to continue...');
    
%     filename= ['demo_ICDEL1_',num2str(iter)];
%     print(gcf, '-depsc2',[filename, '.eps'])
%     print(gcf, '-dpng', [filename, '.png'])

	pause;
                    

    if iter ==1
        res_ = res;
    end
    
	x_RelErr2 = norm(xhat-xinit)/norm(xhat);
    r_RelErr2 = norm(res)/norm(res_);
    
	disp(['=== ICDEL1 iter = ',num2str(iter) ,'  maxcoef= ', num2str(maxcoef),' cosparsity = ',num2str(length(Lambdahat)) ]);
    disp([ 'x_RelErr2= ', num2str(x_RelErr2), '  r_RelErr2= ', num2str(r_RelErr2) ]);
    	
    if iter > params.num_iteration ||  x_RelErr2 < params.stopping_relative_solution_change ...
            || r_RelErr2 < params.stopping_residual_size ...
            || maxcoef < params.stopping_coefficient_size || length(Lambdahat) < 0.5*d
        stop = 1; % convergent
    end
    res_ = res;
    xinit = xhat;

end


