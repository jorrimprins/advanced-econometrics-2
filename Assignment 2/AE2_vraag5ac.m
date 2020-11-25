clc; clear;
data = load('assignmentweakinstruments.mat');
W = [ones(3010,1) data.exper data.exper2 data.south data.smsa data.race];
MW = eye(size(W,1))-W*inv(W'*W)*W';
Y = MW*data.wage;
X = MW*data.ed;
PX = X*inv(X'*X)*X';
MX = eye(size(X,1))-PX;
Z = [data.nearc2 data.nearc4 data.nearc4a data.nearc4b];
Z_select = MW*Z(:,1);
PZ_select = Z_select*inv(Z_select'*Z_select)*Z_select';
MZ_select = eye(size(Z_select,1))-PZ_select;
N = size(Y,1);
k_select = size(Z_select,2);
beta_0=(-3:0.01:3);

%% VRAAG 5A
    for j=1:length(beta_0) 
                b_2sls = inv(X'*PZ_select*X)*X'*PZ_select*Y;
                res_2sls = Y-X*b_2sls;
                var_2sls = res_2sls'*res_2sls/(N-1)*inv(X'*PZ_select*X);
                se_2sls = sqrt(var_2sls);
                
                error_0 = Y-X*beta_0(j);
                sigma_hat_ee = (1/(N-k_select))*error_0'*MZ_select*error_0;
                t_2sls(j) = (b_2sls-beta_0(j))/sigma_hat_ee;
                AR(j) = ((error_0'*PZ_select*error_0)/k_select)/sigma_hat_ee;
    end

index_AR = find(AR>=chi2inv(0.95,k_select)/k_select);
beta_conf_AR = [beta_0(min(index_AR)) beta_0(max(index_AR))];

index_t = find(abs(t_2sls)>=norminv(0.95));
beta_conf_t = [beta_0(min(index_t)) beta_0(max(index_t))];
    
    
%% 5A plots
figure
plot(beta_0,t_2sls, 'blue')
hold on
plot(beta_0, ones(size(beta_0))*norminv(0.95), 'b--')
hold on
plot(beta_0,AR, 'r')
hold on
plot(beta_0,ones(size(beta_0))*(chi2inv(0.95,k_select)/k_select), 'r--')
title('AR and T-statistic for Z=nearc2')
xlabel('beta_0')
legend({'T-statistic', 'CV for T-statistic', 'AR-statistic', 'CV for AR-statistic'})%,'location','northwest')



%% 5C
sigma_hat_VV = (1/(N-k_select))*X'*MZ_select*X;
pi_hat = inv(Z_select'*Z_select)*Z_select'*X;
F = (pi_hat'*Z_select'*Z_select*pi_hat/k_select)/sigma_hat_VV;
error_0 = Y-X*10000;
sigma_hat_ee = (1/(N-k_select))*error_0'*MZ_select*error_0;
AR_large = ((error_0'*PZ_select*error_0)/k_select)/sigma_hat_ee;
