clc; clear;
seed=113049;
data = load('assignmentweakinstruments.mat');
W = [ones(3010,1) data.exper data.exper2 data.south data.smsa data.race];
MW = eye(size(W,1))-W*inv(W'*W)*W';
Y = MW*data.wage;
X = MW*data.ed;
PX = X*inv(X'*X)*X';
MX = eye(size(X,1))-PX;
Z = [data.nearc2 data.nearc4 data.nearc4a data.nearc4b];
Z = MW*Z;
PZ = Z*inv(Z'*Z)*Z';
MZ = eye(size(Z,1))-PZ;
N = size(Y,1);
k = size(Z,2);
beta_0=(-3:0.01:3);
psi_1 = chi2rnd(1);
psi_k = chi2rnd(k-1);
Rep = 5000;

%% Simulate crit values
r_beta_grid=(0:0.1:300);
for j=1:length(r_beta_grid)
 rng(seed);
for i=1:Rep
        psi_1 = chi2rnd(1);
        psi_4 = chi2rnd(4-1);
        psi_11 = chi2rnd(11-1);
        LR4(i,j) = 0.5*(psi_4 + psi_1 - r_beta_grid(j) + sqrt((psi_4 + psi_1 + r_beta_grid(j))^2-4*r_beta_grid(j)*psi_4));
        LR11(i,j) = 0.5*(psi_11 + psi_1 - r_beta_grid(j) + sqrt((psi_11 + psi_1 + r_beta_grid(j))^2-4*r_beta_grid(j)*psi_11));
    end
end

LR4_sort = sort(LR4);
LR4_crit = LR4_sort(0.95*Rep,:);

LR11_sort = sort(LR11);
LR11_crit = LR11_sort(0.95*Rep,:);

%% Determine statistics

    for j=1:length(beta_0) 
                b_2sls = inv(X'*PZ*X)*X'*PZ*Y;
                res_2sls = Y-X*b_2sls;
                var_2sls = res_2sls'*res_2sls/(N-1)*inv(X'*PZ*X);
                se_2sls = sqrt(var_2sls);
                
                error_0 = Y-X*beta_0(j);
                sigma_hat_ee = (1/(N-k))*error_0'*MZ*error_0;
                sigma_hat_eV = (1/(N-k))*error_0'*MZ*X;
                sigma_hat_VV = (1/(N-k))*X'*MZ*X;
        
                rho_hat = sigma_hat_eV/sigma_hat_ee;
                pi_tilde = inv(Z'*Z)*Z'*(X-error_0*rho_hat);
                Zpi_tilde = Z*pi_tilde;
                PZpi = Zpi_tilde*inv(Zpi_tilde'*Zpi_tilde)*Zpi_tilde';

                sigma_hat_VVe = sigma_hat_VV - (sigma_hat_eV^2/sigma_hat_ee);
        
                r_beta(j) = (1/sigma_hat_VVe)*Zpi_tilde'*Zpi_tilde;
                %psi_1 = chi2rnd(1);
                %psi_k = chi2rnd(k-1);
                
                t_2sls(j) = (b_2sls-beta_0(j))/sigma_hat_ee;
                AR(j) = ((error_0'*PZ*error_0)/k)/sigma_hat_ee;
                LM(j) = ((error_0'*PZpi*error_0))/sigma_hat_ee;
                LR(j) = 0.5*(k*AR(j) - r_beta(j) + sqrt((k*AR(j)+r_beta(j))^2 - 4*r_beta(j)*(k*AR(j)-LM(j))));
                if (round(r_beta(j),1)*10+1) <= 0
                  LR_crit(j) = chi2inv(0.95,k-1);
                elseif (round(r_beta(j),1)*10+1) <= 3001
                  LR_crit(j) = LR11_crit(round(r_beta(j),1)*10+1);
                 else LR_crit(j) = chi2inv(0.95,1);
                end
    end
    
index_AR = find(AR<=chi2inv(0.95,k)/k);
beta_conf_AR = [beta_0(min(index_AR)) beta_0(max(index_AR))];

index_t = find(abs(t_2sls)<=norminv(0.95));
beta_conf_t = [beta_0(min(index_t)) beta_0(max(index_t))];

index_LR = find(LR<=LR_crit);
beta_conf_LR = [beta_0(min(index_LR)) beta_0(max(index_LR))];

index_LM = find(LM<=chi2inv(0.95,1));
beta_conf_LM = [beta_0(min(index_LM)) beta_0(max(index_LM(index_LM<300))) beta_0(min(index_LM(index_LM>300))) beta_0(max(index_LM))];
%% PLOTS 5E
figure
plot(beta_0,t_2sls, 'blue')
hold on
plot(beta_0, ones(size(beta_0))*norminv(0.95), 'b--')
hold on
plot(beta_0,AR, 'red')
hold on
plot(beta_0,ones(size(beta_0))*(chi2inv(0.95,k)/k), 'r--')

plot(beta_0,LM, 'k')
hold on
plot(beta_0,ones(size(beta_0))*(chi2inv(0.95,1)), 'k--')
hold on
plot(beta_0,LR, 'g')
hold on
plot(beta_0,crit, 'g--')
title('Statistics for Z=nearc2 nearc4 nearc4a nearc4b')
xlabel('beta_0')
legend({'T-statistic', 'CV for T-statistic', 'AR-statistic', 'CV for AR-statistic', 'LM-statistic', 'CV for LM-statistic', 'LR-statistic', 'CV for LR-statistic'})


