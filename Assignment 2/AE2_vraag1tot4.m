clc; clear;  % start with a clean sheet
seed=113049;
rng(seed); % set seed
Rep = 5000;           % nr MC replications
N = 125;            % nr observations
k = 11;              % nr instruments
e_11=[1 ; zeros(10,1)];      
a = [1.5 0.7 0.5 0.3 0.15 0.07 0.03 0];
rho = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.9 0.95];
Z = normrnd(0,1,[N,k]);
PZ = Z*inv(Z'*Z)*Z';
MZ = eye(N) - PZ;

%% VRAAG 2

r_beta_grid=(0:10:300);
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

%% VRAAG 1/3

    for alpha=1:length(a)            
        pi = a(alpha)*e_11;
        for j=1:length(rho)
            sigma = [1 rho(j); rho(j) 1];
                    rng(seed); % set seed
            for i=1:Rep
                error_mat = mvnrnd([0,0],sigma,N);
                X = Z*pi + error_mat(:,2);
                Y = error_mat(:,1);
                b_2sls = inv(X'*PZ*X)*X'*PZ*Y;
                res_2sls = Y-X*b_2sls;
                var_2sls = res_2sls'*res_2sls/(N-k)*inv(X'*PZ*X);
                se_2sls = sqrt(var_2sls);
                
                error_0 = Y;
                sigma_hat_ee = (1/(N-k))*error_0'*MZ*error_0;
                sigma_hat_eV = (1/(N-k))*error_0'*MZ*X;
                sigma_hat_VV = (1/(N-k))*X'*MZ*X;
        
                rho_hat = sigma_hat_eV/sigma_hat_ee;
                pi_tilde = inv(Z'*Z)*Z'*(X-error_0*rho_hat);
                Zpi_tilde = Z*pi_tilde;
                PZpi = Zpi_tilde*inv(Zpi_tilde'*Zpi_tilde)*Zpi_tilde';

                sigma_hat_VVe = sigma_hat_VV - (sigma_hat_eV^(2)/sigma_hat_ee);
        
                r_beta(i,j) = (1/sigma_hat_VVe)*pi_tilde'*Z'*Z*pi_tilde;
                psi_1 = chi2rnd(1);
                psi_k = chi2rnd(k-1);
                
                t_2sls(i,j) = b_2sls/se_2sls;
                AR(i,j) = ((error_0'*PZ*error_0)/k)/sigma_hat_ee;
                LM(i,j) = (error_0'*PZpi*error_0)/sigma_hat_ee;
                LR(i,j) = 0.5*(k*AR(i,j) - r_beta(i,j) + sqrt((k*AR(i,j) + r_beta(i,j))^2-4*r_beta(i,j)*(k*AR(i,j)-LM(i,j))));
                if (round(r_beta(i,j)/10)+1) <= 0
                  LR_crit(i,j) = chi2inv(0.95,k-1);
                elseif (round(r_beta(i,j)/10)+1) <= 31
                  LR_crit(i,j) = LR11_crit(round(r_beta(i,j)/10)+1);
                 else LR_crit(i,j) = chi2inv(0.95,1);
                end
            end
        end
        
        rejection_freq_t(:,alpha) = mean(abs(t_2sls)>=tinv(0.975,N-1))';
        rejection_freq_AR(:,alpha) = mean(AR>=chi2inv(0.95,k)/k)';
        rejection_freq_LM(:,alpha) = mean(LM>=chi2inv(0.95,1))';
        rejection_freq_LR(:,alpha) = mean(LR>=LR_crit)';
    end
    
    


 %% PLOTS VRAAG 1 EN 3
    legendcell = {'a=1.5', 'a=0.7', 'a=0.5', 'a=0.3', 'a=0.15', 'a=0.07', 'a=0.03', 'a=0'};
    
    figure
    plot(rho,rejection_freq_t)
    legend(legendcell,'location','northwest')
    title('T-statistic')
    xlabel('rho')
    ylabel('Rejection frequency')
    
    figure
    plot(rho,rejection_freq_AR,'LineWidth',5)
    ylim([0 0.1])
    legend(legendcell,'location','northwest')
    title('AR-statistic')
    xlabel('rho')
    ylabel('Rejection frequency')

    figure
    plot(rho,rejection_freq_LM) 
    legend(legendcell,'location','northwest')
    title('LM-statistic')
     ylim([0 0.15])
    xlabel('rho')
    ylabel('Rejection frequency')

    figure
    plot(rho,rejection_freq_LR)
    legend(legendcell,'location','northwest')
    title('LR-statistic')
     ylim([0 0.15])
    xlabel('rho')
    ylabel('Rejection frequency')

    
%% PLOTS VRAAG 2 EN 4
figure
plot(r_beta_grid,LR11_crit)
hold on
plot(r_beta_grid,LR4_crit)
title('LR critical values for k=4,11')
xlabel('r(beta)')
ylabel('CV')
ylim([0,25])
legend({'k=4','k=11'})

