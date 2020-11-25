% Name student 1: Vien Dinh
% Name student 2: Derek Dirks
% Name student 3: Jorrim Prins
% We (enlited above) declare that
% 1. Our assignment will be our own work.
% 2. We shall not make solutions to the assignment available to anyone else.
% 3. We shall not engage in any other activities that will dishonestly improve my results or dishonestly improve or hurt the results of others.

clear; clc;                               % clear memory & screen
rng(142876);
data=importdata('cholesterol_nohdr.txt'); % read ASCII-file
BOOTREP=4000;                             % number of bootstrap replications
n      =size(data,1);                     % sample size
const  =ones(n,1);                        % constant
x=data(:,1); y=data(:,2);
sig=22;                                   % sigma
Xall=[const x x.^2 x.^3 x.^4 x.^5 x.^6];  % regressor matrix upto 6-th order

%% Introducing vectors etc
Cp_error = Cp(y,Xall,sig); %Calculate Cp_error based on different polynomials

c = -2.2:0.05:2;
c_mat = [ones(1,85) ; c ; c.^2 ; c.^3 ; c.^4 ; c.^5 ; c.^6]'; %Create c_matrix for graphs

N_matrix = zeros(BOOTREP,n);

%% Bootstrap replications

for r=1:BOOTREP
    %Resample x and y
   index=unidrnd(n,n,1);  % select the indices  
   XBall=Xall(index,:);  % resample from data
   yB = y(index); %resample from data
   N_matrix(r,:) = histc(index, 1:n);

   X3 = XBall(:,1:4);
   b3 = X3\yB;
   
   %Calculate Cp_errors per bootstrap sample for adapted method
   Cp_B_error(r,:) = Cp(yB,XBall,sig);
   [min_Cp, min_loc] = min(Cp_B_error(r,:));
   Xadapt = XBall(:,1:min_loc);
   badapt = Xadapt\yB;
   
   %Calculate predictions for different values of c
   for i = 1:length(c_mat)
   theta_hat_fixed(r,i) = c_mat(i,1:4)*b3;
   theta_hat_adapt(r,i) = c_mat(i,1:min_loc)*badapt;
   end

end


%% Frequency of polynomials used (proportion)
[val,min_loc]= min(Cp_B_error');
proportions = [(1:6)', histc(min_loc(:)-1, 1:6)/4000];

%% Calculate standard errors

%Create standard errors of fixed and adapted methods
B_se_fixed_data = std(theta_hat_fixed);
B_se_adapt_data = std(theta_hat_adapt);

Bcov = zeros(length(c),n);
%Bcov = zeros(1,n);
%Calculate standard error of smooth method by creating for loop for
%covariance values
for i = 1:length(c)
    for j = 1:n
        Bcov_single = cov(N_matrix(:,j),theta_hat_adapt(:,i));
        Bcov(i,j) = Bcov_single(1,2);
    end
end

B_se_smooth_data = sqrt(sum(Bcov.^2,2))';

%% Output Question 1
fprintf('Question 1')
table((0:6)',Cp_error-80000,'VariableNames',{'Degree','Cp'})
table((1:6)',proportions(:,2),'VariableNames',{'M','Proportion'})

%% Plot Question 2
%Plot histogram of adapted and fixed method predictions
figure
histogram(theta_hat_fixed(:,5),40,'BinLimits',[-20,20],'DisplayStyle','stairs','LineWidth',0.99,'EdgeColor','black')

hold on
histogram(theta_hat_adapt(:,5),40,'BinLimits',[-20,20],'FaceColor','cyan')
legend('Fixed degree (3)', 'Adaptive degree','Location','northeast')
xlabel('Estimated cholesterol decrease')
ylabel('Frequency')

%% Plot Question 3
%Plot standard errors
figure
plot(c,B_se_fixed_data,'--','Color','red','LineWidth',2)
xlim([-2.5 2.25])
ylim([0 7])
hold on
plot(c,B_se_adapt_data,'Color','black','LineWidth',2)
hold on
plot(c,B_se_smooth_data,':','Color','blue','LineWidth',2)
legend('Fixed', 'Adaptive','Smoothed','Location','southeast')
xlabel('Adjusted compliance c')
ylabel('Standard errors of theta hat')


%% Output Question 4
fprintf('Question 4 \n')
fprintf('(c=2) Standard error of fixed approach %4.5f \n', B_se_fixed_data(5))
fprintf('(c=2) Standard error of adaptive approach %4.5f \n', B_se_adapt_data(5))

fprintf('Standard error of fixed approach %4.5f \n', mean(B_se_fixed_data))
fprintf('Standard error of adaptive approach %4.5f \n', mean(B_se_adapt_data))
fprintf('Standard error of smoothed approach %4.5f \n', mean(B_se_smooth_data))

fprintf('Adaptive/fixed ratio %4.5f \n', mean(B_se_adapt_data)/mean(B_se_fixed_data))
fprintf('Adaptive/smooth ratio %4.5f', mean(B_se_adapt_data)/mean(B_se_smooth_data))

%% Functions

function [Cp_error] = Cp(y,Xall,sig)

Cp_error = zeros(7,1);

for i = 1:7
    X_select = Xall(:,1:i);
    b=X_select\y;
    Cp_error(i,1) = sum((y - X_select*b).^2) + 2*i*sig^2;
end

end