%% IMPORT DATA AND SET PARAMETERS
clear;
metroboom = importdata('metroboomin_proc.csv');
metroboom = metroboom(metroboom(:,1)>0,:);
trf = metroboom(:,5);    % hours worked variable
range_trf = range(trf);
trf1mid=linspace(0,range_trf,100)';

%temp = round(metroboom(:,1)); % age variable
temp = metroboom(:,1); % age variable

temp1mid=linspace(min(temp),max(temp),100)';  % generate 100 midpoints where the regression is to be evaluated
%% KERNEL DENSITY ESTIMATION
[trfpoints,pdf_est,bndw] = npdensity_jp(trf,trf1mid,-2); % plug in bandwidth for y
[trfpoints_silver,pdf_silver,bndw_silver] = npdensity_jp(trf,trf1mid,-1); %silverman bandwidth for y

[temppoints,pdf_est_x,bndw_x] = npdensity_jp(temp,temp1mid,-2); % plug in bandwidth for y
[temppoints_silver,pdf_silver_x,bndw_silver_x] = npdensity_jp(temp,temp1mid,-1); %silverman bandwidth for y

% check that the density integrates to 1
binsize = (trfpoints(2)-trfpoints(1)); % since equally spaced 
fprintf('\n \n The density for plug-in bandwidth integrates to %g with bandwidth %g and binsize %g \n',sum(pdf_est).*binsize,bndw,binsize);
fprintf('The density for Silverman bandwidth integrates to %g with bandwidth %g and binsize %g \n ',sum(pdf_silver).*binsize,bndw_silver,binsize);

fprintf('Plug-in bandwith for the temperature variable (use for CV in regression later) is %g \n',bndw_x);
fprintf('Silverman bandwith for the temperature variable (use for CV in regression later) is %g \n',bndw_silver_x);

% Compare with Matlab's own density estimator and mle estimation
k_trf2=ksdensity(trf,trfpoints);  %matlab functiom
mle_est1 = mle(trf(trf>3000)); %MLE
mle_est2 = mle(trf(trf<=3000)); %MLE
pdf_mle = 0.5*normpdf(trf1mid,mle_est1(1),mle_est1(2)) + 0.5*normpdf(trf1mid,mle_est2(1),mle_est2(2));
%We use a bimodal distribution as the density of the traffic variable has 2
%peaks

%% KERNEL REGRESSION
% do a nonparametric regression of hrs on age
[Xpoints,m_reg,s2_reg,bndw, PV] = npregress_jp(trf,temp,temp1mid,bndw_x,0);
CI_lower_m = m_reg'-1.96*sqrt(s2_reg);
CI_upper_m = m_reg'+1.96*sqrt(s2_reg);

[Xpoints,m_reg_opt,s2_reg_opt,bndw_opt, PV_opt] = npregress_jp(trf,temp,temp1mid,0.1,0);
CI_lower_m_opt = m_reg_opt'-1.96*sqrt(s2_reg_opt);
CI_upper_m_opt = m_reg_opt'+1.96*sqrt(s2_reg_opt);
%% OLS regression
X_ols = [ones(length(temp),1) temp temp.^2 temp.^3];
Xpoints_ols = [ones(length(Xpoints),1) Xpoints Xpoints.^2 Xpoints.^3];
[b_ols,SE,res] = ols(trf, X_ols);
Y_pred = Xpoints_ols*b_ols; 
s2 = 1/(length(temp)-3)*(res'*res);
XX = inv(X_ols'*X_ols); 
for i=1:100
    Y_high(i) = Y_pred(i) + 1.96*sqrt(s2*Xpoints_ols(i,:)*XX*Xpoints_ols(i,:)');
    Y_low(i) = Y_pred(i) - 1.96*sqrt(s2*(Xpoints_ols(i,:)*XX*Xpoints_ols(i,:)'));
end

%% CROSS-VALIDATION

h=[0.1:0.1:2.5];
for i=1:length(h)
   [Xpoints,m_reg_cv,s2_reg_cv,bndw_cv, PV1] = npregress_jp(trf,temp,temp1mid,h(i),0);
   PV_vec(i) = PV1
   [Xpoints,m_reg_cv,s2_reg_cv,bndw_cv, PV2] = npregress_jp(trf,temp,temp1mid,h(i),1);
   PV_vec_end(i) = PV2
end

%% PLOTS DENSITY ESTIMATION
figure
plot(trfpoints,pdf_est,':k');
hold on
plot(trfpoints_silver,pdf_silver,'--blue');
hold on
plot(trfpoints,k_trf2,'Color','c');
hold on
plot(trf1mid,pdf_mle,'m') %check x-variable
legend('Plug-in bandwidth', 'Silverman bandwidth','Optimal matlab bandwidth','MLE')
xlabel('Traffic Volume')
ylabel('Pdf(X)')

%% PLOT REGRESSION ESTIMATES FOR PLUG-IN BANDWIDTH

figure
plot(Xpoints,m_reg','m');
hold on
plot(Xpoints,CI_lower_m,'r')
hold on
plot(Xpoints,Y_pred,'b')
hold on
plot(Xpoints, Y_low','c')
plot(Xpoints,CI_upper_m,'r')
hold on
plot(Xpoints, Y_high','c')
legend('Kernel regression', 'CI Kernel', 'OLS regression', 'CI OLS', 'Location', 'NorthWest')
xlim([243 310])
xlabel('Temperature')
ylabel('Traffic Volume')

%% PLOT REGRESSION ESTIMATES FOR "OPTIMAL" BANDWIDTH

figure
plot(Xpoints,m_reg_opt','m');
hold on
plot(Xpoints,CI_lower_m_opt,'r')
hold on
plot(Xpoints,Y_pred,'b')
hold on
plot(Xpoints, Y_low','c')
plot(Xpoints,CI_upper_m_opt,'r')
hold on
plot(Xpoints, Y_high','c')
legend('Kernel regression', 'CI Kernel', 'OLS regression', 'CI OLS', 'Location', 'NorthWest')
xlim([243 310])
xlabel('Temperature')
ylabel('Traffic Volume')

%% PLOT CROSS-VALIDATION
figure % plot cross val without endpoint bias correction
plot(h,PV_vec,'k')
hold on
plot(h,PV_vec_end,'r')
legend('PV(h)', 'PV(h), with bias correction')
xlabel('Bandwidth')
ylabel('PV')
%% FUNCTIONS
 
function [b_ols,SE,res]=ols(y,X)
% Calculates OLS coefficients, their standard errors and the residuals

N=size(X,1);
k=size(X,2);
XXi=(X'*X)\eye(k);
b_ols=XXi*X'*y;
res=y-X*b_ols;
s2=res'*res/(N-k);
SE=sqrt(s2*diag(X*XXi*X'));
end


function [Xmidpoints_used,k_pdf,bandwidth_used] = npdensity_jp(X,Xmidpoints,bandwidth)

   [nr,nc]  = size(X);
   
   X = X(:,1);         % take only first column of X
   meanX = mean(X);
   stdX  = std(X)';
   iota = ones(nr,1); 
   X = sortrows(X);   % order the observations for column 1
   
    
  %  nrbins is the number of bins/midpoints/gridpoints where the density is calculated  
    if Xmidpoints==0
        nrbins = 20;         % our DEFAULT number of bins (evaluation points) when it is NOT user defined
        Xmidpoints_used = linspace(X(floor(0.01*nr+1),1),X(floor(0.99*nr),1),nrbins)';  
        % this creates an equally spaced set of X midpoints between the 1 and 99 percentile 
    else
        [nrbins,ncbins] = size(Xmidpoints);
        Xmidpoints_used = Xmidpoints;
    end;
    
    firstb = Xmidpoints_used(1,:);
    lastb  = Xmidpoints_used(nrbins,:);
    binsize = ((lastb - firstb)/nrbins); % average binsize in this case

% SPECIFY THE BANDWITDTH(S)
delta =0.7764;   % change if kernel not normal

if bandwidth>=0; % use the bandwidth given in the proc argument
  bandwidth_used = bandwidth;
elseif bandwidth == -1
   iqr_trf = min(std(X),iqr(X)/1.349);
  bandwidth_used = 1.3643*delta*nr^(-1/5)*iqr_trf;
else  %use if plugin bandwidth is optimal
  bandwidth_used = 1.3643*delta*stdX*nr^(-1/5);
end;    % endif 

k_pdf  = zeros(nrbins,1);  %  this will contain the density estimates at the midpoints_used
                           %  bin by bin calculations since reduces the amount of workspace used vs k_pdf at each observation %
  
    for J=1:nrbins;        % for each bin       %
       Xb = Xmidpoints_used(J,1);          % one bin at a time  %
       Z = (iota*Xb - X)./bandwidth_used;
       KX = pdf('Normal',Z,0,1)/bandwidth_used;     % CHANGE if YOU WANT A DIFFERENT KERNEL  %
       k_pdf(J,1) = mean(KX);
    end    % for; 

end


function [Xmidpoints_used,m_regress, s2_regress,bandwidth_used, PV] = npregress_jp(Y,X,Xmidpoints,bandwidth, endpoint)

   [nrows,ncols]  = size(X);
   YX = horzcat(Y,X);    % keep Xi and Yi together when sorting
   YX = sortrows(YX,2);  %  sort observations according to the X variable (second column)
   
   Y = YX(:,1);
   X = YX(:,2);   % take only first column of X to determine the largest and smallest values
  
   meanX = mean(X);
   stdX  = std(X)';
   iota  = ones(nrows,1);   
    
  %  nrbins is the number of bins/midpoints/gridpoints where the density is calculated  
    if Xmidpoints==0       % then use default
        nrbins = 100;         %our default number of bins (evaluation points) when it is not user defined
        Xmidpoints_used = linspace(X(floor(0.01*nrows+1),1),X(floor(0.99*nrows),1),nrbins)';  
        % this creates an equally spaced set of X midpoints between the
        %                                         1 and 99 percentile  of X
    else
        [nrbins,ncbins] = size(Xmidpoints);
        Xmidpoints_used = Xmidpoints;
    end;
    
    firstb = Xmidpoints_used(1,:);
    lastb  = Xmidpoints_used(nrbins,:);
    binsize = ((lastb - firstb)/nrbins); % average binsize in this case
      delta =0.7764;  % see (9.11) and table 9.1 Cameron & Trivedi


if bandwidth>=0; % use the bandwidth given in the proc argument
  bandwidth_used = bandwidth;
  elseif bandwidth == -1
   iqr_trf = min(std(X),iqr(X)/1.349);
  bandwidth_used = 1.3643*delta*nrows^(-1/5)*iqr_trf;
 else 
  % SPECIFY THE BANDWITDTH(S)
  bandwidth_used = 1.3643*delta*stdX*nrows^(-1/5);
                  % use normal
                  % use Cross Validation
end;    %endif 

m_regress= zeros(nrbins,1);
s2_regress= zeros(nrbins,1);
weights= zeros(nrows,nrbins);

      %  bin by bin due to workspace limitations %
    
    for J=1:nrbins                        % for each bin       %
       Xb = Xmidpoints_used(J,1);          % one bin at a time  %
       Z = (iota*Xb - X)/bandwidth_used;
       KX = pdf('Normal',Z,0,1)/bandwidth_used;     % CHANGE if YOU WANT A DIFFERENT KERNEL  %
       YKX= Y.*KX;
       m_regress(J,1) = mean(YKX)/mean(KX);
       weights(:,J) = KX/sum(KX);
       errors(:,J) = Y-m_regress(J,1);
    end    % for; 
    s2_regress = sum(weights.^2.*errors.^2);
    
 if endpoint == 1 %set endpoint to 1 to find crossvalidations with endpoint-bias correction
  lowX = round(nrows*0.1);
  highX = round(nrows*0.9);
 X = X(lowX:highX);
  Y = Y(lowX:highX);
  nrows = length(X);
  iota = ones(length(X),1); 
 end
 
    for i=1:nrows
       Z = (X - iota*X(i))/bandwidth_used;
       KX = pdf('Normal',Z,0,1)/bandwidth_used;     % CHANGE if YOU WANT A DIFFERENT KERNEL  %
       YKX= Y.*KX;
       m_regress_2(i,1) = mean(YKX)/mean(KX);
       errors_2(i) = Y(i)-m_regress_2(i,1);
       weights_2(:,i) = KX/sum(KX);
    end
    
    weights_2 = weights_2(:,1);
    errorsquared = errors_2.^2;
    
   if endpoint == 1
       PV = sum(errorsquared'.*(1-weights_2).^2)/0.8;
   else
       PV = sum(errorsquared'.*(1-weights_2).^2);
    end
    


end
