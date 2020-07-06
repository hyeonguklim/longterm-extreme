clear;close all;clc
% original codes provided by Prof. Ying Min Low in 2017
% modified by HyeongUk Lim in 2020

% Number of MCS samples
N_T = 1e5;

% set random seed
randn('state',1);

% start MCS
for i = 1:N_T
    
    % normal random numbers (Hs and Tp in normal space)
    Q1 = normrnd(0,1);
    Q2 = normrnd(0,1);   
    
    % get maximum
    Z(i) = surge_max_mex(Q1,Q2,0);
    
end

% get empirical cdf
[cdf,z] = ecdf(Z);
exprob = 1 - cdf;

% exceedance probability plot
figure();box on;
h1 = plot(z,exprob,'r-','linewidth',2);
xlabel('$z~(m)$','interpreter','latex')
ylabel('$G_Z(z)$','interpreter','latex')
leg1 = legend([h1],'MCS');legend boxoff
set(leg1,'interpreter','latex')
xtickformat('%.1f')
set(gca,'yscale','log',...
    'ticklabelinterpreter','latex',...
    'tickdir','in', ...
    'ticklength',[.02 .02],...
    'xminortick','off',...
    'yminortick','on',...
    'xgrid','off',...
    'ygrid','off',...
    'ytick',[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1] , ...
    'fontsize',15)
save(['mat_files\MCS_' num2str(N_T) '.mat'],'Z','Q1','Q2')
