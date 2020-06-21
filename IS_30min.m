clear;close all;clc
% original MCS codes provided by Prof. Ying Min Low in 2017
% modified for IS by HyeongUk Lim in 2020

% Number of MCS samples
N_T = 100000;

% set random seed
randn('state',1);

% start MCS
for i = 1:N_T
    
    % normal random numbers
    q1 = normrnd(0,1);
    q2 = normrnd(0,1);
    
    % pack inputs
    input.q1 = q1;
    input.q2 = q2;
    
    % Hs and Tp
    h = incdfHs(q1);
    t = incdfTp(q2,q1);
    input.Hs = h;
    input.Tp = t;
    
    % weight
    weight(i) = pdfHs(h)/pdfHsIS(h);
    
    % get maximum
    Z(i) = surge_max(input);
    
end

% calculate exceedance probabilities
z = unique(Z);
len = length(z);
for ii = 1:len
    idx = Z > z(ii);
    exprob(ii) = sum(weight(idx))/len;
    clear idx
end

% exceedance probability plot
figure();box on;
h1 = plot(z,exprob,'b--','linewidth',2);
xlabel('$z~(m)$','interpreter','latex')
ylabel('$G_Z(z)$','interpreter','latex')
leg1 = legend([h1],'IS');legend boxoff
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
    'ytick',[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1],...
    'ylim',[1/N_T 1],...
    'fontsize',15)
save(['mat_files\IS_' num2str(N_T) '.mat'],'Z','weight')

function pdf = pdfHs(h)

% lognormal parameters
zeta_h = 0.6565;
lamda_h = 0.77;

% Weibull parameters
gamma = 1.503;
rho = 2.691;

% pdf
eta = 2.9;
for i = 1:length(h)
    if h(i) <= eta
        pdf(i) = 1./(sqrt(2*pi)*zeta_h*h(i)).*exp(-0.5*((log(h(i))-lamda_h)/zeta_h).^2);
    else
        pdf(i) = gamma/rho*(h(i)/rho).^(gamma-1).*exp(-(h(i)/rho).^gamma);
    end
end

end

function pdf = pdfHsIS(h)

% lognormal parameters
zeta_h = 0.6565;
lamda_h = 1.0482;

% lognormal pdf
pdf = 1./(sqrt(2*pi)*zeta_h*h).*exp(-0.5*((log(h)-lamda_h)/zeta_h).^2);

end

function h = incdfHs(q1)

% importance sampling density: lognormal
% lognormal parameters
lamda_h = 1.0482;
zeta_h = 0.6565;

% normal cdf
u1 = normcdf(q1);

% transform to physical space
h = logninv(u1,lamda_h,zeta_h);

end

% inverse cdf function of Tp given Hs
function t = incdfTp(q2,q1)

% conditional lognormal distribution
h = incdfHs(q1);
lamda_th = 1.134+0.892*h.^0.225;
zeta_th = sqrt(0.005+0.12*exp(-.455*h));

% transform to physical space
u2 = normcdf(q2);
t = logninv(u2,lamda_th,zeta_th);

end