clear;clc;close all;

% load mat files
mat1 = load('mat_files\MCS_100000.mat');
mat2 = load('mat_files\PCE_100000.mat');

Z1 = mat1.Z; % MCS
Z2 = mat2.Z; % PCE

figure();hold on;
for ii = 1:10
    
    % divide into a set
    z1 = Z1( (ii-1)*10000+1:ii*10000 );
    z2 = Z2( (ii-1)*10000+1:ii*10000 );
    
    % calculate exceedance probabilities for MCS
    [cdf1,z1rng] = ecdf(z1);
    exprob1 = 1-cdf1;
    
    % calculate exceedance probabilities for PCE
    [cdf2,z2rng] = ecdf(z2);
    exprob2 = 1-cdf2;
    
    % plot
    plot(z1rng,exprob1,'r-','linewidth',2);
    plot(z2rng,exprob2,'b--','linewidth',2);
end
leg = legend('MCS','PCE');legend boxoff
set(leg,'interpreter','latex')
set(gca,'yscale','log',...
    'box','on',...
    'yscale','log',...
    'ticklabelinterpreter','latex',...
    'ticklength',[.02 .02],...
    'fontsize',15,...
    'xlim',[0 6],...
    'ylim',[1e-4 1])
xlabel('$z~(m)$','interpreter','latex')
ylabel('$G_z(z)$','interpreter','latex')
xtickformat('%.1f')
print('figures\exprob_pce','-dpng','-r0')