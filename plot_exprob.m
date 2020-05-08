clear;clc;close all;

MCS = load('mat_files\MCS_100000.mat');
IS = load('mat_files\IS_100000.mat');

figure();hold on;
for ii = 1:10
    Z1 = MCS.Z( (ii-1)*10000+1:ii*10000 );
    Z2 = IS.Z( (ii-1)*10000+1:ii*10000 );
    [cdf1,z1] = ecdf(Z1);
    [cdf2,z2] = ecdf(Z2);
    plot(z1,1-cdf1,'r-','linewidth',1);
    plot(z2,1-cdf2,'b--','linewidth',1);
end
set(gca,'yscale','log',...
    'box','on',...
    'yscale','log',...
    'ticklabelinterpreter','latex',...
    'ticklength',[.02 .02],...
    'fontsize',15,...
    'ylim',[1e-4 1])
xlabel('$z~(m)$','interpreter','latex')
ylabel('$G_z(z)$','interpreter','latex')
xtickformat('%.1f')