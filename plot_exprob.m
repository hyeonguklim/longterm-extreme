clear;clc;close all;

% load mat files
mat1 = load('mat_files\MCS_100000.mat');
mat2 = load('mat_files\IS_100000.mat');

Z1 = mat1.Z; % MCS
Z2 = mat2.Z; % IS
Weight = mat2.weight; % weight for IS

figure();hold on;
for ii = 1:10
    
    % divide into a set
    z1 = Z1( (ii-1)*10000+1:ii*10000 );
    z2 = Z2( (ii-1)*10000+1:ii*10000 );
    weight = Weight( (ii-1)*10000+1:ii*10000 );
    
    % calculate exceedance probabilities for MCS
    z1uniq = unique(z1);
    len = length(z1);
    for jj = 1:len
        idx = z1 > z1uniq(jj);
        exprob1(jj) = sum(idx)/len;
        clear idx
    end
    
    % calculate exceedance probabilities for IS
    z2uniq = unique(z2);
    len = length(z2);
    for jj = 1:len
        idx = z2 > z2uniq(jj);
        exprob2(jj) = sum(weight(idx))/len;
        clear idx
    end
    
    % plot
    plot(z1uniq,exprob1,'r-','linewidth',2);
    plot(z2uniq,exprob2,'b--','linewidth',2);
end
leg = legend('MCS','MCS-IS');legend boxoff
set(leg,'interpreter','latex')
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
print('figures\exprob','-dpng','-r0')