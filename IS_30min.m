clear;close all;clc

% load input and pack
input = load('data');

% frequency response function and response amplitude operator
input.H = 1./( -(input.w).^2*input.M + i*input.w*input.C + input.K ); % FRF (m/N)
input.RAO = input.H.*input.TF; % RAO (m/m)

% Number of MCS samples
N_T = 100000;

% set random seed
randn('state',1);

% start MCS
for i = 1:N_T
    
    % normal random numbers
    q1 = normrnd(0,1);
    q2 = normrnd(0,1);
    
    % significant wave height, Hs, and spectral peak period, Tp
    h = incdfHs(q1);
    t = incdfTp(q2,q1);
    
    % weight
    weight(i) = pdfHs(h)/pdfHsIS(h);
    
    % pack inputs
    input.Hs = h;
    input.Tp = t;
    input.short_term_var = normrnd(0,1,1,input.N*2); % short-term random variables: (1) random amplitudes and (2) phases
    
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
save(['IS_' num2str(N_T) '.mat'],'Z')

%% subfuctions
% incdfHs and incdfTp depend on a site of interest
% inverse cdf function of Hs
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

% wave spectrum depends on a site of interest
% JONSWAP specturm
function Sn = Jonswap(h,t,w)

wp = 2*pi/t; % peak frequency (rad/s)
Sn = w.^-5.*exp(-1.25*(w/wp).^-4);
Sn(1) = 0;

gamma = 3.3; % shape parameter (peak enhance coefficient)
sig = w;
sig(:) = 0.07;
sig(w>wp) = 0.09;

peak = gamma.^exp(-.5*((w/wp-1)./sig).^2);
Sn = Sn.*peak;
alpha = 5.058*(h/t^2)^2*(1-0.287*log(gamma));
Sn = alpha*9.81^2*Sn;

end

function Z = surge_max(input)

% unpack input
M = input.M; % mass
C = input.C; % damping
K = input.K; % stiffness
Hs = input.Hs; % significant wave height
Tp = input.Tp; % spectral peak period
N = input.N; % number of harnomics
dw = input.dw; % width of a frequency bin
wmin = input.wmin; % minimum frequency (rad/sec)
w = input.w; % frequency range
w_LF = input.w_LF; % frequency range for low frequency response
N_LF = input.N_LF; % number of harmonics for low frequency response
Diag_surge = input.Diag_surge; % diagonal terms in QTF
RAO = input.RAO; % RAO
short_term_var = input.short_term_var; % short-term random variables

% Jonswap spectrum
Snn = Jonswap(Hs,Tp,w);

% number of frequency bins for FFT
Nt = 125664;

% quadratic transfer function (QTF)
Hx = (-w_LF.^2*M + 1i*w_LF*C + K).^-1; % frequency response function
Hx(N) = 0;

% comlex amplitudes (Rayleigh distributed)
Re = short_term_var(1:N);
Im = short_term_var(N+1:N*2);
A = (Re+1i*Im).*sqrt(dw*Snn);

% WF response (wave frequency)
Z = [zeros(1,wmin/dw),RAO.*A]; % 0 for low freq.

% LF response (low frequency); Newman's approximation
for k = N_LF:-1:1
    A_Aconj = A(1:N-k).*conj(A(k+1:N));
    X1(k+1) = Hx(k+1)*sum(0.5*(Diag_surge(1:N-k)+Diag_surge(k+1:N)).*A_Aconj);
end

% combine WF and LF responses
Z(1:length(X1)) = Z(1:length(X1)) + 2*X1;

% IFFT for suge motion
X = -Nt*real(ifft(Z,Nt,2));

% 30-min maximum
Z = max(X(1:41888));

end

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