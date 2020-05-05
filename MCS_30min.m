clear;close all;clc

% load data
load data

% frequency response fuction
H = 1./( -w.^2*M+i*w*C+K );

% response amplitude operator (RAO)
RAO = H.*TF;

% Number of MCS samples
N_T = 1;

% set random seed
randn('state',1); % or rng(1)

% start MCS
for i = 1:N_T
    
    % normal random numbers
    q1 = normrnd(0,1);
    q2 = normrnd(0,1);
    
    % significant wave height, Hs, and spectral peak period, Tp
    h = incdfHs(q1);
    t = incdfTp(q2,q1);
    
    % pack inputs
    input(1,1:2) = [h t];
    input(1,3:1922) = normrnd(0,1,1,960*2);
    
    % get maximum
    G(i) = Glimitmax(input);
end

%% subfuctions
% inverse cdf function of Hs
function h = incdfHs(q1)

% lognormal parameters
lamda_h = 0.77;
zeta_h = 0.6565;

% Weibull parameters
gamma = 1.503;
rho = 2.691;

% normal cdf
u1 = normcdf(q1);

% transform to physical space
eta = 2.9; % for LONOWE distribution
threshold = logncdf(eta,lamda_h,zeta_h);
if u1 <= threshold
    h = logninv(u1,lamda_h,zeta_h);
else
    h = wblinv(u1,rho,gamma);
end

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

% JONSWAP spectrum
function Sn = Jonswap(h,t)
w = 0.2:0.00125:(1.4-0.00125);
wp = 2*pi/t;
Sn = w.^-5.*exp(-1.25*(w/wp).^-4);
Sn(1) = 0;

peaked = 3.3;
sig = w;
sig(:) = 0.07;
sig(w>wp) = 0.09;

peak = peaked.^exp(-.5*((w/wp-1)./sig).^2);
Sn = Sn.*peak;
alpha = 5.058*(h/t^2)^2*(1-0.287*log(peaked));
Sn = alpha*9.81^2*Sn;
end

function Gmax = Glimitmax(xin)
% load data
load data

% quadratic transfer function (QTF)
Hx = (-w_LF.^2*M + 1i*w_LF*C + K).^-1; % frequency response function
Hx(N) = 0;

% ------------------
% Environmental data
% ------------------
Hs=xin(1);
Tp=xin(2);

Snn=Jonswap(Hs,Tp);

N=960;
Nt=125664;

% ----------
% Parameters
% ----------

Re=xin(3:962);
Im=xin(963:1922);
A=(Re+1i*Im).*sqrt(dw*Snn);

% -----------
% WF response (wave frequency)
% -----------

Z=[zeros(1,wmin/dw),RAO(1,:).*A]; % 0 for low freq.

% -----------
% LF response (low frequency)
% -----------

for xx = 160:-1:1
    A_Aconj=A(1:N-xx).*conj(A(xx+1:N)); % A(1:N-xx) - r, A(xx+1:N) - s
    X1(xx+1)=Hx(xx+1)*sum(0.5*(Diag_surge(1:N-xx)+Diag_surge(xx+1:N)).*A_Aconj);
end

Z(1:length(X1))=Z(1:length(X1))+2*X1;

X=-Nt*real(ifft(Z,Nt,2));

Gmax=max(X(1:41888));

end
