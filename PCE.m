clear;close all;clc

% initial parameters
po      = 2;      % order of the polynomial chaos
nrv     = 2;      % the number of random variables
Ns      = 1e5;    % the number of samples to propagate uncertainty
ef      = 10;      % multiplication factor; noe = ef*P; noe: the number of least-square fit, P: the number of basis functions

% load datasets
dataset = load('mat_files\MCS_dataset_LR.mat');

% load QoI and associated input values
QoI = dataset.G;
X(1,:) = dataset.u1;
X(2,:) = dataset.u2;

% Hermite PCE
% define and allocate symbolic random variables
for s = 1:nrv
    eval(['syms xi_' num2str(s) ' real'])
    eval(['sym_idx(s) = xi_' num2str(s) ';'])
    sym_used{s} = sym_idx(s);
end

% obtain Hermite polynomial expressions
alpha = multi_index(nrv,po);
[PhiSym,PhiPol,PhiSqNorm,P] = Hermite_PC(nrv,po,alpha);  % Hermite polynomial basis functions

noe = ef*P;% the number of experimental design; (ef) x (the number of basis functions)

for jj = 1:10
    jj
    % expansion coefficients; using linear regression (least-square fit)
    % Monte Carlo sample without replacement
    [QoI_ex, ex_idx] = datasample(QoI,noe,'Replace',false);% QoI at experimental points
    
    % input space (Gaussian)
    u_used{:,1} = X(1,ex_idx)';
    u_used{:,2} = X(2,ex_idx)';
    
    PhiEx = double(subs(PhiSym,sym_used,u_used));% Polynomials at experimental points
    c = double(lscov(PhiEx,QoI_ex')); % minimize least-square
    Y = c'*PhiSym'; % response surface by PCE

    randn('state',jj);
    % create chaos
    Q = mvnrnd([0 0],[1 0;0 1],Ns);
    for ch = 1:nrv
        chaos{:,ch} = Q(:,ch);
    end
    
    % propagate uncertainty
    QoI_PCE = double(subs(Y,sym_used,chaos));
    
    PCE = QoI_PCE;
    %  calculate Pf(G_pce> b)
    sort_PCE = unique(QoI_PCE);
    for ii = 1:length(unique(PCE))
        Nf = find(PCE>sort_PCE(ii)); % the number of failure points
        pf_PCE(:,ii) = length(Nf)/length(PCE);
    end
    
    save(['mat_files\PCE-LR-order' num2str(po) '-' num2str(noe) 'runs-' num2str(jj) 'seed_' num2str(Ns) 'propa.mat']);
end