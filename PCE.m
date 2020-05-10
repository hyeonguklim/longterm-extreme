clear;close all;clc

% initial parameters
p = 3; % PCE order
n = 2; % number of random variables
N_T = 1e5; % number of samples for uncertainty propagation
fac = 10; % multiplication factor; noe = fac*P

% Hermite PCE
alpha = multi_index(n,p);
[PhiSym,PhiPol,PhiSqNorm,P] = Hermite_PC(n,p,alpha);  % Hermite polynomial basis functions

% number of samples for PCE
N_E = fac*P;

% training samples for PCE
rng('default')
load('mat_files\MCS_100000.mat','Q1','Q2','Z');
[Z_train, idx] = datasample(Z,N_E,'Replace',false);% QoI at experimental points
Z_train = Z_train';
Q1_train = Q1(idx)';
Q2_train = Q2(idx)';

% PCE coefficients
PhiEx = double(subs(PhiSym,{'xi_1','xi_2'},{Q1_train,Q2_train}));
c = lscov(PhiEx,Z_train);

% PCE
Y = PhiSym*c; % response surface by PCE
Yfun = matlabFunction(Y);

% propagate uncertainty
randn('state',1)
Q = normrnd(0,1,[N_T,n]);
Z = Yfun(Q(:,1),Q(:,2));

save(['mat_files\PCE_' num2str(N_T) '.mat'],'Z')

%% subfunctions
% Hermite PCE
function [Psi_sim,Psi_p,PsiSqNorm,P] = Hermite_PC(M,p,alpha)
% Input:
% M = number of random variables
% p = PCE order
% alpha = multi-index sequence
%
% Output:
% Psi_s = basis functions in symbolic form
% Psi_p = basis functions in polynomial form
% PsiSqNorm = squared-norm of basis functions
% P = number of basis functions
%
% Based on:
% Felipe Uribe
% furibec@unal.edu.co
% Universidad Nacional de Colombia
% Manizales Campus

% M-dimensional Hermite polynomials
% number of basis functions
P = 1;
for s = 1:p
    P = P + (1/factorial(s))*prod(M+(0:s-1));
end

% 1-dimension
% symbolic form
syms xi;
He_s  = cell(p,1);
He_s{1} = sym(1);
He_s{2} = xi;
for j = 2:p+1
    He_s{j+1} = expand(xi*He_s{j} - (j-1)*He_s{j-1});
end
% polynomial form
He_p = cell(p,1);
He_p{1} = 1; % H_1 = 1
He_p{2} = [1 0]; % H_2 = x
for n = 2:p+1
    He_p{n+1} = [He_p{n} 0] - (n-1)*[0 0 He_p{n-1}]; % recursive formula
end

% define the number of random variables
x = cell(1,M);
H_s = cell(p,M);
H_p = cell(p,M);
for j = 1:M
    x{j} = sym(sprintf('xi_%d',j));
    for i = 1:p+1
        H_s{i,j} = subs(He_s{i},xi,x{j});
        H_p{i,j} = He_p{i};
    end
end

% M-dimension
Psi_s = cell(P,1); % symbolic
Psi_p = cell(P,1); % polynomial
for i = 2:P+1
    mult_s = 1;
    mult_p = 1;
    for j = 1:M
        mult_s = mult_s*H_s{alpha(i-1,j)+1,j};
        mult_p = conv(mult_p,H_p{alpha(i-1,j)+1,j});
    end
    Psi_s{i-1} = mult_s;
    Psi_p{i-1} = mult_p;
end

% squared-norm
PsiSqNorm = prod(factorial(alpha),2);
for k = 1:P
    Psi_sim(k) = (Psi_s{k});
end

end

% multi-indices
function alpha = multi_index(M,p)
% Input:
% M = number of random variables
% p = PCE order
%
% Output:
% alpha = multi-index sequence
%
% Based on:
% Felipe Uribe
% furibec@unal.edu.co
% Universidad Nacional de Colombia
% Manizales Campus

% multi-index sequence
alpha = cell(p+1,1); % multi-index
alpha{1} = zeros(1,M); % multi-index for length 0

switch M
    case 1 % dimension = 1
        for q = 1:p
            alpha{q+1} = q;
        end
    otherwise % dimension > 1
        for q = 1:p
            s = nchoosek(1:M+q-1,M-1);
            s1 = zeros(size(s,1),1);
            s2 = (M+q)+s1;
            alpha{q+1} = flipud(diff([s1 s s2],1,2))-1; % -1 due to MATLAB indexing
            if sum(alpha{q+1},2) ~= q*ones(nchoosek(M+q-1,M-1),1)
                error('The sum of each row has to be equal to q-th order');
            end
        end
end

alpha = cell2mat(alpha);

end