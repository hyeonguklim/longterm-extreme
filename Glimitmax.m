function Gmax=Glimitmax(xin)

global Diag_surge N N_LF Snn TFv dw dw_LF w w_LF RAO wmin wmax Hx

% quadratic transfer function (QTF)
Hx = (-w_LF.^2*M + i*w_LF*C + K).^-1; % frequency response function
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

for xx=160:-1:1
A_Aconj=A(1:N-xx).*conj(A(xx+1:N)); % A(1:N-xx) - r, A(xx+1:N) - s
X1(xx+1)=Hx(xx+1)*sum(0.5*(Diag_surge(1:N-xx)+Diag_surge(xx+1:N)).*A_Aconj);
end

Z(1:length(X1))=Z(1:length(X1))+2*X1;

X=-Nt*real(ifft(Z,Nt,2));

Gmax=max(X(1:41888));

%[Hs Tp]
%plot(X(1:41888))


