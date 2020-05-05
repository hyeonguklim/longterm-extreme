function Tp=incdfTp(U2,U1)


Hs=incdfHs(U1);

miu=1.134+0.892*Hs.^0.225;
sig=0.005+0.12*exp(-.455*Hs);
sig=sqrt(sig);

cdf=normcdf(U2);
Tp=logninv(cdf,miu,sig);


