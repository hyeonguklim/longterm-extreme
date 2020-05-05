function Hs=incdfHs(U1)

a=0.6565;
B=1.503;
eta=2.9;
rho=2.691;
deta=0.77;

cdf=normcdf(U1);
thre=logncdf(eta,deta,a);
if cdf<=thre
    Hs=logninv(cdf,deta,a);
else
    Hs=wblinv(cdf,rho,B);
end



