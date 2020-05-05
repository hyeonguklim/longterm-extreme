function Snn=Jonswap(Hs,Tp)

% global w

w=0.2:0.00125:(1.4-0.00125);
wp=2*pi/Tp;
Snn=w.^-5.*exp(-1.25*(w/wp).^-4);
Snn(1)=0;

peaked=3.3;
sig=w;
sig(:)=0.07;
sig(w>wp)=0.09;


peak=peaked.^exp(-.5*((w/wp-1)./sig).^2);


Snn=Snn.*peak;

alpha=5.058*(Hs/Tp^2)^2*(1-0.287*log(peaked));


Snn=alpha*9.81^2*Snn;




