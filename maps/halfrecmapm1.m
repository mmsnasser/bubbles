function map = halfrecmapm1(et,etp,alpha,n,thetk)
%
%
%
Psi   = @(z)(i*(2./(1+z)-1)); %i(1-z)/(1+z)
Psip  = @(z)(-2i./((1+z).^2));
Psiv  = @(w)(2./(1-i*w)-1);
Psivp = @(w)(2i./((1-i*w).^2));
%
%%
m             =   length(et)/n-1;
thet(1:n,1)   =   0;
for k=1:m
    Jk = 1+k*n:(k+1)*n;
    thet(Jk,1)  =  thetk(k);
end
%%


%%
A            =  exp(i.*(pi/2-thet)).*(et-alpha);
gam          =  imag(exp(-i.*thet).*Psi(et));
gam(1:n,1)   =  0;
%
[mun , h ]  =  fbie(et,etp,A,gam,n,5,[],1e-14,100);
fet         = (gam+h+i.*mun)./A;
h0          =  mean(h(1:n));
zet         = (Psi(et)+(et-alpha).*fet+i*h0-real(Psi(alpha)))./(imag(Psi(alpha))+h0);
%%



%%
map.zet    =  zet;
map.fet    =  fet;
map.h0     =  h0;
% 
end