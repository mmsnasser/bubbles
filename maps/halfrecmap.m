function map = halfrecmap(et,etp,alpha,n,thetk)
%
%
%%
Psi   = @(z)(2./(z-i)-i); %i(i+z)/(i-z)
Psip  = @(z)(-2./((z-i).^2));
Psiv  = @(w)(2i./(i*w-1)+i);
Psivp = @(w)(2./((i*w-1).^2));
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
zet         =  Psi(et)+(et-alpha).*fet+i*h0;
%%



%%
map.zet    =  zet;
map.fet    =  fet;
map.h0     =  h0;
% 
end