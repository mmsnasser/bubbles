function map = chanmap(et,etp,alpha,n,thetk)
%
%
%
Phi   = @(z)(log((1+z)./(1-z)));
Phip  = @(z)(2./(1-z.^2));
Phiv  = @(z)(tanh(z/2));
Phivp = @(z)((2.*exp(z)./((exp(z)+1).^2)));
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
gam          =  imag(exp(-i.*thet).*Phi(et));
gam(1:n,1)   =  0;
%
[mun , h ]  =  fbie(et,etp,A,gam,n,5,[],1e-14,100);
fet         = (gam+h+i.*mun)./A;
foi         =  fet(n/4+1);
zet         =   Phi(et)+(et-alpha).*fet-(i-alpha)*foi;
%%



%%
map.zet    =  zet;
map.fet    =  fet;
map.foi    =  foi;
% 
end