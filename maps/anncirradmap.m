function map = anncirradmap(et,etp,alpha,n,thetk)
%
%
%
%
%%
m             =   length(et)/n-1;
for k=1:m+1
    Jk = 1+(k-1)*n:k*n;
    thet(Jk,1)  =  thetk(k);
end
%%


%%
A            =  exp(i.*(pi/2-thet)).*(et-alpha);
for k=1:m+1
    Jk = 1+(k-1)*n:k*n;
    gam(Jk,1)   =  imag(exp(-i.*thetk(k)).*clog(et(Jk)/alpha));
end
%
[mun , h ]  =  fbie(et,etp,A,gam,n,5,[],1e-14,100);
fet         = (gam+h+i.*mun)./A;
h0          =  mean(h(1:n));
h1          =  mean(h(n+1:2*n));
c           =  exp(-h0);
zet         =  c*(et/alpha).*exp((et-alpha).*fet);
%%



%%
map.zet    =  zet;
map.fet    =  fet;
map.c      =  c;
map.h0     =  h0;
map.h1     =  h1;
% 
end