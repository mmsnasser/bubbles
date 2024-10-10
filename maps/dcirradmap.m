function map = dcirradmap(et,etp,alpha,n,thetk)
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
    gam(Jk,1)   =  imag(exp(-i.*thetk(k)).*clog(et(Jk)-alpha));
end
%
[mun , h ]  =  fbie(et,etp,A,gam,n,5,[],1e-14,100);
fet         = (gam+h+i.*mun)./A;
c           =  exp(-mean(h(1:n)));
zet         =  c*(et-alpha).*exp((et-alpha).*fet);
%%



%%
map.zet    =  zet;
map.fet    =  fet;
map.c      =  c;
% 
end