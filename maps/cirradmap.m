function map = cirradmap(et,etp,alpha,beta,thetk,n)
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
    gam(Jk,1)   =  imag(exp(-i.*thetk(k)).*clog((et(Jk)-beta)./((et(Jk)-alpha)*(alpha-beta))));
end
%
[mun , h ]  =  fbie(et,etp,A,gam,n,5,[],1e-14,100);
fet         = (gam+h+i.*mun)./A;
zet         = ((et-beta)./((et-alpha)*(alpha-beta))).*exp((et-alpha).*fet);
%%



%%
map.zet    =  zet;
map.fet    =  fet;
% 
end