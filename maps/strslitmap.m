function map = strslitmap(et,etp,alpha,n,thetk)
%
%
%
%
m       =   length(et)/n-1;
for k=1:m+1
    Jk = 1+(k-1)*n:k*n;
    thet(Jk,1)  =  thetk(k);
end
%
if abs(alpha)<1e6 
    A       =  exp(i*(pi/2-thet)).*(et-alpha);
    gam     =  imag(exp(-i.*thet)./(et-alpha));
    [mun,h] =  fbie(et,etp,A,gam,n,5,[],1e-14,200);
    fet     = (gam+h+i.*mun)./A;
    zet     =  1./(et-alpha)+(et-alpha).*fet;
    % 
else  
    A        =  exp(i*(pi/2-thet));
    gam      =  imag(exp(-i.*thet).*et);
    [mun,h]  =  fbie(et,etp,A,gam,n,5,[],1e-14,200);
    fet      = (gam+h+i.*mun)./A;
    zet      =  et+fet;
    %
end
% 
%
map.zet =  zet;
map.fet =  fet;
% 
end