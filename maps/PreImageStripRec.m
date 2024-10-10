function map = PreImageStripRec(Lc,Lk,thetk,alphas,r,n,itol,imaxit)
%
%
%
%%
t             =   (0:2*pi/n:2*pi-2*pi/n).'; 
m             =   length(thetk)-1;
thet(1:n,1)   =   0;
for k=2:m+1
    thet(1+(k-1)*n:k*n,1)  =  thetk(k);
end
cent   =  Lc;
ak     = (1-0.5*r).*Lk;
%%

%%
et(1:n,1)   =   exp(i.*t);
etp(1:n,1)  =   i.*exp(i.*t);
et(1)=1; et(n/4+1)=i; et(n/2+1)=-1;
%%
for k=2:m+1
    xit(1+(k-1)*n:k*n,1)  = cent(k)+0.5*exp(i*thetk(k))*ak(k).*( cos(t)-i*r.*sin(t));
    xitp(1+(k-1)*n:k*n,1) =         0.5*exp(i*thetk(k))*ak(k).*(-sin(t)-i*r.*cos(t));    
end
%%
for k=2:m+1
    et(1+(k-1)*n:k*n,1)    =  Phiv(xit(1+(k-1)*n:k*n,1));
    etp(1+(k-1)*n:k*n,1)   =  Phivp(xit(1+(k-1)*n:k*n,1)).*xitp(1+(k-1)*n:k*n,1);
end
%%
alpha   =  Phiv(alphas);
%%

%%
err = inf;
itr = 0;
while (err>itol)
    itr  =itr+1;  
    %%
    for k=2:m+1
        xit(1+(k-1)*n:k*n,1)    =  cent(k)+0.5.*exp(i*thetk(k)).*ak(k).*(+cos(t)-i*r.*sin(t));
        xitp(1+(k-1)*n:k*n,1)   =          0.5.*exp(i*thetk(k)).*ak(k).*(-sin(t)-i*r.*cos(t));    
    end
    %%
    for k=2:m+1
        et(1+(k-1)*n:k*n,1)    =  Phiv(xit(1+(k-1)*n:k*n,1));
        etp(1+(k-1)*n:k*n,1)   =  Phivp(xit(1+(k-1)*n:k*n,1)).*xitp(1+(k-1)*n:k*n,1);
    end
    %%
    A            =  exp(i.*(pi/2-thet)).*(et-alpha);
    gam          =  imag(exp(-i.*thet).*Phi(et));
    gam(1:n,1)   =  0;
%     JJ = 1+n:(m+1)*n;
%     gam(JJ,1)    =  imag(exp(-i.*thet(JJ)).*Phi(et(JJ)));
%     for k=2:m+1
%         gam(1+(k-1)*n:k*n,1)   =  imag(exp(-i.*thet(1+(k-1)*n:k*n)).*(Phi(et(1+(k-1)*n:k*n))));
%     end
    %%
    [mun , h ]  =  fbie(et,etp,A,gam,n,5,[],1e-14,100);
    fmet        = (gam+h+i.*mun)./A;
    foi         =  fmet(n/4+1);
%     for k=1:m+1
%         wn(1+(k-1)*n:k*n,1) = Phi(et(1+(k-1)*n:k*n,1));
%     end
    wn      =   Phi(et)+(et-alpha).*fmet-(i-alpha)*foi;
    %%
    rotwn           =  exp(-i.*thet).*wn;
    for k=2:m+1
        wnL         =  rotwn((k-1)*n+1:k*n,1);
        centk(k,1)  =  exp(i.*thetk(k)).*((max(real(wnL))+min(real(wnL)))/2+i.*(max(imag(wnL))+min(imag(wnL)))/2);     
        radk(k,1)   =  max(real(wnL))-min(real(wnL)); 
    end
    cent   =  cent-1.0.*(centk-Lc);
    ak   =  ak-(1-0.5*r).*(radk-Lk) ;
    err    = (norm(centk-Lc,1)+norm(radk-Lk,1))/(m);
    [itr err]
    error (itr,1) = err;
    %%
    if itr>=imaxit
        'No convergence after Maximunm number of iterations'
        break;
    end
end
%%
for k=2:m+1
    xit(1+(k-1)*n:k*n,1)    =  cent(k)+0.5.*exp(i*thetk(k)).*ak(k).*(+cos(t)-i*r.*sin(t));
    xitp(1+(k-1)*n:k*n,1)   =          0.5.*exp(i*thetk(k)).*ak(k).*(-sin(t)-i*r.*cos(t));    
end
%
alphav = Phiv(cent(2:end));
%
for k=2:m+1
    et(1+(k-1)*n:k*n,1)    =  Phiv(xit(1+(k-1)*n:k*n,1));
    etp(1+(k-1)*n:k*n,1)   =  Phivp(xit(1+(k-1)*n:k*n,1)).*xitp(1+(k-1)*n:k*n,1);
end
% 
A           =  exp(i.*(pi/2-thet)).*(et-alpha);
gam         =  imag(exp(-i.*thet).*Phi(et));
gam(1:n,1)  =  0;
%
[mun , h ]  =  fbie(et,etp,A,gam,n,5,[],1e-14,100);
fmet        = (gam+h+i.*mun)./A;
foi         =  fmet(n/4+1);
% for k=1:m+1
%     zet(1+(k-1)*n:k*n,1) = Phi(et(1+(k-1)*n:k*n,1));
% end
zet     =   Phi(et)+(et-alpha).*fmet-(i-alpha)*foi;
%%


%%
map.m      =  m;
map.n      =  n;
map.et     =  et;
map.etp    =  etp;
map.zet    =  zet;
% map.zetp  =  zetp;
map.alphav =  alphav;
map.alpha  =  alpha;
map.fmet   =  fmet;
map.foi    =  foi;
% map.ak  =  ak;
% map.zphi  = (et-alpha).*fmet;
% map.zphip =  zphietp;
map.error =  error;
map.itr   =  itr;
% 
end