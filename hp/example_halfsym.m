clear
clc
addpath ../bie; addpath ../fmm; addpath ../maps;
%%
n         =   2^10
t         =   (0:2*pi/n:2*pi-2*pi/n).';
%% 
rad       =   [0.5;0.1];
cen       =   [0; 0.8i ];
theth     =   [ 0  ; 0  ];
thetv     =   [ pi/2 ; pi/2];
%%
m = length(rad);
% 
et(1:n,1)   =   exp(i.*t);et(1)=1;et(n/4+1)=i;et(n/2+1)=-1;
etp(1:n,1)  =   i.*exp(i.*t);
%%
for k=1:m
    Jk = 1+k*n:(k+1)*n;
    et(Jk,1)    =  cen(k)+rad(k)*exp(-i*t);
    etp(Jk,1)   =      -i*rad(k)*exp(-i*t);
end
%%
figure;
hold on; box on
k=1; crv = et(1+(k-1)*n:k*n);crv(n+1)=crv(1);
plot(real(crv),imag(crv),'-k','LineWidth',1.5);
% 
for k=2:m+1
    crv = et(1+(k-1)*n:k*n);crv(n+1)=crv(1);
    plot(real(crv),imag(crv),'-m','LineWidth',1.5);
end 
axis equal
grid on
axis([-1  1.0  -1.0   1.0])
%%
alpha =  0.7;
mapv = halfrecmap(et,etp,alpha,n,thetv);
maph = halfrecmap(et,etp,alpha,n,theth);
%
U = 3.0;
% 
zetho   =  maph.zet; 
zetvo   =  mapv.zet; 
% 
zetv =  zetvo;
zeth = (1-U).*zetho;
zmap =  (zetv-zeth)/U;
%%
zetb1  =  zmap(n+1:2*n);
zetb1p =  derfft(real(zetb1))+i*derfft(imag(zetb1));
Area1  = -(2*pi/n)*sum(real(zetb1).*imag(zetb1p))
%%
Psi   = @(z)(2./(z-i)-i); %i(i+z)/(i-z)
Psiv  = @(w)(2i./(i*w-1)+i);
%
zmapb = Psiv(zmap);
% 
etb   = []; etbp = [];
for k=0:m
    Jk = 1+k*n:(k+1)*n;
    zet1  =  zmapb(Jk);
    zet1p =  derfft(real(zet1))+i*derfft(imag(zet1));
    etb   = [etb;zet1];
    etbp  = [etbp;zet1p];
end
% 
mapb  = halfrecmap(etb,etbp,alpha,n,theth);
zet   =  mapb.zet;
fetb  =  mapb.fet;
fh0   =  mapb.h0;
%%




%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
% 
% 
k=1; crv = zmap(1+(k-1)*n:k*n);crv(n+1)=crv(1);
plot(real(crv),imag(crv),'-k','LineWidth',2);
% 
for k=2:m+1
    crv = zmap(1+(k-1)*n:k*n);crv(n+1)=crv(1);
    plot(real(crv),imag(crv),'-k','LineWidth',2);
end 
set(gca,'FontSize',14)
axis equal
axis([-15 15 -1 20])
% xticks([-2:1:2])
% yticks([0:1:4])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
% print -depsc halffigsym
%%
