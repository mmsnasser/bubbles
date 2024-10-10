clear
clc
addpath ../bie; addpath ../fmm; addpath ../maps;
%%
n         =   2^10
t         =   (0:2*pi/n:2*pi-2*pi/n).';
%% 
rad       =   [ 0.1012     ; 0.1359 ];
cen       =   [-0.5831i    ;-0.8421i ];
theth     =   [ 0          ; 0    ];
thetv     =   [ pi/2       ; pi/2 ];
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
alpha = 0.85i;
mapv = halfrecmap(et,etp,alpha,n,thetv);
maph = halfrecmap(et,etp,alpha,n,theth);
%
U = 4.0;
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
zetb2  =  zmap(2*n+1:3*n);
zetb2p =  derfft(real(zetb2))+i*derfft(imag(zetb2));
Area2  = -(2*pi/n)*sum(real(zetb2).*imag(zetb2p))
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
%
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
axis([-0.5  0.5  -0.1  0.5])
xticks([-0.5:0.1:0.5])
yticks([-0.0:0.1:0.5])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
%%
[xh  , yh]  =  meshgrid([-0.5:0.001:0.5],[0.01:0.001:0.5]);
zho         =  xh+i.*yh;
[mh,nh]     =  size(zho);
%
for j=2:m+1
    Jk = 1+(j-1)*n:j*n;
    [in,on] = inpolygon(real(zho),imag(zho),real(zmap(Jk)),imag(zmap(Jk)));
    zho(in) = NaN+i*NaN;
    zho(on) = NaN+i*NaN;
end
zhovo  =  zho(:);
zhov   =  zhovo(abs(zhovo)>=0).';
% 
zhovb  =  Psiv(zhov);
% 
whov   =  Psi(zhovb)+(zhovb-alpha).*fcau(etb,etbp,fetb,zhovb)+i*fh0;
whovo  =  (NaN+i*NaN).*ones(size(zhovo));
whovo(abs(zhovo)>=0) = whov.';
wzo    =  (NaN+i*NaN).*ones(size(zho));
wzo(:) =   whovo;
%%
figure;
hold on; box on
k=1; crv = zmap(1+(k-1)*n:k*n);crv(n+1)=crv(1);
plot(real(crv),imag(crv),'-k','LineWidth',1.5);
% 
for k=2:m+1
    crv = zmap(1+(k-1)*n:k*n);crv(n+1)=crv(1);
    plot(real(crv),imag(crv),'-m','LineWidth',1.5);
end 
plot(real(zhov),imag(zhov),'.r','LineWidth',1.5);
axis equal
grid on
axis([-2  2.0  -1.0   2.0])
%%
figure;
hold on; box on
% 
for k=1:m+1
    crv = etb(1+(k-1)*n:k*n);crv(n+1)=crv(1);
    plot(real(crv),imag(crv),'-k','LineWidth',1.5);
end 
plot(real(zhovb),imag(zhovb),'.r','LineWidth',1.5);
axis equal
grid on
axis([-1.1  1.1  -1.1   1.1])
%%
figure;
hold on; box on
% 
for k=1:m+1
    crv = zet(1+(k-1)*n:k*n);crv(n+1)=crv(1);
    plot(real(crv),imag(crv),'-k','LineWidth',1.5);
end 
plot(real(whov),imag(whov),'.r','LineWidth',1.5);
axis equal
grid on
axis([-3  3  -1.0   3.0])
%%
figure(10);clf
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
% 
cntv  =  [0.02,0.045,0.05795,0.073,0.1,0.13,0.145915,0.16,...
          0.19,0.22,0.25,0.28,0.31,0.34,0.37];
[cnt1,cnt2] = contour(real(zho),imag(zho),imag(wzo),cntv,'b','LineWidth',1.0);
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
axis([-0.5  0.5  -0.05  0.4])
xticks([-0.4:0.2:0.4])
yticks([-0.0:0.2:0.4])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
print -depsc halffig2DU4LcW
%%
