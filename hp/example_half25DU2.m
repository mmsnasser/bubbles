clear
clc
addpath ../bie; addpath ../fmm; addpath ../maps;
%%
n         =   2^12
t         =   (0:2*pi/n:2*pi-2*pi/n).';
%% 
Coe       =  [
              -0.50+0.35i      0.1
              -0.20+0.40i      0.09
              -0.00+0.30i      0.1
               0.00+0.50i      0.05
               0.30+0.40i      0.1
               0.55+0.40i      0.1
 
              -0.15+0.58i      0.05
              -0.45+0.60i      0.075
               0.45+0.63i      0.045
               0.10+0.57i      0.05

              -0.40+0.73i      0.025
               0.38+0.70i      0.025
               0.20+0.63i      0.025

               0.70-0.05i      0.23
               0.00-0.50i      0.45
              -0.50-0.05i      0.15

               0.30+0.60i      0.05              
               0.15+0.47i      0.04              
              -0.15+0.10i      0.10              
              -0.30+0.25i      0.05              
              -0.35+0.50i      0.04              
              -0.75+0.20i      0.14              
              -0.30+0.60i      0.05    

               0.78+0.37i      0.10              
              -0.25+0.70i      0.025              
              ];
rad       =   Coe(:,2);
cen       =   Coe(:,1);
theth     =   zeros(size(rad));
thetv     =   pi/2*ones(size(rad));
%%
m = length(rad)
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
alpha =   0.2+0.1i;
mapv = halfrecmap(et,etp,alpha,n,thetv);
maph = halfrecmap(et,etp,alpha,n,theth);
%
%%

%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
k = 1;
crv = et((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'k-','LineWidth',1.5);
for k=2:m+1
    crv = et((k-1)*n+1:k*n); crv(n)=crv(1);
    plot(real(crv),imag(crv),'b-','LineWidth',1.5);
end
%
% plot(real(alpha),imag(alpha),'pr','MarkerFaceColor','r','MarkerSize',8);
%
set(gca,'FontSize',14)
axis square
axis([-1.05 1.05 -1.05 1.05])
xticks([-1:0.5:1])
yticks([-1:0.5:1])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
%%
zetho   =  maph.zet; 
zetvo   =  mapv.zet; 
%%
U    = 2.0;
zetv =  zetvo;
zeth = (1-U).*zetho;
zmap =  (zetv-zeth)/U;
% 
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
k = 1;
crv = zmap((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'k-','LineWidth',2);
for k=2:m+1
    crv = zmap((k-1)*n+1:k*n); crv(n)=crv(1);
    plot(real(crv),imag(crv),'b-','LineWidth',1.5);
end
%
%
set(gca,'FontSize',14)
axis equal
% axis square
axis([-5  5  -0.5   4.5])
xticks([-5:1:5])
% yticks([-1:0.5:2])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
% print -depsc halffigU2
%%

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
[xh  , yh]  =  meshgrid([-5.50:0.005:4.50],[0.01:0.005:5.0]);
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
figure(10);
clf
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
% 
cntv  =  [0.05:0.15:5];
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
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'LooseInset',get(gca,'TightInset'))
xticks([-5.5:1:4.5])
yticks([0:1:5])
axis([-5.5  4.5  -0.5   5])
print -depsc halffig25Lc2
%%
