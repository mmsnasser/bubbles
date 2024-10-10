clear
clc
addpath ../bie; addpath ../fmm; addpath ../maps;
%%
n         =   2^12
t         =   (0:2*pi/n:2*pi-2*pi/n).';
%% 
Coe       =  [
               0.00-0.60i      0.29
%                0.00+0.60i      0.29
               
              -0.60+0.30i      0.24              
              -0.60-0.30i      0.24

               0.65+0.35i      0.20
               
               0.31+0.26i      0.09
               0.29-0.24i      0.09
              -0.33-0.00i      0.09
              
              -0.24+0.26i      0.06
              -0.26-0.24i      0.06
               0.26+0.05i      0.06

               0.00+0.24i      0.04
               0.13+0.17i      0.04
              -0.13+0.16i      0.04
               0.09-0.23i      0.04
              -0.08-0.22i      0.04

              -0.17-0.03i      0.03
              -0.17+0.07i      0.03
              -0.17-0.15i      0.03
              
               0.16-0.02i      0.02
               0.16+0.08i      0.022
               0.16-0.10i      0.025
               
              -0.14+0.02i      0.01
               0.11+0.08i      0.01
               0.12-0.15i      0.012
            
              ];
rad       =   Coe(:,2);
cen       =   Coe(:,1);
theth     =   zeros(size(rad)+1);
thetv     =   pi/2*ones(size(rad)+1);
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
% 
alpha =   0.00+0.00i;
%%
figure;
hold on; box on
for k=1:m+1
    crv = et(1+(k-1)*n:k*n);crv(n+1)=crv(1);
    plot(real(crv),imag(crv),'-r','LineWidth',1.5);
end 
plot(real(alpha),imag(alpha),'dm')
axis equal
axis([-1.05  1.05  -1.05   1.05])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
% print -dpdf  fig1_im_one
%%
%  aaa
mapv  =  strslitmap(et,etp,alpha,n,thetv);
maph  =  strslitmap(et,etp,alpha,n,theth);
%
zetvo =  mapv.zet; 
zetho =  maph.zet; 
%%
% figure;
% hold on; box on
% % 
% for k=1:m+1
%     crv = zetho(1+(k-1)*n:k*n);crv(n+1)=crv(1);
%     plot(real(crv),imag(crv),'-b','LineWidth',1.5);
%     crv = zetvo(1+(k-1)*n:k*n);crv(n+1)=crv(1);
%     plot(real(crv),imag(crv),'-r','LineWidth',1.5);
% end 
% axis equal
% axis([-3  3  -3   3])
%%
% 
%     
U    =  6.0;
zetv =  zetvo;
zeth = (1-U).*zetho;
zmap =  (zetv-zeth)/U;
%%
etb   = []; etbp = [];
for k=0:m
    Jk = 1+k*n:(k+1)*n;
    zet1  =  zmap(Jk);
    zet1p =  derfft(real(zet1))+i*derfft(imag(zet1));
    etb   = [etb;zet1];
    etbp  = [etbp;zet1p];
end
% 
mapb  =  strslitmap(etb,etbp,inf,n,theth);
zet   =  mapb.zet;
fetb  =  mapb.fet;
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
%
for k=1:m+1
    crv = etb(1+(k-1)*n:k*n);crv(n+1)=crv(1);
    plot(real(crv),imag(crv),'-k','LineWidth',2);
end 
axis([-8  8  -8   8])
% 
%
set(gca,'FontSize',14)
axis square
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
%%



%% 
[xh  , yh]  =  meshgrid([-10.0:0.025:10.0],[-6.0:0.025:6.0]);
zho         =  xh+i.*yh;
[mh,nh]     =  size(zho);
%
for j=1:m+1
    Jk = 1+(j-1)*n:j*n;
    [in,on] = inpolygon(real(zho),imag(zho),real(etb(Jk)),imag(etb(Jk)));
    zho(in) = NaN+i*NaN;
    zho(on) = NaN+i*NaN;
end
zhovo  =  zho(:);
zhov   =  zhovo(abs(zhovo)>=0).';
%
whov   =  zhov+fcau(etb,etbp,fetb,zhov,n,0);
whovo  =  (NaN+i*NaN).*ones(size(zhovo));
whovo(abs(zhovo)>=0) = whov.';
wzo    =  (NaN+i*NaN).*ones(size(zho));
wzo(:) =   whovo;
%%
figure;clf
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
%
cntv   =  [-5.5,-5.0,-4.5,-4,-3.5,-2.9,-2.4,-1.9,-1.4,-1,-0.64,-0.25,...
           -0.02,0.3,0.67,1,1.45,2.05,2.5,3,3.5,4,4.5,5,5.5];
[cnt1,cnt2] = contour(real(zho),imag(zho),imag(wzo),cntv,'b','LineWidth',1.0);
% 
for k=1:m+1
    crv = etb(1+(k-1)*n:k*n);crv(n+1)=crv(1);
    plot(real(crv),imag(crv),'-k','LineWidth',2);
end 
% plot(real(zho),imag(zho),'.r','LineWidth',1.5);
% 
%
set(gca,'FontSize',14)
axis equal
axis([-10  10  -6  6])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
print -depsc spacefig25Lc6
%%
