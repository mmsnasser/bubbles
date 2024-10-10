clear
clc
addpath ../bie; addpath ../fmm; addpath ../maps;
%%
n         =   2^12
t         =   (0:2*pi/n:2*pi-2*pi/n).';
%% 
rad       =   [0.4 ];
cen       =   [0.0 ];
theth     =   [0      ;   0    ];
thetv     =   [pi/2   ;   pi/2 ];
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
%%
figure;
hold on; box on
for k=1:m+1
    crv = et(1+(k-1)*n:k*n);crv(n+1)=crv(1);
    plot(real(crv),imag(crv),'-r','LineWidth',1.5);
end 
axis equal
axis([-1.05  1.05  -1.05   1.05])
% print -dpdf  fig1_im_one
%%
alpha =  i*sqrt(rad);
mapv  =  strslitmap(et,etp,alpha,n,thetv);
maph  =  strslitmap(et,etp,alpha,n,theth);
%
zetvo =  mapv.zet; zetvo = zetvo-mean(zetvo);
zetho =  maph.zet; zetho = zetho-mean(zetho); 
%%
figure;
hold on; box on
% 
for k=1:m+1
    crv = zetho(1+(k-1)*n:k*n);crv(n+1)=crv(1);
    plot(real(crv),imag(crv),'-b','LineWidth',1.5);
    crv = zetvo(1+(k-1)*n:k*n);crv(n+1)=crv(1);
    plot(real(crv),imag(crv),'-r','LineWidth',1.5);
end 
axis equal
axis([-3  3  -3   3])
%%
% 
%     
U    =  4;
zetv =  zetvo;
zeth = (1-U).*zetho;
zmap =  (zetv-zeth)/U;
%% 
zeto  =  zmap(1:n);
zetop =  derfft(real(zeto))+i*derfft(imag(zeto));
Area  = -(2*pi/n)*sum(real(zeto).*imag(zetop));
%
zmap =  sqrt(pi/Area)*zmap;
zet1  =  zmap(1:n);
zet1p =  derfft(real(zet1))+i*derfft(imag(zet1));
zet2  =  zmap(n+1:2*n);
zet2p =  derfft(real(zet2))+i*derfft(imag(zet2));
%
etb   =  [zet1  ; zet2  ];
etbp  =  [zet1p ; zet2p ];
% 
mapb  =  strslitmap(etb,etbp,inf,n,theth);
zet   =  mapb.zet;
fetb  =  mapb.fet;
%%

%% 
[xh  , yh]  =  meshgrid([-4.0:0.005:4.0],[-3.0:0.005:3.0]);
zho         =  xh+i.*yh;
[mh,nh]     =  size(zho);
%
for j=1:2
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
figure(3);clf
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
%
cntv   =  [-2.5,-2.2,-1.9,-1.6,-1.3,-1,-0.7,-0.5277,-0.35,-0.00,0.35,0.528,...
            0.7,1.0,1.3,1.6,1.9,2.2,2.5];
[cnt1,cnt2] = contour(real(zho),imag(zho),imag(wzo),cntv,'b','LineWidth',1.0);
% 
for k=1:m+1
    crv = etb(1+(k-1)*n:k*n);crv(n+1)=crv(1);
    plot(real(crv),imag(crv),'-k','LineWidth',2);
end 
% plot(real(zho),imag(zho),'.r','LineWidth',1.5);
axis([-3  3  -3.0   3.0])
% 
%
set(gca,'FontSize',14)
axis square
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
print -depsc FigCrowdy4
%%
