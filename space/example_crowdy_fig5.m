clear
clc
addpath ../bie; addpath ../fmm; addpath ../maps;
%%
n         =   2^12
t         =   (0:2*pi/n:2*pi-2*pi/n).';
%% 
rad       =   [0.0725  ;  0.45   ];
cen       =   [0.1575i ; -0.405i ];
theth     =   [0       ;   0     ;   0    ];
thetv     =   [pi/2    ;   pi/2  ;   pi/2 ];
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
alpha =  0.4005i;
mapv  =  strslitmap(et,etp,alpha,n,thetv);
maph  =  strslitmap(et,etp,alpha,n,theth);
%
zetvo =  mapv.zet; 
zetho =  maph.zet; 
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
U    =  6;
zetv =  zetvo;
zeth = (1-U).*zetho;
zmap =  (zetv-zeth)/U;
%%
zet1  =  zmap(1:n);
zet1p =  derfft(real(zet1))+i*derfft(imag(zet1));
zet2  =  zmap(n+1:2*n);
zet2p =  derfft(real(zet2))+i*derfft(imag(zet2));
zet3  =  zmap(2*n+1:3*n);
zet3p =  derfft(real(zet3))+i*derfft(imag(zet3));
%
etb   =  [zet1  ; zet2  ; zet3  ];
etbp  =  [zet1p ; zet2p ; zet3p ];
% 
mapb  =  strslitmap(etb,etbp,inf,n,theth);
zet   =  mapb.zet;
fetb  =  mapb.fet;
%%


%% 
[xh  , yh]  =  meshgrid([-4.0:0.005:4.0],[-2.0:0.005:6.0]);
zho         =  xh+i.*yh;
[mh,nh]     =  size(zho);
%
for j=1:3
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
cntv   =  [-1.5,-1.15,-0.8,-0.4,0.0,0.4,0.8,0.969,1.15,1.5,...
            1.695,1.9,2.3,2.7,3.1,3.325,3.5,3.9,4.3,4.7,...
            5.1,5.5];
[cnt1,cnt2] = contour(real(zho),imag(zho),imag(wzo),cntv,'b','LineWidth',1.0);
% 
for k=1:m+1
    crv = etb(1+(k-1)*n:k*n);crv(n+1)=crv(1);
    plot(real(crv),imag(crv),'-k','LineWidth',2);
end 
% plot(real(zho),imag(zho),'.r','LineWidth',1.5);
axis([-4  4  -2.0   6.0])
% 
%
set(gca,'FontSize',14)
axis square
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
print -depsc FigCrowdy5
%%
