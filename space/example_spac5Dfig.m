clear
clc
addpath ../bie; addpath ../fmm; addpath ../maps;
%%
n         =   2^12
t         =   (0:2*pi/n:2*pi-2*pi/n).';
%% 
rad       =   [ 0.20       ; 0.20       ; 0.20       ; 0.20       ; 0.20       ];
cen       =   [ 0.30-0.40i ;-0.25-0.45i ;-0.55+0.10i ; 0.45+0.10i ; 0.00+0.00i ];
theth     =   0+zeros(length(rad)+1,1);
thetv     =   pi/2+zeros(length(rad)+1,1);
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
alpha =  0.15+0.45i;
mapv  =  strslitmap(et,etp,alpha,n,thetv);
maph  =  strslitmap(et,etp,alpha,n,theth);
%
zetvo =  mapv.zet; 
zetho =  maph.zet; 
%
U    =  2;
zetv =  zetvo;
zeth = (1-U).*zetho;
zmap =  (zetv-zeth)/U;
%%
clr=['r','g','b','k','m','c'];
clr=['b','b','b','b','b','b'];
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
for k=1:m+1
    crv = et((k-1)*n+1:k*n); crv(n)=crv(1);
    plot(real(crv),imag(crv),'color',clr(k),'LineWidth',1.5);
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
print -depsc spacfigd
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
for k=1:m+1
    crv = zetv((k-1)*n+1:k*n); crv(n)=crv(1);
    plot(real(crv),imag(crv),'color',clr(k),'LineWidth',2);
end
%
%
set(gca,'FontSize',14)
axis square
axis([-0.5 1.5 -5 7])
% xticks([-2:1:2])
% yticks([-2:1:3])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
print -depsc spacfigv
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
for k=1:m+1
    crv = zeth((k-1)*n+1:k*n); crv(n)=crv(1);
    plot(real(crv),imag(crv),'color',clr(k),'LineWidth',2);
end
%
%
set(gca,'FontSize',14)
axis square
axis([-6 5 -2 0])
% xticks([-2:1:2])
% yticks([-2:1:2])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
print -depsc spacfigh
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
for k=1:m+1
    crv = zmap((k-1)*n+1:k*n); crv(n)=crv(1);
    plot(real(crv),imag(crv),'color',clr(k),'LineWidth',1.5);
end
%
%
set(gca,'FontSize',14)
axis equal
% axis square
axis([-2.5  3.5  -2   4])
% xticks([-2:1:2])
% yticks([-1:0.5:2])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
print -depsc spacfigmap
%%