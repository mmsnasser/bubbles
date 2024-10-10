clear
clc
addpath ../bie; addpath ../fmm; addpath ../maps;
%%
n         =   2^12
t         =   (0:2*pi/n:2*pi-2*pi/n).';
%% 
Cof       =  [ 
               0.10             0.70+0.05i
               0.10            -0.67+0.00i
               
               0.12             0.65+0.28i
               0.15             0.50+0.60i 
               0.15            -0.50+0.50i 
               0.12            -0.73+0.25i
               0.12             0.71-0.23i
               0.15             0.60-0.55i 
               0.15            -0.50-0.60i 
               0.12            -0.65-0.25i
               
               0.14             0.17+0.68i
               0.14            -0.20+0.60i 
               0.14             0.22-0.70i
               0.14            -0.12-0.60i 
               
               0.10             0.40+0.29i
               0.10            -0.30+0.29i
               0.10             0.35-0.35i
               0.10            -0.40-0.29i
               
               0.10             0.45-0.05i
               0.10            -0.45+0.12i
               0.15             0.10+0.35i
               0.15             0.00-0.30i
               
               0.09             0.27+0.10i
               0.09             0.19-0.10i
               0.10            -0.20-0.10i
               ];

rad       =   Cof(:,1);
cen       =   Cof(:,2);
theth     =   zeros(size(rad));
thetv     =   pi/2+zeros(size(rad));
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
alpha =  0.00+0.00i;
%%
figure(1);clf
clf;hold on; box on
% 
for k=1:m+1
    crv = et(1+(k-1)*n:k*n);crv(n+1)=crv(1);
    plot(real(crv),imag(crv),'-b','LineWidth',1.5);
end 
plot(real(alpha),imag(alpha),'or','LineWidth',1.5);
axis equal
axis square
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
%%



%%
maph = chanmap(et,etp,alpha,n,theth);
mapv = chanmap(et,etp,alpha,n,thetv);
%%
U = 1.5;
Ch      =  2*(1-U)/pi;
Cv      =  2/pi;
% 
zeth    =  maph.zet; zeth = Ch.*zeth;
zetv    =  mapv.zet; zetv = Cv.*zetv;
% zetp   =  map.zetp;
%
zmap = (zetv-zeth)/U;
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
axis([-2  2  -1.25   1.25])
xticks([-3:1:3])
yticks([-1.5:0.5:1.5])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
% set(gcf,'Renderer','zbuffer')
% print -depsc chanfig25Lcs
%%