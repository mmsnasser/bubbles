clear
clc
addpath ../bie; addpath ../fmm; addpath ../maps;
%%
n         =   2^13
t         =   (0:2*pi/n:2*pi-2*pi/n).';
%% 
rad       =   [ 0.31       ;   0.31       ;  0.33];
cen       =   [-0.55+0.33i ;  -0.55-0.33i ;  0.05+0.18i];
theth     =   [ 0          ;   0          ;  0   ];
thetv     =   [ pi/2       ;   pi/2       ;  pi/2];
%%
m = length(rad);
% 
et(1:n,1)   =   exp(i.*t);et(1)=1;et(n/4+1)=i;et(n/2+1)=-1;
etp(1:n,1)  =   i.*exp(i.*t);
%%
k=1; Jk = 1+k*n:(k+1)*n;
et(Jk,1)    =  cen(k)+rad(k)*( 1.0*cos(t)-i*sin(t));
etp(Jk,1)   =         rad(k)*(-1.0*sin(t)-i*cos(t));
k=2; Jk = 1+k*n:(k+1)*n;
et(Jk,1)    =  cen(k)+rad(k)*( 1.0*cos(t)-1.0i*sin(t));
etp(Jk,1)   =         rad(k)*(-1.0*sin(t)-1.0i*cos(t));
k=3; Jk = 1+k*n:(k+1)*n;
et(Jk,1)    =  cen(k)+rad(k)*( 0.5*cos(t)-0.5i*sin(t));
etp(Jk,1)   =         rad(k)*(-0.5*sin(t)-0.5i*cos(t));
%%
alpha = 0.5-0.5i;
mapv = chanmap(et,etp,alpha,n,thetv);
maph = chanmap(et,etp,alpha,n,theth);
%

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

%%
U    = 2;  Ch = 2*(1-U)/pi;  Cv = 2/pi;
zeth = Ch.*zetho;
zetv = Cv.*zetvo;
zmap =  (zetv-zeth)/U;
% 
zetb1  =  zmap(n+1:2*n);
zetb1p =  derfft(real(zetb1))+i*derfft(imag(zetb1));
Area1  = -(2*pi/n)*sum(real(zetb1).*imag(zetb1p))
zetb2  =  zmap(2*n+1:3*n);
zetb2p =  derfft(real(zetb2))+i*derfft(imag(zetb2));
Area2  = -(2*pi/n)*sum(real(zetb2).*imag(zetb2p))
%
%%
% 
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
%
k=2;  crv = zmap((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'-b','LineWidth',1.5);
%
%
set(gca,'FontSize',14)
axis equal
k = 1;
crv = zmap((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'-k','LineWidth',2);
k=3; crv = zmap((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'-b','LineWidth',1.5);
k=4;
crv = zmap((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'-b','LineWidth',1.5);
% axis square
axis([-4  1  -1.25  1.25])
% xticks([-2:1:2])
% yticks([-1:0.5:2])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
print -depsc chanfig2F1B
%%