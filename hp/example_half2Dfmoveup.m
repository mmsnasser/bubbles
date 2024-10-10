clear
clc
addpath ../bie; addpath ../fmm; addpath ../maps;
%%
n         =   2^12
t         =   (0:2*pi/n:2*pi-2*pi/n).';
%% 
centv = [-0.88;-0.80;-0.5 ;-0.2 ; 0];
rv    = [ 0.1 ; 0.1 ; 0.1 ; 0.082 ; 0.064];
for jj=1:length(centv)
cent =centv(jj); r = rv(jj);
rad       =   [ r               ; r ];
cen       =   [-(r+0.01)+cent*i ; (r+0.01)+cent*i ];
theth     =   [ 0               ; 0    ];
thetv     =   [ pi/2            ; pi/2 ];
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
alpha = 0.5;
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

%%
U   = 2;
zetv =  zetvo;
zeth = (1-U).*zetho;
zmap{jj} =  (zetv-zeth)/U;
% 
zetb1  =  zmap{jj}(n+1:2*n);
zetb1p =  derfft(real(zetb1))+i*derfft(imag(zetb1));
Area1  = -(2*pi/n)*sum(real(zetb1).*imag(zetb1p))
zetb2  =  zmap{jj}(2*n+1:3*n);
zetb2p =  derfft(real(zetb2))+i*derfft(imag(zetb2));
Area2  = -(2*pi/n)*sum(real(zetb2).*imag(zetb2p))
%
Area = mean([Area1 Area2]);
zmap{jj} = sqrt(pi/Area)*zmap{jj};
end
%%
cc{1}=[1    0    0];
cc{2}=[0    0    0];
cc{3}=[0    0.5  0.0]; 
cc{4}=[1    0    1];
cc{5}=[0    0    1]; 
% 
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
%
for kk=1:length(centv)
    k=2;
    crv = zmap{kk}((k-1)*n+1:k*n); crv(n)=crv(1);
    plot(real(crv),imag(crv),'color',cc{kk},'LineWidth',1.5);
end
%
%
set(gca,'FontSize',14)
axis equal
lgd = legend;
lgd.NumColumns = 2;
legend('Location','north')
Leg=legend({'$r=0.1, b=-0.88$','$r=0.1, b=-0.8$','$r=0.1, b=-0.5$','$r=0.082, b=-0.2$','$r=0.064, b=0$'},...
        'Interpreter','LaTeX');
Leg.AutoUpdate = 'off';
k = 1;
crv = zmap{1}((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'k-','LineWidth',2);
for kk=1:length(centv)
    k=3;
    crv = zmap{kk}((k-1)*n+1:k*n); crv(n)=crv(1);
    plot(real(crv),imag(crv),'color',cc{kk},'LineWidth',1.5);
end
% axis square
axis([-6  6  -0.5  11.5])
% xticks([-2:1:2])
% yticks([-1:0.5:2])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
print -depsc halffig2D
%%