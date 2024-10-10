clear
clc
addpath ../bie; addpath ../fmm; addpath ../maps;
%%
n         =   2^10
t         =   (0:2*pi/n:2*pi-2*pi/n).';
%% 
rad       =   [ 0.0558     ; 0.1003  ];
cen       =   [ 0.50091i   ; 0.32767i ];
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
alpha = -0.5i;
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
Uv   = [1.2,1.5,2,4,6,8];
% 
clrr  = [1    0    1
         0.6  0.1  0.2
         0    0.5  0.0
         0    0    0
         0    0    1
         1    0    0];
%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
%
for kk=1:length(Uv)
    U    = Uv(kk);
    zetv =  zetvo;
    zeth = (1-U).*zetho;
    zmap =  (zetv-zeth)/U;
    % 
    zetb1  =  zmap(n+1:2*n);
    zetb1p =  derfft(real(zetb1))+i*derfft(imag(zetb1));
    Area1  = -(2*pi/n)*sum(real(zetb1).*imag(zetb1p))
    zetb2  =  zmap(2*n+1:3*n);
    zetb2p =  derfft(real(zetb2))+i*derfft(imag(zetb2));
    Area2  = -(2*pi/n)*sum(real(zetb2).*imag(zetb2p))
    %
    k=2;
    crv = zmap((k-1)*n+1:k*n); crv(n)=crv(1);
    plot(real(crv),imag(crv),'color',clrr(kk,:),'LineWidth',1.5);
end
%
%
set(gca,'FontSize',14)
axis equal
Leg=legend({'$U=1.2$','$U=1.5$','$U=2$','$U=5$','$U=6$','$U=8$'},...
        'Interpreter','LaTeX','location','northwest');
Leg.AutoUpdate = 'off';
k = 1;
crv = zmap((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'k-','LineWidth',2);
for kk=1:length(Uv)
    U    = Uv(kk);
    zetv =  zetvo;
    zeth = (1-U).*zetho;
    zmap =  (zetv-zeth)/U;
    % 
    k=3;
    crv = zmap((k-1)*n+1:k*n); crv(n)=crv(1);
    plot(real(crv),imag(crv),'color',clrr(kk,:),'LineWidth',1.5);
end
% axis square
axis([-2.5  2.0  -0.5  4.0])
xticks([-2:1:2])
yticks([0:1:4])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
print -depsc halffig2Uv
%%
