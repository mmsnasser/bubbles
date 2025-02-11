clear
clc
addpath ../bie; addpath ../fmm; addpath ../maps;
%%
n         =   2^14
t         =   (0:2*pi/n:2*pi-2*pi/n).';
%% 
rad       =   [ 0.42   ;   0.42 ];
cen       =   [-0.35i-0.35  ;   0.35i+0.35];
theth     =   [ 0      ;   0    ];
thetv     =   [ pi/2   ;   pi/2 ];
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
Uv   = [1.2,1.5,2,4,6,8];
for kk=1:length(Uv)
    U    = Uv(kk);  Ch = 2*(1-U)/pi;  Cv = 2/pi;
    zeth = Ch.*zetho;
    zetv = Cv.*zetvo;
    zmap{kk} =  (zetv-zeth)/U;
    % 
    zetb1  =  zmap{kk}(n+1:2*n);
    zetb1p =  derfft(real(zetb1))+i*derfft(imag(zetb1));
    Area1  = -(2*pi/n)*sum(real(zetb1).*imag(zetb1p))
    zetb2  =  zmap{kk}(2*n+1:3*n);
    zetb2p =  derfft(real(zetb2))+i*derfft(imag(zetb2));
    Area2  = -(2*pi/n)*sum(real(zetb2).*imag(zetb2p))
    %
%     Area = mean([Area1 Area2]);
%     zmap{kk} = sqrt(pi/Area)*zmap{kk};
end
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
%
kk=1;
k=2;
crv = zmap{kk}((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'-b','LineWidth',1.5);
k=3;
crv = zmap{kk}((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'-b','LineWidth',1.5);
%
%
set(gca,'FontSize',14)
axis equal
title({'$U=1.2$'},'Interpreter','LaTeX');
k = 1;
crv = zmap{1}((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'k-','LineWidth',2);
% axis square
axis([-2  4  -1.5  1.5])
% xticks([-2:1:2])
% yticks([-1:0.5:2])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
% print -depsc chanfigUDsym
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
%
kk=2;
k=2;
crv = zmap{kk}((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'-b','LineWidth',1.5);
k=3;
crv = zmap{kk}((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'-b','LineWidth',1.5);
%
%
set(gca,'FontSize',14)
axis equal
title({'$U=1.5$'},'Interpreter','LaTeX');
k = 1;
crv = zmap{1}((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'k-','LineWidth',2);
% axis square
axis([-2  4  -1.5  1.5])
% xticks([-2:1:2])
% yticks([-1:0.5:2])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
% print -depsc chanfigUDsym
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
%
kk=3;
k=2;
crv = zmap{kk}((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'-b','LineWidth',1.5);
k=3;
crv = zmap{kk}((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'-b','LineWidth',1.5);
%
%
set(gca,'FontSize',14)
axis equal
title({'$U=2$'},'Interpreter','LaTeX');
k = 1;
crv = zmap{1}((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'k-','LineWidth',2);
% axis square
axis([-2  4  -1.5  1.5])
% xticks([-2:1:2])
% yticks([-1:0.5:2])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
% print -depsc chanfigUDsym
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
%
kk=4;
k=2;
crv = zmap{kk}((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'-b','LineWidth',1.5);
k=3;
crv = zmap{kk}((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'-b','LineWidth',1.5);
%
%
set(gca,'FontSize',14)
axis equal
title({'$U=4$'},'Interpreter','LaTeX');
k = 1;
crv = zmap{1}((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'k-','LineWidth',2);
% axis square
axis([-2  4  -1.5  1.5])
% xticks([-2:1:2])
% yticks([-1:0.5:2])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
% print -depsc chanfigUDsym
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
%
kk=5;
k=2;
crv = zmap{kk}((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'-b','LineWidth',1.5);
k=3;
crv = zmap{kk}((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'-b','LineWidth',1.5);
%
%
set(gca,'FontSize',14)
axis equal
title({'$U=6$'},'Interpreter','LaTeX');
k = 1;
crv = zmap{1}((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'k-','LineWidth',2);
% axis square
axis([-2  4  -1.5  1.5])
% xticks([-2:1:2])
% yticks([-1:0.5:2])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
% print -depsc chanfigUDsym
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
%
kk=6;
k=2;
crv = zmap{kk}((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'-b','LineWidth',1.5);
k=3;
crv = zmap{kk}((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'-b','LineWidth',1.5);
%
%
set(gca,'FontSize',14)
axis equal
title({'$U=8$'},'Interpreter','LaTeX');
k = 1;
crv = zmap{1}((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'k-','LineWidth',2);
% axis square
axis([-2  4  -1.5  1.5])
% xticks([-2:1:2])
% yticks([-1:0.5:2])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
% print -depsc chanfigUDsym
%%
cc{1}=[1    0    1];
cc{2}=[0.6  0.1  0.2];
cc{3}=[0    0.5  0.0]; 
cc{4}=[0    0    0];
cc{5}=[0    0    1]; 
cc{6}=[1    0    0]; 
% 
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
%
for kk=1:length(Uv)
    k=2;
    crv = zmap{kk}((k-1)*n+1:k*n); crv(n)=crv(1);
    plot(real(crv),imag(crv),'color',cc{kk},'LineWidth',1.5);
end
%
%
set(gca,'FontSize',14)
axis equal
lgd = legend;
lgd.NumColumns = 5;
legend('Location','north')
Leg=legend({'$U=1.2$','$U=1.5$','$U=2$','$U=4$','$U=6$','$U=8$'},...
        'Interpreter','LaTeX');
Leg.AutoUpdate = 'off';
k = 1;
crv = zmap{1}((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'k-','LineWidth',2);
for kk=1:length(Uv)
    k=3;
    crv = zmap{kk}((k-1)*n+1:k*n); crv(n)=crv(1);
    plot(real(crv),imag(crv),'color',cc{kk},'LineWidth',1.5);
end
% axis square
axis([-2.5  4.5  -1.5  2.0])
% xticks([-2:1:2])
% yticks([-1:0.5:2])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
% set(gcf,'Renderer','painters')
set(gcf,'Renderer','zbuffer')
print  -depsc -r1000 chanfig3mdU3
%%


