clear
clc
addpath ../bie; addpath ../fmm; addpath ../maps;
%%
n         =   2^12
t         =   (0:2*pi/n:2*pi-2*pi/n).';
%% 
rad       =   [ 0.53578263];
cen       =   [ 0.33-0.0i ];
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
alpha = -0.5+0.00i;
%
figure;
hold on; box on
for k=1:m+1
    crv = et(1+(k-1)*n:k*n);crv(n+1)=crv(1);
    plot(real(crv),imag(crv),'-b','LineWidth',1.5);
end 
plot(real(alpha),imag(alpha),'pr')
axis equal
axis([-1.05  1.05  -1.05   1.05])
% print -dpdf  fig1_im_one
%%
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
figure(10);clf
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
% 
% 
Uv    = [1.2;1.5;2;4;6;8];
% Uv    = [2];
clrr  = [1    0    1
         0.6  0.1  0.2
         0    0.5  0.0
         0    0    0
         0    0    1
         1    0    0];
%
for kk=1:length(Uv)
%%     
U    =  Uv(kk);
zetv =  zetvo;
zeth = (1-U).*zetho;
zmap =  (zetv-zeth)/U;
% 
zet  =  zmap(1:n);
zetp =  derfft(real(zet))+i*derfft(imag(zet));
Area1 = -(2*pi/n)*sum(real(zet).*imag(zetp))
%
zet2 =  zmap(n+1:2*n);
zet2p =  derfft(real(zet2))+i*derfft(imag(zet2));
Area2 = -(2*pi/n)*sum(real(zet2).*imag(zet2p))
%
zmap =  sqrt(pi/Area1)*zmap;
zet  =  zmap(1:n);
zetp =  derfft(real(zet))+i*derfft(imag(zet));
Area = -(2*pi/n)*sum(real(zet).*imag(zetp))
%% 
% 
k=1;
crv = zmap(1+(k-1)*n:k*n);crv(n+1)=crv(1);
plot(real(crv),imag(crv),'color',clrr(kk,:),'LineWidth',2);
grid on
%
end
set(gca,'FontSize',14)
axis equal
axis([-4.5  4.55  -4   5])
% xticks([-1:0.5:1])
% yticks([-1:0.5:1])
Leg=legend({'$U=1.2$','$U=1.5$','$U=2$','$U=4$','$U=6$',...
        '$U=8$'},'Interpreter','LaTeX',...
        'location','northeast');
%
Leg.AutoUpdate = 'off';
for kk=1:length(Uv)
%     
U    =  Uv(kk);
zetv =  zetvo;
zeth = (1-U).*zetho;
zmap =  (zetv-zeth)/U;
% 
zet  =  zmap(1:n);
zetp =  derfft(real(zet))+i*derfft(imag(zet));
Area = -(2*pi/n)*sum(real(zet).*imag(zetp));
%
zmap =  sqrt(pi/Area)*zmap;
zet  =  zmap(1:n);
zetp =  derfft(real(zet))+i*derfft(imag(zet));
Area = -(2*pi/n)*sum(real(zet).*imag(zetp))
% 
% 
k=2;
crv = zmap(1+(k-1)*n:k*n);crv(n+1)=crv(1);
plot(real(crv),imag(crv),'color',clrr(kk,:),'LineWidth',2); 
%
end
%
set(gca,'FontSize',14)
axis square
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
print -depsc FigCrowdy2new
%%