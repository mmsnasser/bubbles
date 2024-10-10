clear
clc
addpath ../bie; addpath ../fmm; addpath ../maps;
%%
n         =   2^12
t         =   (0:2*pi/n:2*pi-2*pi/n).';
%% 
Coe       =  [
               0.00            0.25
               0.40i           0.05
               
               1i/1.575         1.25/15.75            
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
alpha =    0.5i;
%%
figure;
symcir = 0.5.*exp(i.*t);
hold on; box on
for k=1:m+1
    crv = et(1+(k-1)*n:k*n);crv(n+1)=crv(1);
    plot(real(crv),imag(crv),'-r','LineWidth',1.5);
end 
plot(real(symcir),imag(symcir),'--m')
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
U    =  3.0;
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
figure;clf
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
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
% axis([-10  10  -6  6])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
% print -depsc spacefigsym
%%
