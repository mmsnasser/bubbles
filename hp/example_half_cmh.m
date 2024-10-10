clear
clc
addpath ../bie; addpath ../fmm; addpath ../maps;
%%
n         =   2^10
t         =   (0:2*pi/n:2*pi-2*pi/n).';
%% 
rad       =   [ 0.17  ;   0.22 ];
cen       =   [-0.50  ;   0.50 ];
theth     =   [ 0     ;   0    ];
%%
m = length(rad);
U = 2;
% 
et(1:n,1)   =   exp(i.*t);
etp(1:n,1)  =   i.*exp(i.*t);
%%
for k=1:m
    Jk = 1+k*n:(k+1)*n;
    et(Jk,1)    =  cen(k)+rad(k)*exp(-i*t);
    etp(Jk,1)   =      -i*rad(k)*exp(-i*t);
end
%%
alpha = 0.5i;
maph = halfrecmap(et,etp,alpha,n,theth);
%
zeth    =  maph.zet; 
h0      =  maph.h0;
fet     =  maph.fet;
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
figure;
hold on; box on
k=1; crv = zeth(1+(k-1)*n:k*n);crv(n+1)=crv(1);
plot(real(crv),imag(crv),'-k','LineWidth',1.5);
% 
for k=2:m+1
    crv = zeth(1+(k-1)*n:k*n);crv(n+1)=crv(1);
    plot(real(crv),imag(crv),'-b','LineWidth',1.5);
end 
axis equal
axis([-4  4  -1   3])
%%


%%
[xh  , yh]  =  meshgrid([-1.0:0.001:1.0],[-1.0:0.1:1.0]);
zho         =  xh+i.*yh;
[mh,nh]     =  size(zho);
zho(abs(zho)>1-1e-3)=NaN+i*NaN;
for k=1:m
    zho(abs(zho-cen(k))./rad(k)<=1+1e-3)=NaN+i*NaN;
end
z_ind    =   1;
for k=1:mh
    for j=1:nh
        if (abs(zho(k,j))>=0)
            zh(z_ind) = zho(k,j);
            z_ind  = z_ind+1;
        end
    end
end
%%
[xv  , yv]  =  meshgrid([-1.0:0.1:1.0],[-1.0:0.001:1.0]);
zvo         =  xv+i.*yv;
[mv,nv]     =  size(zvo);
zvo(abs(zvo)>=1-1e-3)  =  NaN+i*NaN;
for k=1:m
    zvo(abs(zvo-cen(k))./rad(k)<=1+1e-3)=NaN+i*NaN;
end
z_ind    =   1;
for k=1:mv
    for j=1:nv
        if (abs(zvo(k,j))>=0)
            zv(z_ind) = zvo(k,j);
            z_ind  = z_ind+1;
        end
    end
end
%%
Psi   = @(z)(2./(z-i)-i); %i(i+z)/(i-z)
% 
fzh       =  fcau (et,etp,fet,zh);
wzh       =  Psi(zh)+(zh-alpha).*fzh+i*h0;
% 
fzv       =  fcau (et,etp,fet,zv);
wzv       =  Psi(zv)+(zv-alpha).*fzv+i*h0;
%%
z_ind    =   1;
for k=1:mh
    for j=1:nh
        if (abs(zho(k,j))>=0)
            wzho(k,j)= wzh(z_ind) ;
            z_ind  = z_ind+1;
        else
            wzho(k,j) = NaN+i*NaN;            
        end
    end
end
z_ind    =   1;
for k=1:mv
    for j=1:nv
        if (abs(zvo(k,j))>=0)
            wzvo(k,j)= wzv(z_ind) ;
            z_ind  = z_ind+1;
        else
            wzvo(k,j) = NaN+i*NaN;            
        end
    end
end
%%

%%
figure;
hold on
box on
for k=1:mh
    plot(real(zho(k,:)),imag(zho(k,:)),'r')
end
for k=1:nv
    plot(real(zvo(:,k)),imag(zvo(:,k)),'b')
end
for k=1:m+1
    c_cr    =  et((k-1)*n+1:k*n,1); c_cr(n+1)  =  c_cr(1);
    plot(real(c_cr),imag(c_cr),'k')
end
axis equal
axis([-1.05  1.05  -1.05   1.05])
% print -dpdf   fig1_cm_or
% print -depsc  fig_disk1_or
%%
figure;
hold on
box on
for k=1:mh
    plot(real(wzho(k,:)),imag(wzho(k,:)),'r')
end
for k=1:nv
    plot(real(wzvo(:,k)),imag(wzvo(:,k)),'b')
end
for k=1:m+1
    c_cr    =  zeth((k-1)*n+1:k*n,1); c_cr(n+1)  =  c_cr(1);
    plot(real(c_cr),imag(c_cr),'k')
end
axis equal
axis([-4  4  -1   3])
% print -dpdf  fig1_cm_im
% print -depsc  fig_disk1_im
%%