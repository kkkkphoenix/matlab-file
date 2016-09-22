%单个高斯波束传播
c=3e8;f=3e10;
lambda=1;
w0x=5;w0y=w0x;
z0x=10;z0y=z0x;
n1=1;n2=2;
R1=60;R2=R1;
C11=1/R1;C22=1/R2;
k1=2*pi/lambda;k2=n2/n1*k1;
Ex=10;Ey=0;
theta=(-90:1:90)*pi/180;
distance=30*lambda;

%入射波的参数
z =0;
zrx=w0x^2*pi*n1/lambda;
Rx=(z-z0x)*(1+(zrx/(z-z0x))^2);
wx=w0x*sqrt(1+((z-z0x)/zrx)^2);
Qi11=1/Rx-1i*lambda/pi/n1/wx^2;
zry=w0y^2*pi*n1/lambda;
Ry=(z-z0y)*(1+(zry/(z-z0y))^2);
wy=w0y*sqrt(1+((z-z0y)/zry)^2);
Qi22=1/Ry-1i*lambda/pi/n1/wy^2;
Q=[Qi11,0;0,Qi22];
etai=0.5*atan2(z-z0x,zrx)+0.5*atan2(z-z0y,zry);
E0i=sqrt(w0x*w0y/wx/wy)*exp(-1i*k1*z0x+1i*etai);

%坐标系的定义
x_hat=[1,0,0];
y_hat=[0,1,0];
z_hat=[0,0,1];
thetai=30*pi/180;
thetat=asin(n1*thetai/n2);
phi=0;
%这里采用的是另外的定义方式，因为原文的方法在theta=0处不连续
zi_hat=[sin(thetai)*cos(phi),sin(thetai)*sin(phi),cos(thetai)];
yi_hat=cross(zi_hat,x_hat);
yi_hat=yi_hat/norm(yi_hat);
xi_hat=cross(yi_hat,zi_hat);
xi_hat=xi_hat/norm(xi_hat);

zt_hat=[sin(thetat)*cos(phi),sin(thetat)*sin(phi),cos(thetat)];
yt_hat=cross(zt_hat,x_hat);
yt_hat=yt_hat/norm(yt_hat);
xt_hat=cross(yt_hat,zt_hat);
xt_hat=xt_hat/norm(xt_hat);

c11=dot(xt_hat,x_hat);c12=dot(yt_hat,x_hat);c13=dot(zt_hat,x_hat);
c21=dot(xt_hat,y_hat);c22=dot(yt_hat,y_hat);c23=dot(zt_hat,y_hat);
c31=dot(xt_hat,z_hat);c32=dot(yt_hat,z_hat);c33=dot(zt_hat,z_hat);
C_refract2principle=[c11,c12,c13;c21,c22,c23;c31,c32,c33];

%反射波和折射波的曲率矩阵
Qr11=Qi11-2*C11/cos(thetai);
Qr22=Qi22+2*C22*cos(thetai);
Qr=[Qr11,0;0,Qr22];
Qt11=sec(thetat)/k2*(k1*Qi11*cos(thetai)+C11*(k2*cos(thetat)-k1*cos(thetai)));
Qt22=Qi22+C22*(cos(thetat)-k1/k2*cos(thetai));
Qt=[Qt11,0;0,Qt22];

%引入参数g使得计算简单，无实际意义
gx=-real(Qt11)/imag(Qt11);
Rx=1/real(Qt11);
wx=sqrt(-lambda/pi/n2/imag(Qt11));
w0x=sqrt(wx^2/(1+gx^2));
z0x=-Rx*gx^2/(1+gx^2);
zrx=w0x^2*n2*pi/lambda;
gy=-real(Qt22)/imag(Qt22);
Ry=1/real(Qt22);
wy=sqrt(-lambda/pi/n2/imag(Qt22));
w0y=sqrt(wy^2/(1+gy^2));
z0y=-Ry*gy^2/(1+gy^2);
zry=w0y^2*n2*pi/lambda;
etat=0.5*atan2(z-z0x,zrx)+0.5*atan2(z-z0y,zry);

%计算幅度
RTE=sin(thetai-thetat)/sin(thetat+thetai);
TTE=2*cos(thetai)*sin(thetat)/sin(thetat+thetai);
TTM=2*cos(thetai)*sin(thetat)/sin(thetai+thetat)/cos(thetai-thetat);
At=sqrt(w0x*w0y/wx/wy)*exp(1i*etat);
Et0x=Ex*TTE*E0i/At;
% Et0y=Ey*TTE*E0i/At;
H=zeros(length(theta),3);

for kk=1:length(theta)
    Y=0;
    Z=distance*cos(theta(kk));
    X=distance*sin(theta(kk));
    point=[X,Y,Z];
%     pp=point*C_refract2principle;%main coordinate
    Xt=X;Yt=Y;Zt=Z;
    Rx=(Zt-z0x)*(1+(zrx/(Zt-z0x))^2);
    wx=w0x*sqrt(1+((Zt-z0x)/zry)^2);
    qx=1/(1/Rx-1i*lambda/pi/n2/wx^2);
    Ry=(Zt-z0y)*(1+(zry/(Zt-z0y))^2);
    wy=w0y*sqrt(1+((Zt-z0y)/zry)^2);
    qy=1/(1/Ry-1i*lambda/pi/n2/wy^2);
    temp=sqrt(w0x*w0y/wx/wy)*exp(-1i*k2*(Zt+0.5*(Xt^2/qx+Yt^2/qy)+1i*etat));
    H(kk,:)=[Ex*Et0x*temp,Ey*Et0y*temp,0];
end
plot(theta*180/pi,(abs(H(:,1))/100000),'--r')
hold on
plot(theta*180/pi,(abs(H(:,1))/TTE/100000),'-b')
hold on
plot(theta*180/pi,(abs(H(:,1))/TTE*RTE/100000),'-g')
title('X-direction ')
xlabel('theta(degrees)')
ylabel('amplitude')
legend('transimitted','incident','reflected')
% figure(2)
% plot(theta*180/pi,(abs(H(:,2))),'--r')
% title('Y-direction')
% xlabel('theta(degrees)')
% ylabel('amplitude')