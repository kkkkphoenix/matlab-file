function [ w0x,w0y,wx,wy,z0x,z0y,Et0x,Et0y ] = single_GB_refraction( w0x,z0x,w0y,z0y,thetai,thetat,C )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明


%入射波的参数
z=0;
zrx=w0x^2*pi*n1/lambda;
Rx=(z-z0x)*(1+(zrx/(z-z0x))^2);
wx=w0x*sqrt(1+((z-z0x)/zrx)^2);
Qi11=1/Rx-1i*lambda/pi/n1/wx^2;
zry=w0y^2*pi*n1/lambda;
Ry=(z-z0y)*(1+(zry/(z-z0y))^2);
wy=w0y*sqrt(1+((z-z0y)/zry)^2);
Qi22=1/Ry-1i*lambda/pi/n1/wy^2;
Q=[Qi11,0;0,Qi22];
% etai=0.5*atan2(z-z0x,zrx)+0.5*atan2(z-z0y,zry);
etai=0;%在入射点处，额外的相位为0
E0i=sqrt(w0x*w0y/wx/wy)*exp(-1i*k1*z0x+1i*etai);%也即入射点Eoi

C11=C(1);C22=C(4);
Qt11=sec(thetat)/k2*(k1*Qi11*cos(thetai)+C11*(k2*cos(thetat)-k1*cos(thetai)));
Qt22=Qi22+C22*(cos(thetat)-k1/k2*cos(thetai));

gx=-real(Qt11)/imag(Qt11);
Rx=1/real(Qt11);
wx=sqrt(-lambda/pi/n2/imag(Qt11));
w0x=sqrt(wx^2/(1+gx^2));
z0x=-Rx*gx^2/(1+gx^2);
gy=-real(Qt22)/imag(Qt22);
Ry=1/real(Qt22);
wy=sqrt(-lambda/pi/n2/imag(Qt22));
w0y=sqrt(wy^2/(1+gy^2));
z0y=-Ry*gy^2/(1+gy^2);
etat=0;

TTE1=2*cos(thetai)*sin(thetat)/sin(thetat+thetai);
TTM1=2*cos(thetai)*sin(thetat)/sin(thetai+thetat)/cos(thetai-thetat);

At=sqrt(w0x*w0y/wx/wy)*exp(1i*etat);
Et0x=Ex*TTE1*E0i/At;Et0y=Ey*TTM1*E0i/At;

X=0;Y=0;Z=sqrt(ss1^2+ss2^2);
qx=1/Qt(1);qy=1/Qt(4);
temp=sqrt(w0x*w0y/wx/wy)*exp(-1i*k2*(Z+0.5*(X^2/qx+Y^2/qy)+1i*etat));
H1=[Ex*Et0x*temp,Ey*Et0y*temp,0];























end

