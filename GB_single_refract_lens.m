%------------------------------------------------------------%
%�˶δ������������Ƶ�͸��������
%�ɵ����ĵ�����˹���������д�������д����Ľ��ƺ������������Ĵ���
%���Ҫ�����������Σ�Ҳ���Ը�д�������鷳���
%2016.7.5 ��С�� BUPT
%------------------------------------------------------------%
%��ʼ������
%30GHz 10mm���� λ��-10mm ������1.5
%------------------------------------------------------------%
c=3e8;f=3e10;
lambda=c/f;
w0x=0.01;w0y=w0x;
z0x=-0.010;z0y=z0x;
n1=1;n2=1.5;
% R1=60;R2=R1;
% C11=1/R1;C22=1/R2;
k1=2*pi/lambda;k2=n2/n1*k1;
Ex=1;Ey=0;
theta=(-90:1:90)*pi/180;
distance=50*lambda;


%����ȡ����͸����һ��8.197,0,2cm
%����ϵΪ�˷���ʹ��͸��x�����㣬�����ھ���ԭ����ԭ��12cm��͸����8cm��
%��ʵ�ϣ�������ѡȡ��Ӱ�����ļ���
point=[0.08197,0,0.02];
[n,kk1,kk2]=cal_curvature(point);
% ����ת������������ϵ
% nn1=n(1);nn2=n(2);nn3=n(3);
% n=[nn3,nn2,nn1];

C=[abs(1/kk1),0;0,abs(1/kk2)];

%��Դλ��-12cm����Ϊԭ���ƶ���
feed=[-0.12,0,0];
incident=point-feed;
%��һ�������
thetai=acos(dot(incident,n)/(norm(incident)*norm(n)));
thetat=asin(n1*sin(thetai)/n2);

%���䲨�Ĳ���
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
etai=0;%������㴦���������λΪ0
E0i=sqrt(w0x*w0y/wx/wy)*exp(-1i*k1*z0x+1i*etai);%Ҳ�������Eoi

C11=C(1);C22=C(4);
Qt11=sec(thetat)/k2*(k1*Qi11*cos(thetai)+C11*(k2*cos(thetat)-k1*cos(thetai)));
Qt22=Qi22+C22*(cos(thetat)-k1/k2*cos(thetai));%��������˹�����ļ���

%�������g�ǵļ���򵥣�����ʵ�ʺ���
%���������λ��
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

%������ϵ��
TTE1=2*cos(thetai)*sin(thetat)/sin(thetat+thetai);
TTM1=2*cos(thetai)*sin(thetat)/sin(thetai+thetat)/cos(thetai-thetat);

At=sqrt(w0x*w0y/wx/wy)*exp(1i*etat);
Et0x=Ex*TTE1*E0i/At;Et0y=Ey*TTM1*E0i/At;

%������䲨����ˮƽ�ߵļнǣ����ڼ���
theta1=atan2(incident(3),incident(1));
theta2=theta1-(thetai-thetat);
%��������e-6���ϣ�����Ϊ��ƽ�в���
%������Ƶĳ���
s1=[1,0,tan(theta2)];
s1=s1/norm(s1);
d=0.02;%houdu 
ss1=d-(incident(1)-0.20);
ss2=ss1*tan(theta2);%�����Ĵ�ֱ���룬����Ϊ0
out_point=[0.10,0,incident(3)+ss2];
z0x=-(z0x+ss1);z0y=-(z0y+ss1);
%��������ǣ���Ϊԭ��������һЩ��z0x�Ǹ�ֵ���Ҳ࣬��absС��ss1����Ϊ���������ڣ�Ϊ��


%����ڶ�������ĳ�ǿ
X=0;Y=0;Z=sqrt(ss1^2+ss2^2);
qx=1/Qt11;qy=1/Qt22;
temp=sqrt(w0x*w0y/wx/wy)*exp(-1i*k2*(Z+0.5*(X^2/qx+Y^2/qy)+1i*etat));
H1=[Et0x*temp,Et0y*temp,0];%�ڵڶ������������䲨��ǿ


%�ڶ�������
zrx=w0x^2*pi*n1/lambda;
Rx=(z-z0x)*(1+(zrx/(z-z0x))^2);
wx=w0x*sqrt(1+((z-z0x)/zrx)^2);
Qi11=1/Rx-1i*lambda/pi/n1/wx^2;
zry=w0y^2*pi*n1/lambda;
Ry=(z-z0y)*(1+(zry/(z-z0y))^2);
wy=w0y*sqrt(1+((z-z0y)/zry)^2);
Qi22=1/Ry-1i*lambda/pi/n1/wy^2;
Q=[Qi11,0;0,Qi22];

thetai=atan2(ss2,ss1);
thetat=asin(n2*sin(thetai)/n1);

C11=0;C22=0;%�ڶ������� ƽ�� ��ֱ����
Qt11=sec(thetat)/k1*(k2*Qi11*cos(thetai)+C11*(k1*cos(thetat)-k2*cos(thetai)));
Qt22=Qi22+C22*(cos(thetat)-k2/k1*cos(thetai));
Qt=[Qt11,0;0,Qt22];

%������ڶ������䲨�Ĳ���
gx=-real(Qt11)/imag(Qt11);
Rx=1/real(Qt11);
wx=sqrt(-lambda/pi/n2/imag(Qt11));
w0x=sqrt(wx^2/(1+gx^2));
z0x=-Rx*gx^2/(1+gx^2);
zrx=w0x^2*pi*n1/lambda;
gy=-real(Qt22)/imag(Qt22);
Ry=1/real(Qt22);
wy=sqrt(-lambda/pi/n2/imag(Qt22));
w0y=sqrt(wy^2/(1+gy^2));
z0y=-Ry*gy^2/(1+gy^2);
zry=w0y^2*pi*n1/lambda;
% etat=0.5*atan2(z-z0x,zrx)+0.5*atan2(z-z0y,zry);

%�������
TTE2=1;
TTM2=1;

At=sqrt(w0x*w0y/wx/wy)*exp(1i*etat);
Et0x=TTM2*H1(1)/At;Et0y=TTE2*H1(2)/At;
H=zeros(length(theta),3);
x=(-200*lambda:5*lambda:200*lambda);


 for kk=1:length(theta)
    Yt=0;
    Zt=distance*cos(theta(kk));
    Xt=distance*sin(theta(kk));
%     point=[X,Y,Z];
%     pp=point*C_refract2principle;%main coordinate
%     Xt=pp(1);Yt=pp(2);Zt=pp(3);

    Rx=(Zt-z0x)*(1+(zrx/(Zt-z0x))^2);
    wx=w0x*sqrt(1+((Zt-z0x)/zry)^2);
    qx=1/(1/Rx-1i*lambda/pi/n2/wx^2);
    Ry=(Zt-z0y)*(1+(zry/(Zt-z0y))^2);
    wy=w0y*sqrt(1+((Zt-z0y)/zry)^2);
    qy=1/(1/Ry-1i*lambda/pi/n2/wy^2);
    temp=sqrt(w0x*w0y/wx/wy)*exp(-1i*k2*(Zt+0.5*(Xt^2/qx+Yt^2/qy)+1i*etat));
    H(kk,:)=[Et0x*temp,Et0y*temp,0];
 end


plot(theta/pi*180,abs(H(:,1))/max(abs(H(:,1))),'-b')
title('X-direction')
xlabel('angle')
ylabel('amplitude')
grid on

