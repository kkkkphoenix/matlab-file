clear all
close all;
clc
freq=30*1e9;
lamda=3e8/freq;	
k0=2*pi/lamda;	
%-1: Zogbi
%1: DGBA
sign=1;
% R1=-100*lamda;
% R2=-100*lamda;
R1= -2.311307969496022;
R2= -3.086835306513838;
% bi1=100*lamda;
% bi2=100*lamda;
bi1=0.245;
bi2=0.245;
wi1=sqrt(2*bi1/k0);
wi2=sqrt(2*bi2/k0);
wi1/lamda
wi2/lamda



%field coordinate
field_origin=[1.154700000000000 0 0.333333022500000];
field_x_hat=[1,0,0];
field_y_hat=[0,1,0];
field_z_hat=[0,0,1];
x_start=-0.5;
x_end=0.5;
np_x=100;
y_start=-0.5;
y_end=0.5;
np_y=100;

% rhoi1=25*lamda;
% rhoi2=35*lamda;
% rc1=90*lamda;
% rc2=90*lamda;
% rc1=50*lamda;
% rc2=50*lamda;
rc1= 0.938357776433826;
rc2= 0.938357776433826;
if(-1==sign)
    rhoi1=rc1+1i*bi1;
    rhoi2=rc2+1i*bi2;
else
    rhoi1=rc1-1i*bi1;
    rhoi2=rc2-1i*bi2;
end
zc=pi*wi1^2/lamda;
w=wi1*(1+(rc1/zc)^2)^0.5;
w/lamda

theta_i=10*pi/180;
theta=pi+theta_i;
phi=-10*pi/180;

% x_hat=[1,0,0];
% y_hat=[0,1,0];
% z_hat=[0,0,1];

x_hat=[-0.114547024531772  -0.993417827085319 -0.000000000000001];
y_hat=[0.859615282944689  -0.099118789917685   0.501235504342196];
z_hat=[-0.497936285581639   0.057415035612080   0.865310908972506];

% zi_hat=[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)];
% yi_hat=cross(z_hat,zi_hat);
% yi_hat=yi_hat/norm(yi_hat);
% xi_hat=cross(yi_hat,zi_hat);
% xi_hat=xi_hat/norm(xi_hat);

xi_hat=[-0.500000349687582 0 -0.866025201892125];
yi_hat=[0.122474458587274   0.989949493661167  -0.070710727571947];
zi_hat=[0.857321210110919  -0.141421356237310  -0.494975093003628];
c11=dot(x_hat,xi_hat);c12=dot(y_hat,xi_hat);c13=dot(z_hat,xi_hat);
c21=dot(x_hat,yi_hat);c22=dot(y_hat,yi_hat);c23=dot(z_hat,yi_hat);
c31=dot(x_hat,zi_hat);c32=dot(y_hat,zi_hat);c33=dot(z_hat,zi_hat);
C=[c11,c12,c13;c21,c22,c23;c31,c32,c33];
            
%incident GB: page 106
Qi11=1/(rhoi1);
Qi12=0;
Qi21=0;
Qi22=1/(rhoi2)         ;
Qi=[Qi11,Qi12;Qi21,Qi22];

%for reflection as a GB
%zogbi, 6.102-6.104
Qr11=Qi(1,1)-2/C(3,3)*((C(2,2))^2/R1+((C(2,1))^2)/R2);
Qr12=Qi(1,2)+2/C(3,3)*((C(1,2)*C(2,2))/R1+(C(2,1)*C(1,1))/R2);
Qr22=Qi(2,2)-2/C(3,3)*((C(1,2))^2/R1+((C(1,1))^2)/R2);
%zogbi,6.100-6.101
rhor1=2/((Qr11+Qr22)+sqrt((Qr11+Qr22)^2-4*(Qr11*Qr22-Qr12^2)));
rhor2=2/((Qr11+Qr22)-sqrt((Qr11+Qr22)^2-4*(Qr11*Qr22-Qr12^2)));
%反射坐标系(在主轴坐标系中的表示)
XR_hat=xi_hat-2*dot(z_hat,xi_hat)*z_hat;
XR_hat=XR_hat/norm(XR_hat);
YR_hat=yi_hat-2*dot(z_hat,yi_hat)*z_hat;
YR_hat=YR_hat/norm(YR_hat);
%ZR_hat=cross(XR_hat,YR_hat);
ZR_hat=zi_hat-2*dot(z_hat,zi_hat)*z_hat;
ZR_hat=ZR_hat/norm(ZR_hat);

d11=dot(XR_hat,x_hat);d12=dot(YR_hat,x_hat);d13=dot(ZR_hat,x_hat);
d21=dot(XR_hat,y_hat);d22=dot(YR_hat,y_hat);d23=dot(ZR_hat,y_hat);
d31=dot(XR_hat,z_hat);d32=dot(YR_hat,z_hat);d33=dot(ZR_hat,z_hat);
C_reflect2principal=[d11,d12,d13;d21,d22,d23;d31,d32,d33];

%DGBA
% q0=-1i*pi*wi1^2/lamda;
% q=q0+rc1;
% scale=q0/q*exp(1i*k0*rc1);
% Ki1=rc1+zc^2/rc1;
% Ki2=Ki1;
% ki1_hat=xi_hat;
% ki2_hat=yi_hat;
% Kc1=R1;Kc2=R2;
% %[Cr1,Cr2,cr1,cr2] = reflect(z_hat,zi_hat,ki1_hat,ki2_hat,1/Ki1,1/Ki2,kc1_hat,1/Kc1,1/Kc2);
% Kr1=1/Cr1;
% Kr2=1/Cr2;
% XXR_hat=cr1;
% YYR_hat=cr2;
% cos_angle=dot(-zi_hat,z_hat);
% sr_hat=-zi_hat+2*cos_angle*z_hat;
% sr_hat=sr_hat/absv(sr_hat);
% ZZR_hat=sr_hat;
% 
% d11=dot(XXR_hat,x_hat);d12=dot(YYR_hat,x_hat);d13=dot(ZZR_hat,x_hat);
% d21=dot(XXR_hat,y_hat);d22=dot(YYR_hat,y_hat);d23=dot(ZZR_hat,y_hat);
% d31=dot(XXR_hat,z_hat);d32=dot(YYR_hat,z_hat);d33=dot(ZZR_hat,z_hat);
% C_reflect2principal_DGBA=[d11,d12,d13;d21,d22,d23;d31,d32,d33];
% 
% w=wi1*(1+(lamda*rc1/(pi*wi1^2))^2)^0.5;
% Q1R=1/(1/Kr1+1i*lamda/(pi*w^2));
% Q2R=1/(1/Kr2+1i*lamda/(pi*w^2));

%6.15-6.17
% H1=10; %x direction
% H2=50; %y direction
H1=0.001750640768269 - 0.097303996727819i;
H2=1.336226161036131e-019 -7.427003206925225e-018i;
% H2=0.014785392670062-0.001848802452484i;
W1=H1*c11+H2*c21;
W2=H1*c12+H2*c22;
W3=H1*c13+H2*c23;
W=[W1,W2,W3];
delta=0.4*lamda;

point_inc=[1.150884105166089  -0.132703829379098   0.335536132463956];
bound=15*lamda;
xx=(-bound:delta:bound)+point_inc(1);
yy=(-bound:delta:bound)+point_inc(2);
theta=(-90:1:90)*pi/180;
XX=(-150*lamda:1*lamda:150*lamda);
distance=400*lamda;
H=zeros(length(theta),3);
H_aproxmation_Zogbi=H;
H_aproxmation_Chou=H;
H_aproxmation_GB=H;
H_aproxmation_DGBA=H;


% field_origin=[x_center,y_center,z];
% field_x_hat=[1,0,0];
% field_y_hat=[0,1,0];
% field_z_hat=[0,0,1];
% x_start=-0.5;
% x_end=0.5;
% np_x=100;
% ppp=linspace(x_start,x_end,np_x);

% distance=1.333333022500000;

for kk=1:length(theta)
    X=0;
%     if(theta(kk)==0)
%         X=X;
%     end
    fprintf('theta=%d\n',theta(kk)*180/pi);
%     X=XX(kk);
%     Z=distance;
    Z=distance*cos(theta(kk));
    Y=distance*sin(theta(kk));
    point=[X,Y,Z];
    %r=sqrt(X^2+Y^2+Z^2);
    %转到主轴坐标系下
    pp=C_reflect2principal*point';
%     temp1=cal_reflection_Chou(Qi,C,H1,H2,pp,R1,R2,x_hat,y_hat,k0);
%     temp1=(C_reflect2principal)\temp1';
%     temp1=temp1';
%     H_aproxmation_Chou(kk,:)=temp1;
    temp2=cal_reflection_Zogbi(Qi,C,W,pp,R1,R2,x_hat,y_hat,k0,sign);
    temp2=(C_reflect2principal)\temp2';
    temp2=inv(C_reflect2principal)*temp2';
    temp2=temp2';
    H_aproxmation_Zogbi(kk,:)=temp2;
%     %投影到反射坐标系
%     XR=sum(XR_hat.*point);
%     YR=sum(YR_hat.*point);
%     ZR=sum(ZR_hat.*point);
    XR=X;YR=Y;ZR=Z;
    %zogbi,6.107
    if(-1==sign)
        temp=sqrt(rhor1*rhor2/((ZR+rhor1)*(ZR+rhor2)));
    else
        temp=sqrt(rhor1*rhor2/((ZR+rhor1)*(ZR+rhor2)));
    end
    if(imag(temp)<0)
        temp=-temp;
    end
    temp1=XR^2*(ZR+rhor1*rhor2*Qr11)+YR^2*(ZR+rhor1*rhor2*Qr22)+2*XR*YR*rhor1*rhor2*Qr12;
    temp2=(ZR+rhor1)*(ZR+rhor2);
    if(-1==sign)
        phase=exp(-1i*k0*ZR)*exp(-1i*k0/2*(temp1/temp2));
    else
        phase=exp(1i*k0*ZR)*exp(1i*k0/2*(temp1/temp2));
    end
    temp=temp*phase;
    H_aproxmation_GB(kk,:)=-[H1*temp,H2*temp,0];
    %计算得到的是反射坐标中的场分布
    %重新投影到主轴坐标系
    H_temp=C_reflect2principal*H_temp';
    H_aproxmation_GB(kk,:)=H_temp';
    
%     %DGBA
%     %投影到反射坐标系
%     XR=sum(XXR_hat.*point);
%     YR=sum(YYR_hat.*point);
%     ZR=sum(ZZR_hat.*point);
%     temp=sqrt(Q1R*Q2R/(Q1R+ZR)/(Q2R+ZR));
%     if(imag(temp<0))
%         temp=-temp;
%     end
%     FGES=scale*temp*exp(1i*k0*(ZR+XR^2/(2*(Q1R+ZR))+YR^2/(2*(Q2R+ZR))));
%     %计算得到的是反射坐标中的场分布
%     %重新投影到主轴坐标系
%     H_temp=[FGES*co,FGES*cx,0];
%     H_temp=C_reflect2principal*H_temp';
%     H_aproxmation_DGBA(kkk,:)=H_temp';
    
    %物理积分
%     temp=zeros(1,3);
%     X=pp(1);Y=pp(2);Z=pp(3);
%     for m=1:length(xx)
%         for n=1:length(yy)
%             zz=-0.5*(xx(m)^2/R1+yy(n)^2/R2);
%             R=sqrt((X-xx(m))^2+(Y-yy(n))^2+(Z-zz)^2);
%             %6.25
% %             if(-1==sign)
% %                 B_temp1=exp(-1i*k0*(xx(m)*c31+yy(n)*c32));
% %                 B_temp2=exp(-1i*k0/2*xx(m)^2*(Qi11*c11^2+Qi22*c21^2+2*Qi12*c11*c21-c33/R1));
% %                 B_temp3=exp(-1i*k0/2*yy(n)^2*(Qi11*c12^2+Qi22*c22^2+2*Qi12*c12*c22-c33/R2));
% %                 B_temp4=exp(-1i*k0/2*2*xx(m)*yy(n)*(Qi11*c11*c12+Qi22*c21*c22+Qi21*(c11*c22+c21*c12)));
% %                 B=B_temp1*B_temp2*B_temp3*B_temp4;
% %             else
% %                 B_temp1=exp(1i*k0*(xx(m)*c31+yy(n)*c32));
% %                 B_temp2=exp(1i*k0/2*xx(m)^2*(Qi11*c11^2+Qi22*c21^2+2*Qi12*c11*c21-c33/R1));
% %                 B_temp3=exp(1i*k0/2*yy(n)^2*(Qi11*c12^2+Qi22*c22^2+2*Qi12*c12*c22-c33/R2));
% %                 B_temp4=exp(1i*k0/2*2*xx(m)*yy(n)*(Qi11*c11*c12+Qi22*c21*c22+Qi21*(c11*c22+c21*c12)));
% %                 B=B_temp1*B_temp2*B_temp3*B_temp4;
% %             end
%             point_i=[xx(m),yy(n),zz];
%             point_i=C*point_i';
%             xi=point_i(1);
%             yi=point_i(2);
%             zi=point_i(3);
%             phi=-1i*(-sign)*k0/2*((xi^2*(zi+rhoi1*rhoi2*Qi11)+yi^2*(zi+rhoi1*rhoi2*Qi22)+2*xi*yi*rhoi1*rhoi2*Qi12)/((zi+rhoi1)*(zi+rhoi2)));
%             %temp1=sqrt(rhoi1/(zi+rhoi1));
%             %temp2=sqrt(rhoi2/(zi+rhoi2));
%             temp1=sqrt(rhoi1*rhoi2/(zi+rhoi1)/(zi+rhoi2));
% %             if(imag(temp1)<0)
% %                 temp1=-temp1;
% %             end
%             B=temp1*exp(-1i*(-sign)*k0*zi)*exp(phi);
%             vx=-W2+yy(n)*W3/R2;
%             vy=W1-xx(m)*W3/R1;
%             vz=xx(m)*W2/R1-yy(n)*W1/R2;
%             W_bar=x_hat*(vz*(Y-yy(n))-vy*(Z-zz))/R+y_hat*(vx*(Z-zz)-vz*(X-xx(m)))/R+z_hat*(vy*(X-xx(m))-vx*(Y-yy(n)))/R;
%             temp=temp-1i*(-sign)*k0/(2*pi)*W_bar*B*exp(-1i*(-sign)*k0*R)/R*delta*delta;
%             %H(kk,:)=H(kk,:)-1i*k0/(2*pi)*W_bar*B*exp(-1i*k0*R)/R*delta*delta;
%         end
%     end
%     temp=(C_reflect2principal)\temp';
%     temp=temp';
%     H(kk,:)=temp;
%  end

% plot(theta*180/pi,abs(H(:,2)),'linewidth',1)
% hold on
% plot(theta*180/pi,abs(H_aproxmation_Zogbi(:,2)),'--r','linewidth',1)
% hold on
% plot(theta*180/pi,abs(H_aproxmation_Chou(:,2)),'*-g','linewidth',1)
% hold on
% plot(theta*180/pi,abs(H_aproxmation_GB(:,2)),'+k','linewidth',1)
% plot(XX,20*log10(abs(H(:,2))/max(abs(H(:,2)))),'linewidth',1)
% hold on
% plot(XX,20*log10(abs(H_aproxmation_Zogbi(:,2))/max(abs(H_aproxmation_Zogbi(:,2)))),'--r','linewidth',1)
% hold on
% plot(XX,20*log10(abs(H_aproxmation_Chou(:,2))/max(abs(H_aproxmation_Chou(:,2)))),'*-g','linewidth',1)
% hold on
% plot(XX,20*log10(abs(H_aproxmation_GB(:,2))/max(abs(H_aproxmation_GB(:,2)))),'+k','linewidth',1)
index=1;
plot(theta*180/pi,(abs(H(:,index))),'linewidth',1)
hold on
plot(theta*180/pi,(abs(H_aproxmation_Zogbi(:,index))),'--r','linewidth',1)
hold on
plot(theta*180/pi,(abs(H_aproxmation_Chou(:,index))),'*-g','linewidth',1)
hold on
plot(theta*180/pi,(abs(H_aproxmation_GB(:,index))),'+k','linewidth',1)
legend('reference','approximation-Zogbi','aprpoximation-Chou','aprpoximation-GB')
title('H compared(amplitude)')
xlabel('theta(degrees)')
ylabel('Magnetic Field')
set(gca,'Fontsize',12)

figure(2)
plot(theta*180/pi,angle((H(:,index)))*180/pi,'linewidth',1)
hold on
plot(theta*180/pi,angle((H_aproxmation_Zogbi(:,index)))*180/pi,'--r','linewidth',1)
hold on
plot(theta*180/pi,angle((H_aproxmation_Chou(:,index)))*180/pi,'*-g','linewidth',1)
hold on
plot(theta*180/pi,angle((H_aproxmation_GB(:,index)))*180/pi,'+k','linewidth',1)
legend('reference','approximation-Zogbi','aprpoximation-Chou','aprpoximation-GB')
title('H compared(phase)')
xlabel('theta(degrees)')
ylabel('phase(degrees)')
set(gca,'Fontsize',12)

figure(3)
plot(theta*180/pi,abs(angle((H_aproxmation_GB(:,index)))*180/pi-angle((H(:,index)))*180/pi),'-b','linewidth',1)
title('H compared(phase),GB-Reference')
figure(4)
plot(theta*180/pi,abs(angle((H_aproxmation_Zogbi(:,index)))*180/pi-angle((H(:,index)))*180/pi),'-b','linewidth',1)
title('H compared(phase),Zogbi-Reference')