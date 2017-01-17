%test
freq=14e10;lamda=3e8/freq;
k0=2*pi/lamda;
w0=5*lamda;
theta_max=10*pi/180; %maxmimum expansion angle
delta_theta=2*pi/180;

decay_deg=15;decay_value=-16;linear_x=1;

[Cm,directions,GBnumber]=generate_GB_coeffs(freq,w0,theta_max,delta_theta,decay_deg,decay_value,linear_x);

%-----------------------------------------%
%求出恢复场值%
theta=directions(:,1);
phi=directions(:,2);
R=0.001;

E_recover=zeros(length(theta),2);

for i=length(theta)
    %波束坐标系
    xx=R*sin(theta(i))*cos(phi(i));
    yy=R*sin(theta(i))*sin(phi(i));
    zz=R*cos(theta(i));
    point=[xx,yy,zz];
      for j=1:GBnumber
          em=10^(decay_value/20*(theta(i)/decay_deg*pi/180));
          
          %坐标变换
          zi_hat=[sin(theta(i))*cos(phi(i)),sin(theta(i))*sin(phi(i)),cos(theta(i))];
          yi_hat=cross(zi_hat,x_hat);
          yi_hat=yi_hat/norm(yi_hat);
          xi_hat=cross(yi_hat,zi_hat);
          xi_hat=xi_hat/norm(xi_hat);
          
          x_hat=[1,0,0];
          y_hat=[0,1,0];
          z_hat=[0,0,1];
          
          c11=dot(xi_hat,x_hat);c12=dot(yi_hat,x_hat);c13=dot(zi_hat,x_hat);
          c21=dot(xi_hat,y_hat);c22=dot(yi_hat,y_hat);c23=dot(zi_hat,y_hat);
          c31=dot(xi_hat,z_hat);c32=dot(yi_hat,z_hat);c33=dot(zi_hat,z_hat);
          %坐标变换矩阵
          C=[c11,c12,c12;c21,c22,c23;c31,c32,c33];
          %在新坐标中的点
          new_xx=dot(C(1,:),point);
          new_yy=dot(C(2,:),point);
          new_zz=dot(C(3,:),point);
          
          %(5.57)高斯基模函数
          E_i_1=emEy_theta*(1i*b/(new_zz+1i*b))*exp(-1i*k0*(new_zz+1/2*(new_xx^2+new_yy^2)/(new_zz+1i*b)));
          E_i_2=emEy_phi*(1i*b/(new_zz+1i*b))*exp(-1i*k0*(new_zz+1/2*(new_xx^2+new_yy^2)/(new_zz+1i*b)));
          
          E_recover(i,1)=E_recover(i,1)+Cm(j)*E_i_1;
          E_recover(i,2)=E_recover(i,2)+Cm(j)*E_i_2;   
      end
      E_recover_total(i)=sqrt(abs(E_recover(i,1))^2+abs(E_recover(i,2))^2);
end
      