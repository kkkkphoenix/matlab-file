% input variables
%freq: frequency in Hz
%w0 beam waist for individual sub_GB
%theta_max, in radian
%delta_theta: angle spacing in radian 0.02;
% feed---------
%decay_deg : in degrees
%decay_value: in dB
%linear_x: =1, linear x polarizaion, =others: y polarization  
%-------

% return variables
%Cm: GB coeffs
%directions: direction of GB (theta, phi)
%GBnumber: GB numbers

function [Cm,directions,GBnumber]=generate_GB_coeffs(freq,w0,theta_max,delta_theta,decay_deg,decay_value,linear_x)
%test---------------------------------------%
% freq=14e10;lamda=3e8/freq;
% w0=5*lamda;
% theta_max=10*pi/180; %maxmimum expansion angle
% delta_theta=0.018;
% 
% decay_deg=15;decay_value=-16;linear_x=1;
%----------------------------------------%



lambda=3e8/freq;
k0=2*pi/lambda;
bb=k0*w0^2/2;

delta_theta_y=delta_theta;
delta_theta_x=delta_theta*sin(pi/3);
% (5.68)
P=floor(theta_max/delta_theta_x);
p=-P:1:P;

GBnumber=0;
pq_total=zeros(5000,2);
for ii=1:length(p)
	theta_x_p=p(ii)*delta_theta_x;
    %（5.70）
    Qp=floor( sqrt(theta_max^2-theta_x_p^2)/delta_theta_y )+1;
    % 判断p的奇偶数
    temp=mod(p(ii),2);
    if(temp==0)
        q=-Qp:1:Qp;
        % 从左到右，从下往上取值
        for jj=1:1:length(q)
            theta_y_q=q(jj)*delta_theta_y;  % 弧度单位
%             if( sqrt(theta_x_p^2+theta_y_q^2) <= theta_max )
                GBnumber=GBnumber+1;
                pq_total(GBnumber,1)=theta_x_p;     % 弧度单位
                pq_total(GBnumber,2)=theta_y_q;     % 弧度单位
%             end
        end
    elseif (temp==1)
        q=(-Qp+1):1:Qp;
        % 从左到右，从下往上取值
        for jj=1:1:length(q)
            theta_y_q=(q(jj)-0.5)*delta_theta_y;
%             if( sqrt(theta_x_p^2+theta_y_q^2) <= theta_max )
                GBnumber=GBnumber+1;
                pq_total(GBnumber,1)=theta_x_p;     % 坐标x位置，弧度单位
                pq_total(GBnumber,2)=theta_y_q;     % 坐标y位置，弧度单位
%             end
        end
   end
end

fprintf('Espanded GB is %d\n',GBnumber);
fprintf('Determine GB direction...\n');
directions=zeros(GBnumber,2);
for ii=1:1:GBnumber
    directions(ii,1)=sqrt(pq_total(ii,1)^2+pq_total(ii,2)^2);
	if(directions(ii,1)~=0)
        cosphi=pq_total(ii,1)/directions(ii,1);
        sinphi=pq_total(ii,2)/directions(ii,1);
    elseif (directions(ii,1)==0)
        cosphi=1;
        sinphi=0;
    end
    directions(ii,2)=atan2(sinphi,cosphi);
end
zi_hat=[sin(directions(:,1)).*cos(directions(:,2)),sin(directions(:,1)).*sin(directions(:,2)),cos(directions(:,1))];
% 开始求M*M的矩阵
fprintf('sovling Matrix...\n')
tic
[E_theta,E_phi]=GB_farfield( directions(:,1)*180/pi,decay_deg,decay_value,linear_x);
EV_L=(abs(E_theta)).^2+(abs(E_phi)).^2;
Z11=E_theta* (conj(E_theta)).';
Z22=E_phi* (conj(E_phi)).';
zi_L_1=zi_hat(:,1)* (zi_hat(:,1)).';
zi_L_2=zi_hat(:,2)* (zi_hat(:,2)).';
zi_L_3=zi_hat(:,3)* (zi_hat(:,3)).';
dot_zi_L_M=zi_L_1+zi_L_2+zi_L_3;
norm_zi_L_M=norm(zi_hat)*(norm(zi_hat))';
costheta=dot_zi_L_M./norm_zi_L_M;
%a×b =[a2b3 ? a3b2, a3b1 ? a1b3, a1b2 ? a2b1]
x1=zi_hat(:,2)*(zi_hat(:,3)).'-zi_hat(:,3)*(zi_hat(:,2)).';
x2=zi_hat(:,3)*(zi_hat(:,1)).'-zi_hat(:,1)*(zi_hat(:,3)).';
x3=zi_hat(:,1)*(zi_hat(:,2)).'-zi_hat(:,2)*(zi_hat(:,1)).';
cross_zi_L_M=sqrt(x1.^2+x2.^2+x3.^2);
sintheta=(cross_zi_L_M)./norm_zi_L_M;
E_Matrix=( Z11+Z22 ).*exp(-0.5*k0*bb*(sintheta./costheta).^2);
toc
fprintf('solving condition of matrix....\n')
tic
vv=rcond(E_Matrix);
toc
fprintf('condition of matrix is %e\n',vv)
opts.UT=false;
opts.TRANSA=false;
fprintf('solving matrix....\n')
tic
Im=E_Matrix\EV_L;
toc
% Cm是系数
Cm=Im/(1i*bb);
fprintf('GB expansion finished\n')

% name=strcat('coeffs','_',num2str(freq/1e9),'_',num2str(w0/lambda),'.txt');
% fprintf('writing coeffs to %s\n',name);
% T=[freq bb GBnumber delta_theta];
% dlmwrite(name,T);
% for m=1:GBnumber
%     T=[real(Cm(m)) imag(Cm(m)) directions(m,1) directions(m,2)];
%     dlmwrite(name,T,'-append','delimiter',' ');
% end