% �μ�GBԶ����ʽ���Ƶ�����
%theta in degrees
function [E_theta,E_phi]=GB_farfield(theta,decay_deg,decay_value,linear_x)
    % ����ת��Ϊ�Ƕ�
    % theta=theta*180/pi;
    % Decay_deg=15;
    % Decay_value=-16;
    if(1==linear_x)
        E_theta=10.^(decay_value/20*(theta/decay_deg).^2);
        E_phi=zeros(size(theta));
    else
        E_theta=zeros(size(theta));
        E_phi=10.^(decay_value/20*(theta/decay_deg).^2);
    end
end





