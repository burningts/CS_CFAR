%% ѹ����֪�ع��㷨����
clear; clc;
m = 256;%�۲�ֵ����
n = 1024;%�ź�x�ĳ���
range = [2,5,6,10,11,12,15];                         %Ŀ�����,��λ��m
% range = 0:500:9000;
signal_SNR = 20;                            %�ź������,��λ��dB

[x,s] = genFMCW(2*n,range,signal_SNR);
x = x(:);
Psi = eye(n);%x������ϡ��ģ�����ϡ�����Ϊ��λ��theta=Psi*x
Phi = randn(m,n)/sqrt(m);%��������Ϊ��˹����
A = Phi * Psi;%���о��� 
y = Phi * x;%�õ��۲�����y

%% �ָ��ع��ź�x
phat = mle(real(x));
sigma = phat(1);
TargetPfa = 0.1;
[ Pos_theta,theta_ls1,Ta,Pfa ] = CS_SAMP(x,y,A,4,20,TargetPfa,sigma );

%% ��ͼ
figure;
% plot(x_r,'k.-');%���x�Ļָ��ź�
% hold on;
plot(x,'r');%���ԭ�ź�x
% figure;
% h =  histogram(x,'Normalization','pdf');
% hold on
% y = 0:1:700;
% mu = 5;
% sigma = 27;
% 
% % f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% f = y.*exp(-y.^2/(2*sigma^2))/sigma^2;
% plot(y,f,'LineWidth',1.5);
% disp(h.NumBins);
% hold off;
% legend('Recovery','Original')
% fprintf('\n�ָ��в');
% norm(x_r-x)%�ָ��в�