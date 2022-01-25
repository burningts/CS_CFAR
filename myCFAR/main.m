%% 压缩感知重构算法测试
clear; clc;
m = 256;%观测值个数
n = 1024;%信号x的长度
range = [2,5,6,10,11,12,15];                %目标距离,单位：m
% range = 0:500:9000;
signal_SNR = 20;                            %信号信噪比,单位：dB

[x,s] = genFMCW(2*n,range,signal_SNR);
x = x(:);
Psi = eye(n);%x本身是稀疏的，定义稀疏矩阵为单位阵theta=Psi*x
Phi = randn(m,n)/sqrt(m);%测量矩阵为高斯矩阵
SenseMartix = Phi * Psi;%传感矩阵 
y = Phi * x;%得到观测向量y

%% 恢复重构信号x
phat = mle(real(x));
Sigma = phat(1);
TargetPfa = 0.01;
StepSize = 3;
[ TargetPos,TargetVal,Ta,Pfa,ii] = CS_SAMP(x,y,SenseMartix,StepSize,TargetPfa,Sigma );

%% 绘图
figure;
plot(x,'r');%绘出原信号x
% figure;
% h =  histogram(x,'Normalization','pdf');
% hold on
% y = 0:1:700;
% mu = 5;
% sigma = 27; 
% % f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% f = y.*exp(-y.^2/(2*sigma^2))/sigma^2;
% plot(y,f,'LineWidth',1.5);
% disp(h.NumBins);
% hold off;
% legend('Recovery','Original')
% fprintf('\n恢复残差：');
% norm(x_r-x)%恢复残差