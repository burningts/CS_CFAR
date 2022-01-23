%% 压缩感知重构算法测试
clear; clc;
m = 256;%观测值个数
n = 1024;%信号x的长度
range = [2,5,6,10];                                     %目标距离,单位：m
% range = 0:500:9000;
signal_SNR = 20;                            %信号信噪比,单位：dB

[x,s] = genFMCW(n,range,signal_SNR);
x = x(:);
Psi = eye(n);%x本身是稀疏的，定义稀疏矩阵为单位阵theta=Psi*x
Phi = randn(m,n)/sqrt(m);%测量矩阵为高斯矩阵
A = Phi * Psi;%传感矩阵 
y = Phi * x;%得到观测向量y

%% 恢复重构信号x
% theta = CS_SAMP( y,A,5);
S = 4;
 [y_rows,y_columns] = size(y);
    if y_rows<y_columns
        y = y';%y should be a column vector
    end
    [m,n] = size(A);%传感矩阵A为M*N矩阵
    theta = zeros(n,1);%用来存储恢复的theta(列向量)
    Pos_theta = [];%用来迭代过程中存储A被选择的列序号
    Xmax = [];
    r_n = y;%初始化残差(residual)为y
    L = S;%初始化步长(Size of the finalist in the first stage)
    Stage = 1;%初始化Stage
    IterMax = 10;
    for ii=1:IterMax%最多迭代M次
        %(1)Preliminary Test
        product = A'*r_n;%传感矩阵A各列与残差的内积
        [val,pos]=sort(abs(product),'descend');%降序排列
        Sk = pos(1:L);%选出最大的L个
        Xmaxi = x(Sk);
        Xmax = union(Xmax,Xmaxi);
        %(2)Make Candidate List
        Ck = union(Pos_theta,Sk);
        %(3)Final Test
        if length(Ck)<=m
            At = A(:,Ck);%将A的这几列组成矩阵At
        else
            theta_ls=0;
            break;
        end
        %y=At*theta，以下求theta的最小二乘解(Least Square)
        theta_ls = (At'*At)^(-1)*At'*y;%最小二乘解
        [val,pos]=sort(abs(theta_ls),'descend');%降序排列
        F = Ck(pos(1:L));
        %(4)Compute Residue
        %A(:,F)*theta_ls是y在A(:,F)列空间上的正交投影
        theta_ls1 = (A(:,F)'*A(:,F))^(-1)*A(:,F)'*y;
        r_new = y - A(:,F)*theta_ls1;%更新残差r
        norm_r = norm(r_n);
        norm_r_new = norm(r_new);
        if norm(r_new)<1e-6%halting condition true 
            Pos_theta = F;
            %r_n = r_new;%更新r_n以便输出最新的r_n
            break;%quit the iteration
        elseif norm(r_new)>=norm(r_n)%stage switching 
            Stage = Stage + 1;%Update the stage index 
            L = Stage*S;%Update the size of finalist
            if ii == IterMax%最后一次循环
                Pos_theta = F;%更新Pos_theta以与theta_ls匹配，防止报错
            end
            %ii = ii - 1;%迭代次数不更新
        else
            Pos_theta = F;%Update the finalist Fk
            r_n = r_new;%Update the residue
        end
    end
    theta(Pos_theta)=theta_ls1;%恢复出的theta

%GLRT
phat = mle(real(x));     
% td = sqrt(-2*phat^2*log(0.3));
%% 绘图
figure;
% plot(x_r,'k.-');%绘出x的恢复信号
% hold on;
plot(real(x),'r');%绘出原信号x
% figure;
% h =  histogram(real(x),'Normalization','pdf');
% hold on
% y = 0:1:700;
% mu = 5;
% sigma = 35;
% % f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% f = y.*exp(-y.^2/(2*sigma^2))/sigma^2;
% plot(y,f,'LineWidth',1.5)
% disp(h.NumBins);
% hold off;
% legend('Recovery','Original')
% fprintf('\n恢复残差：');
% norm(x_r-x)%恢复残差