%% ѹ����֪�ع��㷨����
clear; clc;
m = 256;%�۲�ֵ����
n = 1024;%�ź�x�ĳ���
range = [2,5,6,10];                                     %Ŀ�����,��λ��m
% range = 0:500:9000;
signal_SNR = 20;                            %�ź������,��λ��dB

[x,s] = genFMCW(n,range,signal_SNR);
x = x(:);
Psi = eye(n);%x������ϡ��ģ�����ϡ�����Ϊ��λ��theta=Psi*x
Phi = randn(m,n)/sqrt(m);%��������Ϊ��˹����
A = Phi * Psi;%���о��� 
y = Phi * x;%�õ��۲�����y

%% �ָ��ع��ź�x
% theta = CS_SAMP( y,A,5);
S = 4;
 [y_rows,y_columns] = size(y);
    if y_rows<y_columns
        y = y';%y should be a column vector
    end
    [m,n] = size(A);%���о���AΪM*N����
    theta = zeros(n,1);%�����洢�ָ���theta(������)
    Pos_theta = [];%�������������д洢A��ѡ��������
    Xmax = [];
    r_n = y;%��ʼ���в�(residual)Ϊy
    L = S;%��ʼ������(Size of the finalist in the first stage)
    Stage = 1;%��ʼ��Stage
    IterMax = 10;
    for ii=1:IterMax%������M��
        %(1)Preliminary Test
        product = A'*r_n;%���о���A������в���ڻ�
        [val,pos]=sort(abs(product),'descend');%��������
        Sk = pos(1:L);%ѡ������L��
        Xmaxi = x(Sk);
        Xmax = union(Xmax,Xmaxi);
        %(2)Make Candidate List
        Ck = union(Pos_theta,Sk);
        %(3)Final Test
        if length(Ck)<=m
            At = A(:,Ck);%��A���⼸����ɾ���At
        else
            theta_ls=0;
            break;
        end
        %y=At*theta��������theta����С���˽�(Least Square)
        theta_ls = (At'*At)^(-1)*At'*y;%��С���˽�
        [val,pos]=sort(abs(theta_ls),'descend');%��������
        F = Ck(pos(1:L));
        %(4)Compute Residue
        %A(:,F)*theta_ls��y��A(:,F)�пռ��ϵ�����ͶӰ
        theta_ls1 = (A(:,F)'*A(:,F))^(-1)*A(:,F)'*y;
        r_new = y - A(:,F)*theta_ls1;%���²в�r
        norm_r = norm(r_n);
        norm_r_new = norm(r_new);
        if norm(r_new)<1e-6%halting condition true 
            Pos_theta = F;
            %r_n = r_new;%����r_n�Ա�������µ�r_n
            break;%quit the iteration
        elseif norm(r_new)>=norm(r_n)%stage switching 
            Stage = Stage + 1;%Update the stage index 
            L = Stage*S;%Update the size of finalist
            if ii == IterMax%���һ��ѭ��
                Pos_theta = F;%����Pos_theta����theta_lsƥ�䣬��ֹ����
            end
            %ii = ii - 1;%��������������
        else
            Pos_theta = F;%Update the finalist Fk
            r_n = r_new;%Update the residue
        end
    end
    theta(Pos_theta)=theta_ls1;%�ָ�����theta

%GLRT
phat = mle(real(x));     
% td = sqrt(-2*phat^2*log(0.3));
%% ��ͼ
figure;
% plot(x_r,'k.-');%���x�Ļָ��ź�
% hold on;
plot(real(x),'r');%���ԭ�ź�x
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
% fprintf('\n�ָ��в');
% norm(x_r-x)%�ָ��в�