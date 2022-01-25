function [ Pos_theta,theta_ls1,Ta,Pfa ] = CS_SAMP(x,y,A,S,count,TargetPfa,sigma )
%CS_SAMP Summary of this function goes here
%Version: 1.0 written by jbb0523 @2015-05-08
%   Detailed explanation goes here
%   y = Phi * x
%   x = Psi * theta
%	y = Phi*Psi * theta
%   令 A = Phi*Psi, 则y=A*theta
%   现在已知y和A，求theta
%   Reference:Thong T.Do，Lu Gan，Nam Nguyen，Trac D.Tran．Sparsity adaptive
%   matching pursuit algorithm for practical compressed sensing[C]．Asilomar
%   Conference on Signals，Systems，and Computers，Pacific Grove，California，
%   2008，10：581-587.
%   Available at:
%   http://dsp.rice.edu/sites/dsp.rice.edu/files/cs/asilomar08_final.pdf
    [y_rows,y_columns] = size(y);
    if y_rows<y_columns
        y = y';%y should be a column vector
    end
    [m,n] = size(A);%传感矩阵A为M*N矩阵
    theta = zeros(n,1);%用来存储恢复的theta(列向量)
    Pos_theta = [];%用来迭代过程中存储A被选择的列序号
    Xmax = [];
    Ta = [];
    r_n = y;%初始化残差(residual)为y
    L = S;%初始化步长(Size of the finalist in the first stage)
    Stage = 1;%初始化Stage
    IterMax = count;
    for ii=1:IterMax%最多迭代M次
        %(1)Preliminary Test
        product = A'*r_n;%传感矩阵A各列与残差的内积
        [val,pos]=sort(abs(product),'descend');%降序排列s
        Sk = pos(1:L);%选出最大的L个
        Xmaxi = x(Sk);
        Xmax = union(Xmax,Xmaxi);
        Ta = mean(Xmaxi);
    %   Ta = union(Ta,Tai);        
        %(2)Make Candidate List
        Ck = union(Pos_theta,Sk);
        %(3)Final Test
        if length(Ck)<=256
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
        Pfa = exp(-Ta.^2/(2*sigma^2));
        if Pfa >= TargetPfa   %halting condition true
            break;      %quit the iteration    
        % if norm(r_new)<1e-6 
        %     Pos_theta = F;
        %     %r_n = r_new;%更新r_n以便输出最新的r_n
        %     break;
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
end