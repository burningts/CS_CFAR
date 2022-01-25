function [ TargetPos,TargetVal,Ta,Pfa,ii ] = CS_SAMP(x,y,SenseMartix,StepSize,...
                                                  TargetPfa,Sigma )
%CS_SAMP Summary of this function goes here
%   y = Phi * x
%   x = Psi * theta
%	y = Phi*Psi * theta
%   令 A = Phi*Psi, 则y=A*theta
%   现在已知y和A，求theta
    [y_rows,y_columns] = size(y);
    if y_rows<y_columns
        y = y';%y should be a column vector
    end
    [m,n] = size(SenseMartix);%传感矩阵A为M*N矩阵
%     theta = zeros(n,1);%用来存储恢复的theta(列向量)
    TargetPos = [];%用来迭代过程中存储A被选择的列序号
    Ta = [];
    r_n = y;%初始化残差(residual)为y
    L = StepSize;%初始化步长(Size of the finalist in the first stage)
    Stage = 1;%初始化Stage
    IterMax = m;
    for ii=1:IterMax%最多迭代M次
        %(1)Preliminary Test
        product = SenseMartix'*r_n;%传感矩阵A各列与残差的内积
        [val,pos]=sort(abs(product),'descend');%降序排列s
        Sk = pos(1:L);%选出最大的L个
        Xmax = x(Sk);
        Ta = mean(Xmax);
    %   Ta = union(Ta,Tai);        
        %(2)Make Candidate List
        Ck = union(TargetPos,Sk);
        %(3)Final Test
        if length(Ck)<= m
            At = SenseMartix(:,Ck);%将A的这几列组成矩阵At
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
        theta_ls = (SenseMartix(:,F)'*SenseMartix(:,F))^(-1)*SenseMartix(:,F)'*y;
        r_new = y - SenseMartix(:,F)*theta_ls;%更新残差r
        norm_r = norm(r_n);
        norm_r_new = norm(r_new);
        Pfa = exp(-Ta.^2/(2*Sigma^2));
        if Pfa >= TargetPfa   %halting condition true
            TargetPos = F;
            break;      %quit the iteration    
        elseif norm(r_new)>=norm(r_n)%stage switching 
            Stage = Stage + 1;%Update the stage index 
            L = Stage*StepSize;%Update the size of finalist
            if ii == IterMax%最后一次循环
                TargetPos = F;%更新Pos_theta以与theta_ls匹配，防止报错
            end
            %ii = ii - 1;%迭代次数不更新
        else
            TargetPos = F;%Update the finalist Fk
            r_n = r_new;%Update the residue
        end
    end
   TargetVal = x(TargetPos);%恢复出的theta
end