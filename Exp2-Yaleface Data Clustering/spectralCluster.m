function [C,bre] = spectralCluster(X,K,Sigma)
format long; warning off
% 各维度数据归一化
X=X./max(abs(X),[],1:ndims(X));

% 计算相似矩阵即||x_i-x_j||^2
S=pdist2(X,X).^2;

% 高斯核函数RBF计算邻接矩阵
W=exp(-S./(2.*Sigma.*Sigma));

% 计算度矩阵
D=diag(sum(W,2));

% 计算拉普拉斯矩阵
L=D-W;
L_norm=D^(-1/2)*L*D^(-1/2);


% 求特征向量 V
% 'SM':绝对值最小特征值排列
[V, bre,~]=eigs0(L_norm,K,'SM');

if bre == 1
    C = X;
    return
end

% 按行标准化
V_norm=V./vecnorm(V.')';

% 对特征向量合并矩阵的行向量求k-means
C = kmedoids(V_norm,K);
end

