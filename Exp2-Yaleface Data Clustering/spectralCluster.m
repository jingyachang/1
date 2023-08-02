function [C,bre] = spectralCluster(X,K,Sigma)
format long; warning off
% ��ά�����ݹ�һ��
X=X./max(abs(X),[],1:ndims(X));

% �������ƾ���||x_i-x_j||^2
S=pdist2(X,X).^2;

% ��˹�˺���RBF�����ڽӾ���
W=exp(-S./(2.*Sigma.*Sigma));

% ����Ⱦ���
D=diag(sum(W,2));

% ����������˹����
L=D-W;
L_norm=D^(-1/2)*L*D^(-1/2);


% ���������� V
% 'SM':����ֵ��С����ֵ����
[V, bre,~]=eigs0(L_norm,K,'SM');

if bre == 1
    C = X;
    return
end

% ���б�׼��
V_norm=V./vecnorm(V.')';

% �����������ϲ��������������k-means
C = kmedoids(V_norm,K);
end

