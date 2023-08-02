    function [idx] = main(X0,Gr,E)

    epsi = 1e-5;   
    [n,k] = size(X0);    
    G      = Gr.edge;    
    [~,r]   = size(G);  
    n2      = max(G,[],'all'); 
    if n~=n2
        fprintf('The vertex numbers in X0 and G are not equal\n');
        return;
    end
    Lambda = 0.1;      
    delta = 1;          
    alpha  = 0.8;          
    iterNum = 1000;       
     fbest  = 1.0000e+10;  
     l   = 0;              
     L   = 1;            

    deg_ik = Gr.degree.^(1/r);
    [U] = SparseU(G,X0,Gr.weights,deg_ik);


    %% Initialization.
    
       
    
    X_E = X_G(X0,E,n,k);
    [Nablaf0,~,~] =  Grad(X0,X_E,E,Lambda,r,U); 
      X1 = X0- Nablaf0;
     X1 = Proj(X1,n,k);
     X_E = X_G(X1,E,n,k);
     [Nablaf01] =  Grad(X1,X_E,E,Lambda,r,U);     
     fr = Objf(U,X1,r,Lambda,E,X_E);
      
     iter = 0;


    %% 
    while iter < iterNum
         if norm(Nablaf01,'fro')<epsi
             break
         end
         
         s = X1-X0;
         y = -Nablaf01+Nablaf0;
       if mod(iter,2)==1
          t0 = trace(s'*s)/abs(trace(s'*y));
       else
          t0 = abs(trace(s'*y))/trace(y'*y);
       end
  
         
    % Armijo search :  
      for m=0:20

         
          Xk = -t0*alpha^m*Nablaf01+X1;
          Xk = Proj(Xk,n,k);
          X_E = X_G(Xk,E,n,k);      
          f1 = Objf(U,Xk,r,Lambda,E,X_E);

       if f1 <= fr-delta*alpha^m*t0*norm(Nablaf01)^2
          break;   
       end     
       
       end
    % Armijo search finished
    
          X0 = X1;
          X1 = Xk;
          Nablaf0 = Nablaf01;
         [Nablaf01,~, ~] = Grad(X1,X_E,E,Lambda,r,U);   %  
          iter = iter+1;

    %  fr updating :
           if f1<fbest
              fbest = f1; fc = f1; l=0;
          else
              l = l+1;
              fc = max(fc,f1); 
              if l == L
                 l = 0; fr = fc; fc = f1;
              end
           end
% fr updating 
    end


      idx = spectralCluster(X1,k,1);      


    end 
    %%
    function [ProjX] = Proj(X,n,k)

    [U,~,V] = svd(X);
    ProjX = U*eye(n,k)*V';
    end

    %%
    function C = spectralCluster(X,K,Sigma)
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
    [V, ~]=eigs(L_norm,K,'SM');

    % 按行标准化
    V_norm=V./vecnorm(V.')';

    % 对特征向量合并矩阵的行向量求k-means
    C = kmedoids(V_norm,K);
    end








