function [ProjX] = Proj(X,n,k)

    [U,~,V] = svd(X);
    ProjX = U*eye(n,k)*V';

end