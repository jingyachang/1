    function [X] = InitialPoint(Gr,ClusterNumber )

    N      = Gr.sz(1);
    M      = Gr.sz(2);
    K      = Gr.sz(3);
    degree = Gr.degree; 
    weights = Gr.weights; 
    G      = Gr.edge;
    H      = sparse(N,M);
    for j  = 1:M
        H_m = sparse(G(j,:)',ones(K,1)*j,ones(K,1),N,M);
        H   = H+H_m;    
    end

    D_v    = diag(degree);
    W      = diag(weights);
    LapMatr = eye(N) - 1/2*D_v^(-1/2)*H*W*H'*D_v^(-1/2); 


    [V,~]  = eig(LapMatr);
    X     = V(:,1:ClusterNumber);
    for i = 1:ClusterNumber
        X(:,i) = X(:,i)/norm(X(:,i));
    end






