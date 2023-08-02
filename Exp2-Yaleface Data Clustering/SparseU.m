 function [U] = SparseU(G,X,weights,deg_ik)


[m,r]  = size(G);  n = max(G(:));
colN = nchoosek(r,2);
[n2, ~]   =  size(X); 


if n ~= n2
     fprintf('The vertices number n does not equal to the row number of X');
     return;
end


U      =  sparse(n,m*colN); 
v      =  [ones(colN,1)  -ones(colN,1)]; 
weight =  weights.^(1/r); 


for i = 1:m
    combM = nchoosek(G(i,:),2);
    d     = 1./deg_ik(combM);
    A_m   = sparse(combM, repmat((1:colN)'+(i-1)*colN,1, 2), d.*v*weight(i), n, m*colN);
    U     = U+A_m; 
   
end

return;