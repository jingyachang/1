 function Gr = ReadGraphC(str,weights)

 if ischar(str)
    G = load(str);
elseif isnumeric(str)
    G = str;
else
    error('Input STR is invalid.');
end
[M,K] = size(G);  N = max(G(:));
if nargin == 1
    weights = ones(M,1);
else
    weights = weights(:);
end



%% Compute degree and return the hypergraph
degree = full(sum(sparse(G,repmat((1:M)',1,K),repmat(weights,1,K),N,M), 2)); % degree vector

deg_iK = degree.^(1/K);


Gr.sz      = [N,M,K];
Gr.e       = deg_iK;
Gr.edge    = G;
Gr.weights = weights;
Gr.degree  = degree;

return;
