function [X_G] = X_G(X,E,m,n)

[IndNr,IndNc] = find(E);
c = sub2ind([m,n],IndNr,IndNc);
X_G = zeros([m,n]);
X_G(c) = X(c);
end

