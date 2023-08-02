function [Nablaf,gradf0,G0] = Grad(X,X_E,E,Lambda,r,U)


G0 = r* U*(( U'*X).^(r-1));
G0 = G0 + 2*Lambda*(X_E-E);
Nablaf = G0-0.5*X*(X'*G0+G0'*X); 
gradf0 = -norm(Nablaf,'fro')^2;
