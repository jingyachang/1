 function [Nablaf] = Grad(X,X_E,E,Lambda,r,U)

    G0 = r* U*(( U'*X).^(r-1));
    G = G0 + 2*Lambda*(X_E-E);
    Nablaf = G-0.5*X*(X'*G+G'*X); 


 end
