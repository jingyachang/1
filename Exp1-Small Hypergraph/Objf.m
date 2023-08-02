function [Obj] = Objf(U,X,r,Lambda,E,X_E)
  

Obj  = sum(sum(( U'*X).^r));

Obj = Obj +Lambda*norm(X_E-E,'fro')^2;

end

