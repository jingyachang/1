function [iter, idx,bre] = main(X0,Gr,E,ratio)
% Output: X1 --> optimal point   
%         idx --> The clustering results




%% 
    epsi = 1e-5;  
    [n,k] = size(X0);   
    G      = Gr.edge;   
    [~,r]   = size(G);  
    n2      = max(G,[],'all'); 
    if n~=n2
        fprintf('The vertex numbers in X0 and G are not equal\n');
        return;
    end

    delta = 1;         
    alpha = 0.9;          
    iterNum = 1000;      
    fbest  = 1.0000e+10;  
    l   = 0;              
    L   = 1;            

    
    deg_ik = Gr.degree.^(1/r);
    [U] = SparseU(G,X0,Gr.weights,deg_ik);

    if ratio == 0 || ratio == 0.1
        Lambda = 0.3;     
    elseif ratio == 0.2
        Lambda = 0.2;
    else
        Lambda = 0.15;
    end



%% Initialization.
     [X_E,~] = X_G(X0,E,n,k);
     [Nablaf0] =  Grad(X0,X_E,E,Lambda,r,U);   
     X1 = X0- Nablaf0;
     X1 = Proj(X1,n,k);
     X_E = X_G(X1,E,n,k);
     [Nablaf01] =  Grad(X1,X_E,E,Lambda,r,U);     
     [fr] = Objf(U,X1,r,Lambda,E,X_E);
     iter = 0;



%% While loop
    while iter < iterNum

        if norm(Nablaf01)<epsi
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
        %  beta =1e-1;
          Xk = -t0*alpha^m*Nablaf01+X1;
          Xk = Proj(Xk,n,k);
          X_E = X_G(Xk,E,n,k);


          [f1] = Objf(U,Xk,r,Lambda,E,X_E);
       if f1 <= fr-delta*alpha^m*t0*norm(Nablaf01)^2
         break;           
       end     

      end
    %  Armijo finish   
    
          X0 = X1;
          X1 = Xk;
          Nablaf0 = Nablaf01;
          [Nablaf01] = Grad(X1,X_E,E,Lambda,r,U); 
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
% fr updating .
    end


%% clustering
     [idx,bre] = spectralCluster(X1,k,1);
  
     

end 





    


