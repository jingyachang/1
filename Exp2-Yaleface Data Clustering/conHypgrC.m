
function [G,w] = conHypgrC(Y,r)


%% similarity of the trajectory vector

        [~,n] = size(Y);
        Z = NSN(Y,10,10,1e-4); 
        S_m = abs(Z)+abs(Z');


%% Hypergraph incidence matrix

        [~,I] = sort(S_m,'descend');
        G = I(1:r-1,:)';
        v = (1:n)';
        G = [v G];
        % unique rows of matrix G
        for i = 1:n
            G(i,:) = sort(G(i,:));
        end
        % delete the repetated edges and their weights
        [G, ~,~] = unique(G,'rows'); 

%% weight computation
       [m,~] = size(G);
       w = ones(m,1);
       for i = 1:m          
        g = nchoosek(G(i,:),2);
        [lg,~] = size(g);
        for j = 1:lg
           w(i) = w(i)+S_m(g(j,1),g(j,2));
        end       
       end
      
      
end
 


 


