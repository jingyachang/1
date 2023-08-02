


    %% 
    r0          = 2;  
    cNum        = 2;
    r1          = r0+cNum; 
    numEdgesub = 5; 
    numGraphsub = 4; 
    numEdgerand = 7;
    vector = 1:r0*numEdgesub*numGraphsub;
    a = reshape(vector,[r0,numEdgesub*numGraphsub])';
    b = ones(numEdgesub*numGraphsub,cNum);
    for i = 1:numGraphsub
        for j = 1:cNum
        b((i-1)*(numEdgesub)+1:i*numEdgesub,j) = b((i-1)*(numEdgesub)+1:i*numEdgesub,j)+r0*numEdgesub*numGraphsub+i-1+(j-1)*numGraphsub;
        end
    end
    b = [a b];
    numVer = max(b(:));
    c = [36    29     9    44;
        13    33     6     1;
        48    28    37    41;
         5    18    21     4;
         2    17    45    27;
        24    16    46    47;
        15     3    34    42; ];

    G = [b;c];
    label = ones(numVer,1);

    for i = 1:numGraphsub
        label((i-1)*r0*numEdgesub+1:i*r0*numEdgesub,1) = label((i-1)*r0*numEdgesub+1:i*r0*numEdgesub,1)+i-1;
        for j = 1:2
        label(r0*numEdgesub*numGraphsub+i+numGraphsub*(j-1)) = i;
        end
    end

    weights = ones(numEdgesub*numGraphsub+numEdgerand,1);
    Gr = ReadGraphC(G,weights);
    %%  
    m = numGraphsub;    
    ratio1 = 0.1;       
    runtime = 200;      

    Error1 = zeros(runtime,1);
    Error2 = Error1;
    X_1 = zeros(numVer,runtime*numGraphsub);
    X_2 = zeros(numVer,runtime*numGraphsub);
    idx_1 = zeros(numVer,runtime);
    idx_2 = zeros(numVer,runtime);
    
    
    for i = 1:runtime
    
    X0 = InitialPoint(Gr,m );
    k = randsample(numVer,round(numVer*ratio1)); 
    s = label(k);
    d = ones(round(numVer*ratio1),1); 
    E1 = sparse(k,s,d,numVer,m);
    E1 = -full(E1)*0.1;
    E0 = zeros(numVer,m);

    
    idx1 = main(X0,Gr,E1); 
    Error1(i) = computeCE(idx1,label);
    idx_1(:,i) = idx1;

    idx2 = main(X0,Gr,E0); 
    Error2(i) = computeCE(idx2,label);
    idx_2(:,i) = idx2;

    end

    max1 = max(Error1);   min1 = min(Error1);
    max2 = max(Error2);   min2 = min(Error2);

    ave1 = sum(Error1)/runtime; 
    ave2 = sum(Error2)/runtime;


    med1 = median(Error1);
    med2 = median(Error2);

%%
    
    fprintf(' Ratio        max            min           average        median\n')
    fprintf(' %.1f       %f        %f       %f       %f\n', 0, max2, min2,ave2,med2);
    fprintf(' %.1f       %f        %f       %f       %f\n',0.1,  max1, min1,ave1,med1);

  




