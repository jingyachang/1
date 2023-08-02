
 
	r = 4;  
    numPerson = 4;
    numImage = 8;
    [PicMatrix,vecPerson,vecImage] = readimage(numPerson,numImage);
    numVer = numPerson*numImage;  
    label = zeros(numVer,1); 

    for i = 1:numPerson    
            label((i-1)*numImage+1:i*numImage) = label((i-1)*numImage+1:i*numImage)+i;    
    end

    [edges,weights] = conHypgrC(PicMatrix,r);
    Gr = ReadGraphC(edges);



%%  
    ratioNum = 4;
    numRun = 200;
    ratioVector = 0.1*(0:1:ratioNum-1);
    m = numPerson;          

    Error = zeros(numRun,ratioNum);
    iterNum = zeros(numRun,ratioNum);
    runTime = zeros(numRun,ratioNum);


    for i = 1:numRun
        X0 = rand(numVer,m);
        for t = 1:ratioNum
           
            k = randsample(numVer,round(numVer*ratioVector(t))); 
            s = label(k);
            d = ones(round(numVer*ratioVector(t)),1); 
            E = sparse(k,s,d,numVer,m);
            E = full(E);

             tic;  
            [iter,idx,bre] = main(X0,Gr,E,ratioVector(t));   
            if bre == 1
                Error(i,t) = 0;
                runTime(i,t) = 0;
                iterNum(i,t) = 0;
            else
                Error(i,t) = computeCE(idx,label);
                runTime(i,t) = toc;
                iterNum(i,t) = iter;
            end

         end

    end


    averError = sum(Error,1)/numRun;
    averIter = round(sum(iterNum,1)/numRun);
    averTime = sum(runTime,1)/numRun;

%%
    fprintf('    ratio      averError       averIter        averTime\n');
    for t = 1:ratioNum
        fprintf('     %2.1f       %8.5f         %d          %8.5f\n', ratioVector(t), averError(t),averIter(t),averTime(t));
    end







