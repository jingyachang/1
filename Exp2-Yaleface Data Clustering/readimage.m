function [PicMatrix,vecPerson,vecImage] = readimage(numPerson,numImage)

        fprintf('add Data folder to path (readimage.m)\n '); pause
        addpath(genpath('D:Data'));
       
       facename = ['+000E+00.pgm';'+000E+20.pgm';'+000E+45.pgm';'+000E-20.pgm'; ...
                    '+000E-35.pgm';'+005E+10.pgm';'+005E-10.pgm';'+010E+00.pgm'; ...
                    '+010E-20.pgm';'+015E+20.pgm';'+020E+10.pgm';'+020E-10.pgm'; ...
                    '-005E+10.pgm';'-005E-10.pgm';'-010E+00.pgm';'-010E-20.pgm'; ...
                    '-015E+20.pgm';'-020E+10.pgm';'-020E-10.pgm';'-025E+00.pgm'];
                
        vecPerson = [9	12	27	37];
        vecImage = [4	6	18	14	1	20	10	12];
        ratioPic = 0.5; 
        PicMatrix = zeros(192*168*ratioPic^2,numPerson*numImage);
        m = 192; n = 168; bigim = zeros(numPerson*m,numPerson*n); 
        
        for i = 1:numPerson
            for j=1:numImage
                person = vecPerson(i);
                pic = getpgmraw(['yaleB',num2str(person,'%02d'),'\yaleB',num2str(person,'%02d'),'_P00A',facename(vecImage(j),:)]);
                bigim((i-1)*m+(1:m),(j-1)*n+(1:n)) = pic;  
                pic = imresize(pic,ratioPic);
                PicMatrix(:,(i-1)*numImage+j) = pic(:);
            end
        end

        figure(person); imshow(bigim,[]);

end



