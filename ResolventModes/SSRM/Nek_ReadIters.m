function [freq,InputsOutputs] = Nek_ReadIters(reaFile,iIterList,iif,meshSize,readAllFiles)    
    
    fileloc =@(i,runtype,prefix,iif) sprintf('IterAr%02.f/%s/%s%s0.f%05.0f',i,runtype,prefix,reaFile,iif); 
    nCurrIter = length(iIterList);
    % Allocate memory for inputs and outputs
    X =nan(meshSize,nCurrIter);
    if readAllFiles
        Y =nan(meshSize,nCurrIter);
        xx=nan(meshSize,nCurrIter);
        yy=nan(meshSize,nCurrIter);
    else
        Y =nan(meshSize,1);
    end
    
    
    i=0;
    for ii=iIterList
        i=i+1;
        % Read iteration inputs
        file = fileloc(ii-1,'dir','fc1',iif);
        disp(file);
        [data,lr1,~,freq,~,fields,~,~,~,~,~]= readnek(file );
        if lr1(3)>1
            ndim = 3;
        else
            ndim =2;
        end       
        X(:,i) = reshape(data(:,:,(1:ndim)+ndim*(fields(1) == 'X')),[],1,1); 

        file = fileloc(ii-1,'dir','fs1',iif);
        disp(file);
        [data,~,~,~,~,fields,~,~,~,~,~]= readnek(file );
        X(:,i) = X(:,i)+1i*reshape(data(:,:,(1:ndim)+ndim*(fields(1) == 'X')),[],1,1); 
        
        if readAllFiles
            % Read iteration output
            file = fileloc(ii-1,'adj','c01',iif);
            disp(file);
            [data,~,~,~,~,fields,~,~,~,~,~]= readnek(file );
            Y(:,i) = reshape(data(:,:,(1:ndim)+ndim*(fields(1) == 'X')),[],1,1); 

            file = fileloc(ii-1,'adj','s01',iif);
            disp(file);
            [data,~,~,~,~,fields,~,~,~,~,~]= readnek(file );
            Y(:,i) = Y(:,i)+1i*reshape(data(:,:,(1:ndim)+ndim*(fields(1) == 'X')),[],1,1); 
        
            % Read middle normalization factor
            file = fileloc(ii-1,'dir','c01',iif);
            disp(file);
            [data,~,~,~,~,fields,~,~,~,~,~]= readnek(file );
            xx(:,i) = reshape(data(:,:,(1:ndim)+ndim*(fields(1) == 'X')),[],1,1); 

            file = fileloc(ii-1,'dir','s01',iif);
            disp(file);
            [data,~,~,~,~,fields,~,~,~,~,~]= readnek(file );
            xx(:,i) = xx(:,i)+1i*reshape(data(:,:,(1:ndim)+ndim*(fields(1) == 'X')),[],1,1); 

            file = fileloc(ii-1,'adj','fc1',iif);
            disp(file);
            [data,~,~,~,~,fields,~,~,~,~,~]= readnek(file );
            yy(:,i) = reshape(data(:,:,(1:ndim)+ndim*(fields(1) == 'X')),[],1,1); 

            file = fileloc(ii-1,'adj','fs1',iif);
            disp(file);
            [data,~,~,~,~,fields,~,~,~,~,~] = readnek(file );
            yy(:,i) = yy(:,i)+1i*reshape(data(:,:,(1:ndim)+ndim*(fields(1) == 'X')),[],1,1); 
            
        else %only reads last output
            if ii == iIterList(end)
                % Read i    teration output
                file = fileloc(ii-1,'adj','c01',iif);
                disp(file);
                [data,~,~,~,~,fields,~,~,~,~,~]= readnek(file );
                Y(:,1) = reshape(data(:,:,(1:ndim)+ndim*(fields(1) == 'X')),[],1,1); 

                file = fileloc(ii-1,'adj','s01',iif);
                disp(file);
                [data,~,~,~,~,fields,~,~,~,~,~]= readnek(file );
                Y(:,1) = Y(:,1)+1i*reshape(data(:,:,(1:ndim)+ndim*(fields(1) == 'X')),[],1,1); 
            end
        end

    end

    if readAllFiles
        InputsOutputs = {X,xx,yy,Y};
    else
        InputsOutputs = {X,Y};
    end