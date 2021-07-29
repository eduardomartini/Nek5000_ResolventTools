function [freq,X,Y,varargout] = Nek_ReadIters(reaFile,iIterList,iif,meshSize)
    if (nargout==3)
        disp('Only reading X and Y')
        readxxyy=false;
    elseif (nargout==5)
        disp('Reading X, Y, xx and yy')
        readxxyy=true ;
    else
        error('Wrong number of outputs')
    end
    
    fileloc =@(i,runtype,prefix,iif) sprintf('IterAr%02.f/%s/%s%s0.f%05.0f',i,runtype,prefix,reaFile,iif); 
    nCurrIter = length(iIterList);
    X =nan(meshSize,nCurrIter);
    Y =nan(meshSize,nCurrIter);
    if readxxyy
        xx=nan(meshSize,nCurrIter);
        yy=nan(meshSize,nCurrIter);
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
        
        % Read iteration output
        file = fileloc(ii-1,'adj','c01',iif);
        disp(file);
        [data,~,~,~,~,fields,~,~,~,~,~]= readnek(file );
        Y(:,i) = reshape(data(:,:,(1:ndim)+ndim*(fields(1) == 'X')),[],1,1); 

        file = fileloc(ii-1,'adj','s01',iif);
        disp(file);
        [data,~,~,~,~,fields,~,~,~,~,~]= readnek(file );
        Y(:,i) = Y(:,i)+1i*reshape(data(:,:,(1:ndim)+ndim*(fields(1) == 'X')),[],1,1); 
        
        if readxxyy
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
            
        end

    end

    if readxxyy
        varargout{1}= xx;
        varargout{2}= yy;
    end