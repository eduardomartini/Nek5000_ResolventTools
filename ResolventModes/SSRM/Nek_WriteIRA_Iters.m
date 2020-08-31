function  Nek_WriteIRA_Iters(reaFile,nCurrIter,iif,X,Y,xx,yy, ...
                        lr1,elmap,freq,istep,fields,emode,wdsz,etag)
    fileloc =@(i,runtype,prefix,iif) sprintf('IterAr%02.f/%s/%s%s0.f%05.0f',i,runtype,prefix,reaFile,iif);     
    if lr1(3)>1
        ndim = 3;
    else
        ndim =2;
    end       
    for i=1:nCurrIter
        % Read iteration inputs

        file = fileloc(i-1,'dir','fc1',iif);  disp(file);
        data = reshape(real(X(:,i)),[],prod(lr1),ndim); 
        writenek(file,data,lr1,elmap,freq,istep,fields,emode,wdsz,etag );

        file = fileloc(i-1,'dir','fs1',iif);  disp(file);
        data = reshape(imag(X(:,i)),[],prod(lr1),ndim); 
        writenek(file,data,lr1,elmap,freq,istep,fields,emode,wdsz,etag );
        
        % Read iteration output
        file = fileloc(i-1,'adj','c01',iif);
        disp(file);
        data = reshape(real(Y(:,i)),[],prod(lr1),ndim); 
        writenek(file,data,lr1,elmap,freq,istep,fields,emode,wdsz,etag );

        file = fileloc(i-1,'adj','s01',iif);
        disp(file);
        data = reshape(imag(Y(:,i)),[],prod(lr1),ndim); 
        writenek(file,data,lr1,elmap,freq,istep,fields,emode,wdsz,etag );
        
        
        % Read middle normalization factor
        file = fileloc(i-1,'dir','c01',iif);
        disp(file);
        data = reshape(real(xx(:,i)),[],prod(lr1),ndim); 
        writenek(file,data,lr1,elmap,freq,istep,fields,emode,wdsz,etag );

        file = fileloc(i-1,'dir','s01',iif);
        disp(file);
        data = reshape(imag(xx(:,i)),[],prod(lr1),ndim); 
        writenek(file,data,lr1,elmap,freq,istep,fields,emode,wdsz,etag );
        
        file = fileloc(i-1,'adj','fc1',iif);
        disp(file);
        data = reshape(real(yy(:,i)),[],prod(lr1),ndim); 
        writenek(file,data,lr1,elmap,freq,istep,fields,emode,wdsz,etag );

        file = fileloc(i-1,'adj','fs1',iif);
        disp(file);
        data = reshape(imag(yy(:,i)),[],prod(lr1),ndim); 
        writenek(file,data,lr1,elmap,freq,istep,fields,emode,wdsz,etag );
        
    end
