isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

%reads parameters
parameters = dlmread('params.input');
nfreqs = parameters(4);

%checks for parallel run, otherwize run for all frequencies
if exist('iProc') & exist('nProc') 
      iFreqsList=iProc:nProc:nfreqs;
else
      iFreqsList=1:nfreqs;
end
disp(['Frequency list to be processed:' num2str(iFreqsList) ]);

%reads integration weights 
[bm1,lr1,elmap,~,istep,fields,emode,wdsz,etag,header,status] = readnek([ 'bm1' reaFile '0.f00001']);
if lr1(3)>1
    ndim = 3;
else
    ndim =2;
end
XY = bm1(:,:,1:ndim);
M = bm1(:,:,(1:ndim)+ndim);
M=M(:);
m=sqrt(M);

IP = @ (X,Y) X'* (Y.*M);
NORM = @(X) sqrt(IP(X,X));

clear bm1
freq=[];

%Gets nulber of iterations run so far
if isOctave
    nCurrIter=length(glob('IterAr*'));    
else
    nCurrIter=length(dir('IterAr*'));
end
 
%number of modes to be saved
if ~exist('nOutputModes')
    nOutputModes=0;
else
    if nOutputModes<0
        nOutputModes=nCurrIter;
    end
    nOutputModes=min(nOutputModes,nCurrIter);
end


sig= nan(nfreqs,nCurrIter);
for iif =  iFreqsList
    %%Load all Relevant Files
    readAllFiles = nOutputModes > 0; %only read all files if gains and modes are requested
    [freq,InputsOuputs] = Nek_ReadIters(reaFile,1:nCurrIter,iif,numel(XY),readAllFiles);
    if readAllFiles
        X  = InputsOuputs{ 1 };
        xx = InputsOuputs{ 2 };
        yy = InputsOuputs{ 3 };
        Y  = InputsOuputs{ 4 };
    else
        X  = InputsOuputs{1  };
        Y  = InputsOuputs{ 2 };
    end
    clear InputsOuputs
    
    % Scales input to norm 1, and correct for normalization in the adjoint run
    % Not needed if gains are not requested
    if readAllFiles
        scale = 0;
        for i=1:nCurrIter
            %normalize input
            scale = sqrt((X(:,i)'*(M.*X(:,i))) );
            X(:,i)=X(:,i)/scale;
            Y(:,i)=Y(:,i)/scale;

            %correct normalization
            scale = sqrt(  (xx(:,i)'*(M.*xx(:,i)))   /    ( yy(:,i)'*(M.*yy(:,i)) )  );
            Y(:,i)=Y(:,i)*scale;
        end
    end
    
    %% Pre multiply by the square of the norm.
    Xm = X.*m;
    Ym = Y.*m;
    % Computes Arnoldi residual (readallfiles=True). 
    % Gets the component of last output which is ortogonal to the inputs (readallfiles=false). 
    fk = Ym(:,end);
    fk = fk - Xm*(pinv(Xm)*fk);
    fk=fk/sqrt(fk'*fk);

    %Prepares inputs for next run
    fields(1:3)='U  ';
    fexp = (fk./ m) ;
    filelocOut = sprintf('ForceFiles/extHarmForceCos0.f%05.0f',iif);         
    writenek(filelocOut,reshape(real(fexp ),size(XY,1),size(XY,2),size(XY,3)), ...
                        lr1,elmap,freq,istep,fields,emode,wdsz,etag);
    filelocOut = sprintf('ForceFiles/extHarmForceSin0.f%05.0f',iif);         
    writenek(filelocOut,reshape(imag(fexp ),size(XY,1),size(XY,2),size(XY,3)), ...
                        lr1,elmap,freq,istep,fields,emode,wdsz,etag);

    if readAllFiles
        %Computes Arnoldi factorization matrix H
        Vk = Xm;
        Hk = Vk'*Ym;
    
        % Compute Gains for all Iteratins, and modes for the last iteration
        sig= nan(nCurrIter,nCurrIter);
        for j=1:(nCurrIter)
            [psi,S] = eig(Hk(1:j,1:j)); 
            [S,order]=sort(diag(S),'descend');
            sig(1:j,j) = sqrt(S(1:j));
        end
        psi=psi(:,order);
        dlmwrite(sprintf('gains_SS_%03.0f.dat',iif),[ones(nCurrIter,1)*freq,sig(:,:)])

        mkdir('Resolvent');
        RespModes = xx*psi;
        ForcModes =  X*psi;

        for i=1:nCurrIter
            RespModes (:,i) = RespModes(:,i) / sqrt(RespModes(:,i)'*(M.*RespModes (:,i)));
            ForcModes (:,i) = ForcModes(:,i) / sqrt(ForcModes(:,i)'*(M.*ForcModes (:,i)));
        end

        for i=1:nOutputModes
                fields(1:3) = 'XU ';
                data = XY*0;            
                data(:) = real(ForcModes(:,i));
                filelocOut = sprintf('Resolvent/resforce_%02.0f_real0.f%05.0f',i,iif);         
                writenek(filelocOut,cat(3,XY,data), ...
                                        lr1,elmap,freq,istep,fields,emode,wdsz,etag);
                data = XY*0;            
                data(:) = real(RespModes(:,i));
                filelocOut = sprintf('Resolvent/resresp_%02.0f_real0.f%05.0f',i,iif);         
                writenek(filelocOut,cat(3,XY,real(data)), ...
                                    lr1,elmap,freq,istep,fields,emode,wdsz,etag);

                data(:) = imag(ForcModes(:,i));
                filelocOut = sprintf('Resolvent/resforce_%02.0f_imag0.f%05.0f',i,iif);         
                writenek(filelocOut,cat(3,XY,data), ...
                                        lr1,elmap,freq,istep,fields,emode,wdsz,etag);

                data(:) = imag(RespModes(:,i));
                filelocOut = sprintf('Resolvent/resresp_%02.0f_imag0.f%05.0f',i,iif);         
                writenek(filelocOut,cat(3,XY,data), ...
                                    lr1,elmap,freq,istep,fields,emode,wdsz,etag);

        end
    end
end
