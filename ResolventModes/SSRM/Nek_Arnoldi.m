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
end


sig= nan(nfreqs,nCurrIter);
for iif =  iFreqsList
    %%Load all Relevant Files
    [freq,X,Y,xx,yy] = Nek_ReadIters(reaFile,1:nCurrIter,iif,numel(XY));
    % Scales input to norm 1, and correct for normalization in the adjoint run
    scale = 0;
    for i=1:nCurrIter
        %normalize input
        scale = NORM(X(:,i));
        X(:,i)=X(:,i)/scale;
        Y(:,i)=Y(:,i)/scale;
        
        %correct normalization
        scale = NORM(xx(:,i))/NORM(yy(:,i));
        Y(:,i)=Y(:,i)*scale;
    end
        
    % Computes residual.Component of the response ortogonal to previous
    % inputs
    fk = Y(:,end)- X*(pinv(X.*m)*(Y(:,end).*m));
    fk = fk/NORM(fk);
    
    Hk = IP(X,Y);
    
    % Prepares next run
    fields(1:3)='U  ';
    filelocOut = sprintf('ForceFiles/extHarmForceCos0.f%05.0f',iif);         
    disp(['Writting ' filelocOut])
    writenek(filelocOut,reshape(real(fk ),size(XY,1),size(XY,2),size(XY,3)), ...
                        lr1,elmap,freq,istep,fields,emode,wdsz,etag);
    filelocOut = sprintf('ForceFiles/extHarmForceSin0.f%05.0f',iif);         
    disp(['Writting ' filelocOut])
    writenek(filelocOut,reshape(imag(fk ),size(XY,1),size(XY,2),size(XY,3)), ...
                        lr1,elmap,freq,istep,fields,emode,wdsz,etag);
    
    
    % Compute Gains for all Iterations, and save modes found for the last iteration
    if nOutputModes>0
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
            RespModes (:,i) = RespModes(:,i) / NORM(RespModes(:,i));
            ForcModes (:,i) = ForcModes(:,i) / NORM(ForcModes(:,i));
        end

        for i=1:nOutputModes
                fields(1:3)='XU ';

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
