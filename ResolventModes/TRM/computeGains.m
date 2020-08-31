nIters = length(dir('Iter*'));    
ifList = 1:size(dlmread( 'freqList.input'),1);

qhat  ={};
qhat_f={};

[bm1,lr1,elmap,~,istep,fields,emode,wdsz,etag,header,status] = readnek(['bm1' reaFile '0.f00001']);
if lr1(3)>1
    ndim = 3;
else
    ndim =2;
end
XY = bm1(:,:,1:ndim);
M = bm1(:,:,ndim+1);
clear bm1
freq=[];

gains_powIter    = nan(nIters,length(ifList));
gains_powIter_da = nan(nIters,length(ifList));
gains_krylov  = nan(5,nIters,length(freq));
NormList  = nan(nIters,length(ifList));
NormListf = nan(nIters,length(ifList));

for iif = 1:length(ifList)
    for i=0:nIters-1
        file = sprintf('Iter%02.0f/%s%s0.f%05.0f',i,'c01',reaFile,ifList(iif));
        disp(file);
        [data,lr1,elmap,freq(iif),istep,fields,emode,wdsz,etag,header,status] = readnek(file );
        qhat{i+1,iif} = data(:,:,(1:ndim)+ndim*(fields(1) == 'X')); 
        
        file = sprintf('Iter%02.0f/%s%s0.f%05.0f',i,'s01',reaFile,ifList(iif));
        disp(file)
        [data,lr1,elmap,freq(iif),istep,fields,emode,wdsz,etag,header,status] = readnek(file );
        qhat{i+1,iif} = qhat{i+1,iif} + data(:,:,(1:ndim)+ndim*(fields(1) == 'X'))*1i; 

        file = sprintf('Iter%02.0f/FIR/%s%s0.f%05.0f',i,'c01',reaFile,ifList(iif));
        disp(file);
        [data,lr1,elmap,freq(iif),istep,fields,emode,wdsz,etag,header,status] = readnek(file );
        qfhat{i+1,iif} = data(:,:,(1:ndim)+ndim*(fields(1) == 'X')); 
        file = sprintf('Iter%02.0f/FIR/%s%s0.f%05.0f',i,'s01',reaFile,ifList(iif));
        disp(file)
        [data,lr1,elmap,freq(iif),istep,fields,emode,wdsz,etag,header,status] = readnek(file );
        qfhat{i+1,iif} = qfhat{i+1,iif} + data(:,:,(1:ndim)+ndim*(fields(1) == 'X'))*1i; 
        
        % Fix phase (adjoint run running back in time)
%         if mod(i,2)==1
%             qhat{i+1,iif} = - conj(qhat{i+1,iif});
%             qfhat{i+1,iif}= - conj(qfhat{i+1,iif});
%         end
    end


    %% Power Iteration 
    IP = @(x,y) sum(sum(sum(conj(x).*y.*repmat(M,1,1,ndim))));
    N  = @(x) sqrt(IP(x,x));
    costheta = @(x,y)  IP(x,y)/(N(x)*N(y));

    for i=1:nIters
        NormList (i,iif) = N(qhat {i,iif});
        NormListf(i,iif) = N(qfhat{i,iif});
    end
    gains_powIter_da(2:end,iif) = NormList(2:end,iif)./NormListf(1:end-1,iif);
    gains_powIter(3:end,iif) = sqrt(gains_powIter_da(2:end-1,iif).*gains_powIter_da(3:end,iif));
    
    MM = repmat(M,1,1,ndim);
    MM = MM(:);
    mm = sqrt(MM);  

    clear XX YY X Y        
    for i = 1:nIters-2
       scale(i) = IP(qhat{i+1,iif},qfhat{i+1,iif})/N(qhat{i+1,iif})^2;
       X(:,i)   = qfhat{i  ,iif}(:);
       Y(:,i)   = qhat{i+2,iif}(:)/scale(i);
    end
    
    for i=2:2:nIters-2
%         [Xn,R] = gson(X(:,(1:2:i)).*mm);
        [Xn,R] = qr(X(:,(1:2:i)).*mm,0);
        Yn = (Y(:,1:2:i).*mm)*inv(R);
        H = Xn'*Yn ;
        n=size(H,2);
        [psi,gains_tmp]= eig(H);
        [gains_tmp,order]= sort(sqrt(diag(gains_tmp)),'descend'); 
        gains_krylov(1:min(5,n),i+2,iif) = gains_tmp(1:min(5,n));
    end
end

dlmwrite('gains_powerIter.dat',[freq;gains_powIter].');
for i=1:size(gains_krylov,1)
    dlmwrite(sprintf('gains_krylov_%03.0f.dat',i),[freq;squeeze(gains_krylov(1,:,:))].');
end
dlmwrite('run_freqNorms.dat',[freq;NormList].');
dlmwrite('run_freqNorms_filtered.dat',[freq;NormListf].');



