function Tu=getEstimationKernels(HGs)
    df = HGs.freq(2) - HGs.freq(1);
    dt = HGs.t(2)    - HGs.t(1)   ;
    nt = length(HGs.freq);
    freq = HGs.freq;
     FFT = @(x) ifftshift(ifft( fftshift(    x   ,3),[],3)   ,3)*dt*nt; % fft defined for conventions used;
    IFFT = @(x) ifftshift( fft( fftshift(    x   ,3),[],3)   ,3)*df; % fft defined for conventions used;
%     IFFT = @(x)            fft(              x      ,[],3)      *df   ; %ifft defined for conventions used;

    fprintf('Computing Estimation Kernels...');clock=tic();
        
    % non-causal estimation
    nt = length(HGs.t);
    for iw=1:nt
%         Tu.nchat(:,:,iw)= HGs.ghat(:,:,iw)*inv(HGs.Ghat(:,:,iw));
        Tu.nchat(:,:,iw)= HGs.ghat(:,:,iw)*HGs.iGhat(:,:,iw);
    end
    % truncated non-causal estimation
    Tu.nc                =  IFFT(Tu.nchat);
    Tu.tnc                = Tu.nc;
    Tu.tnc(:,:,1:round(end/2)) = 0; 
    Tu.tnchat             = FFT(Tu.tnc);

    % causal estimation
    for iw = 1:nt
%         RHS_hat(:,:,iw) = HGs.ghat(:,:,iw) *inv(HGs.Gminushat(:,:,iw));
        RHS_hat(:,:,iw) = HGs.ghat(:,:,iw) *HGs.iGminushat(:,:,iw);
    end
    RHS = IFFT(RHS_hat);
    RHSm=RHS;
    RHSm(:,:,1:end/2) = 0; % ()_+ operation.
    RHShatm =FFT(RHSm);

    for iw = 1:nt
%         Tu.chat(:,:,iw) = RHShatm(:,:,iw) *inv(HGs.Gplushat(:,:,iw));
        Tu.chat(:,:,iw) = RHShatm(:,:,iw) *HGs.iGplushat(:,:,iw);
    end
    Tu.c = IFFT(Tu.chat);
%     Tu.c (:,:,end/2:end) = 0;

    Tu.t    = HGs.t;
    Tu.freq = HGs.freq;
    disp([' Done in ' num2str(toc(clock)) 's']);
    