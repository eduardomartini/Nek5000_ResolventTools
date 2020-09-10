function Gamma=getControlKernels(HGs)
    
    fprintf('Computing Control Kernels...');clock=tic();
    df = HGs.freq(2) - HGs.freq(1);
    dt = HGs.t(2)    - HGs.t(1)   ;
    nt = length(HGs.t);
    freq=HGs.freq;
    
     FFT = @(x) ifftshift(ifft( fftshift(    x   ,3),[],3)   ,3)*dt*nt; % fft defined for conventions used;
    IFFT = @(x) ifftshift( fft( fftshift(    x   ,3),[],3)   ,3)*df; % fft defined for conventions used;
    
    % non-causal control
    if isfield(HGs,'hghat')
        for iw = 1:nt
            Gamma.nchat(:,:,iw)= HGs.iHhat(:,:,iw)*HGs.hghat(:,:,iw)*HGs.iGhat(:,:,iw);        
        end
    else
        for iw=1:nt
            Gamma.nchat(:,:,iw)= HGs.iHhat(:,:,iw)*HGs.hhat(:,:,iw)*HGs.ghat(:,:,iw)*HGs.iGhat(:,:,iw);
        end
    end
    
    Gamma.nc = IFFT(Gamma.nchat);
    % truncated non-causal control
    Gamma.tnc             = Gamma.nc;
    Gamma.tnc(:,:,freq<0) = 0; 
    Gamma.tnchat          = FFT(Gamma.tnc);

    % causal control
    if isfield(HGs,'hghat')
        for iw = 1:nt
            RHS_hat(:,:,iw)= HGs.iHminushat(:,:,iw)*HGs.hghat(:,:,iw)*HGs.iGminushat(:,:,iw);        
        end
    else
        for iw = 1:nt
            RHS_hat(:,:,iw) = HGs.iHminushat(:,:,iw)*HGs.hhat(:,:,iw)*HGs.ghat(:,:,iw) *HGs.iGminushat(:,:,iw);
    %         RHS_hat(:,:,iw) = inv(HGs.Hminushat(:,:,iw))*HGs.hhat(:,:,iw)*HGs.ghat(:,:,iw) *inv(HGs.Gminushat(:,:,iw));
        end
    end
    RHS = IFFT(RHS_hat);
    RHSm=RHS;
    RHSm(:,:,freq<0) = 0; % ()_+ operation.
    RHSmhat = FFT(RHSm);

    for iw = 1:nt
%         Gamma.chat(:,:,iw) = inv(HGs.Hplushat(:,:,iw))*RHSmhat(:,:,iw) *inv(HGs.Gplushat(:,:,iw));
        Gamma.chat(:,:,iw) = HGs.iHplushat(:,:,iw)*RHSmhat(:,:,iw) *HGs.iGplushat(:,:,iw);
    end
    Gamma.c = IFFT(Gamma.chat);
    
    %% Compute control considering actuator-sensor feedback
    na=size(HGs.Ray,2);    
    for iw = 1:nt
        Gamma.cfbhat  (:,:,iw) = inv(eye(na)-Gamma.chat  (:,:,iw)*HGs.Ray(:,:,iw))*Gamma.chat(:,:,iw);
        Gamma.tncfbhat(:,:,iw) = inv(eye(na)-Gamma.tnchat(:,:,iw)*HGs.Ray(:,:,iw))*Gamma.tnchat(:,:,iw);
    end
    Gamma.cfb   = IFFT(Gamma.cfbhat);
    Gamma.tncfb = IFFT(Gamma.tncfbhat);

    
    
    Gamma.t    = HGs.t;
    Gamma.freq = HGs.freq;
    
    disp(['Done in ' num2str(toc(clock)) 's']);
    disp('');
