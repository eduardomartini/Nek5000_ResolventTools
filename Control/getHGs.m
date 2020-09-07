function [HGs] = getHGs(DATA,N,P,omega0,tol)
% getHGs(DATA,N,P,omega0,tol)
% Computes HGs terms from interpolated DATA.
% If DATA contains RazdRfzRfyd, h*g is computed from it, otherwise h is
% obtained from Raz and g from RfzRfyd, and multipled afterwards.
% By convention all variables ended in "hat" correspond to their frequency
% domain representation. Absensce of "hat" sufix indicates a time domain
% representation.

    disp('Computing G,g,H,h matrices, and Wiener-Hopt decmopositions...');clock=tic();
    if ~exist('tol');tol=1e-6;end
    
    ny = size(DATA.RfyRfydhat,1);
    na = size(DATA.Rayhat,2);
    
    
    nt = length(DATA.t);
    dt = DATA.t(2)-DATA.t(1);
    df = DATA.freq(2)-DATA.freq(1);
     FFT = @(x) ifftshift(ifft( fftshift(    x   ,3),[],3)   ,3)*dt*nt; % fft defined for conventions used;
    IFFT = @(x) ifftshift( fft( fftshift(    x   ,3),[],3)   ,3)*df; % fft defined for conventions used;
    
    % Compute "G"  and "g" 
    if length(N)==1 && N(1)<0 
        N=eye(ny)*max(abs(DATA.RfyRfydhat(:)))*abs(N);
    end
    if length(P)==1 && P(1)<0 
        if isfield(DATA,'Razhat')
            P=eye(na)*max(abs(DATA.Razhat(:)))^2*abs(P);
        elseif isfield(DATA,'RazdRazhat')
            P=eye(na)*max(abs(DATA.RazdRazhat(:)))*abs(P);
        else
            error('Missing Raz|RazdRaz data');
        end
    end
    
    fprintf('\tComputing G from Sensor adj-dir run'); c_tmp=tic;
    for iw = 1:nt
        HGs.Ghat(:,:,iw) =  (DATA.RfyRfydhat(:,:,iw)+DATA.RfyRfydhat(:,:,iw)')/2 + N; % removes the non-hermian part of y_esthat
    end
    c_tmp=toc(c_tmp);fprintf(', completed in %5.2fs \n', c_tmp);
    
    fprintf('\tWiener-Hopf Fact of G'); c_tmp=tic;
    [HGs.Gminushat,HGs.Gplushat]=WienerHopfDec(DATA.freq*2*pi,conj(HGs.Ghat),omega0,[],tol); % reverse factorization
    HGs.Gminushat  = conj(HGs.Gminushat); HGs.Gplushat = conj(HGs.Gplushat);
    c_tmp=toc(c_tmp);fprintf(', completed in %5.2fs \n', c_tmp);
    
    % Compute "H"    
    if isfield(DATA,'RazdRaz')
        fprintf('\tComputing H from Actuator Direct-Adjoint Run'); c_tmp=tic;
        for iw = 1:nt
            HGs.Hhat(:,:,iw)         =  DATA.RazdRazhat(:,:,iw) + P ;
        end        
        c_tmp=toc(c_tmp);fprintf(', completed in %5.2fs \n', c_tmp);
    else
        fprintf('\tComputing H from Actuator Impulse response'); c_tmp=tic;
        for iw = 1:nt
            HGs.Hhat(:,:,iw)         =  DATA.Razhat(:,:,iw)'*DATA.Razhat(:,:,iw) + P ;
        end    
        c_tmp=toc(c_tmp);fprintf(', completed in %5.2fs \n', c_tmp);
    end
    
    fprintf('\tWiener-Hopf Fact of H'); c_tmp=tic;
    [HGs.Hplushat,HGs.Hminushat]=WienerHopfDec(DATA.freq*2*pi,HGs.Hhat,omega0,[],tol);         
    c_tmp=toc(c_tmp);fprintf(', completed in %5.2fs \n', c_tmp);
    
    % Compute "Ray"  and "h" 
    fprintf('\tComputing Actuator-Sensor Feedback'); c_tmp=tic;
    for iw = 1:nt
        HGs.Ray(:,:,iw)  = DATA.Ray(:,:,iw);
    end
    c_tmp=toc(c_tmp);fprintf(', completed in %5.2fs \n', c_tmp);
    
    HGs.G      = IFFT(HGs.Ghat);
    HGs.Gplus  = IFFT(HGs.Gplushat);
    HGs.Gminus = IFFT(HGs.Gminushat);
    HGs.H      = IFFT(HGs.Hhat);
    HGs.Hplus  = IFFT(HGs.Hplushat);
    HGs.Hminus = IFFT(HGs.Hminushat);
    
    if_hg = isfield(DATA,'RazdRfzRfydhat');
    if ~if_hg % if hg is not present, compute g and h
        fprintf('\tComputing g from Sensor adj-dir run'); c_tmp=tic;
        for iw = 1:nt
            HGs.ghat(:,:,iw) =  DATA.RfzRfydhat(:,:,iw);
        end
        c_tmp=toc(c_tmp);fprintf(', completed in %5.2fs \n', c_tmp);
        
        fprintf('\tComputing h from Actuator Impulse response'); c_tmp=tic;
        for iw = 1:nt
            HGs.hhat(:,:,iw)         = -DATA.Razhat(:,:,iw)';
        end
        HGs.h      = IFFT(HGs.hhat);
        HGs.g      = IFFT(HGs.ghat);
        c_tmp=toc(c_tmp);fprintf(', completed in %5.2fs \n', c_tmp);
    else
        fprintf('\tComputing hg term from adj-dir-adj run'); c_tmp=tic;
        for iw = 1:nt
            HGs.hghat(:,:,iw) = -DATA.RazdRfzRfydhat(:,:,iw);
        end
        HGs.hg      = IFFT(HGs.hghat);
        c_tmp=toc(c_tmp);fprintf(', completed in %5.2fs \n', c_tmp);
    end
    
    
    fprintf('\tPre inverting H, H_+-, G and G_-+'); c_tmp=tic();
    for iw = 1:nt
        HGs.iHhat     (:,:,iw) = inv(HGs.Hhat     (:,:,iw));
        HGs.iGhat     (:,:,iw) = inv(HGs.Ghat     (:,:,iw));
        HGs.iHminushat(:,:,iw) = inv(HGs.Hminushat(:,:,iw));
        HGs.iHplushat (:,:,iw) = inv(HGs.Hplushat (:,:,iw));
        HGs.iGminushat(:,:,iw) = inv(HGs.Gminushat(:,:,iw));
        HGs.iGplushat (:,:,iw) = inv(HGs.Gplushat (:,:,iw));
    end
    c_tmp=toc(c_tmp);fprintf(', completed in %5.2fs \n', c_tmp);
    
    HGs.t=DATA.t;
    HGs.freq=DATA.freq;
    fprintf('Total Time  %5.2fs ', toc(clock));

