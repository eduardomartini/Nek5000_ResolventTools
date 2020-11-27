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
        HGs.Glhat(:,:,iw) =  (DATA.RfyRfydhat(:,:,iw)+DATA.RfyRfydhat(:,:,iw)')/2 + N; % removes the non-hermian part of y_esthat
    end
    c_tmp=toc(c_tmp);fprintf(', completed in %5.2fs \n', c_tmp);
    
    fprintf('\tWiener-Hopf Fact of G'); c_tmp=tic;
    [HGs.Glminushat,HGs.Glplushat]=WienerHopfDec(DATA.freq*2*pi,conj(HGs.Glhat),omega0,[],tol); % reverse factorization
    HGs.Glminushat  = conj(HGs.Glminushat); HGs.Glplushat = conj(HGs.Glplushat);
    c_tmp=toc(c_tmp);fprintf(', completed in %5.2fs \n', c_tmp);
    
    % Compute "H"    
    if isfield(DATA,'RazdRaz')
        fprintf('\tComputing H from Actuator Direct-Adjoint Run'); c_tmp=tic;
        for iw = 1:nt
            HGs.Hlhat(:,:,iw)         =  DATA.RazdRazhat(:,:,iw) + P ;
        end        
        c_tmp=toc(c_tmp);fprintf(', completed in %5.2fs \n', c_tmp);
    else
        fprintf('\tComputing H from Actuator Impulse response'); c_tmp=tic;
        for iw = 1:nt
            HGs.Hlhat(:,:,iw)         =  DATA.Razhat(:,:,iw)'*DATA.Razhat(:,:,iw) + P ;
        end    
        c_tmp=toc(c_tmp);fprintf(', completed in %5.2fs \n', c_tmp);
    end
    
    fprintf('\tWiener-Hopf Fact of H'); c_tmp=tic;
    [HGs.Hlplushat,HGs.Hlminushat]=WienerHopfDec(DATA.freq*2*pi,HGs.Hlhat,omega0,[],tol);         
    c_tmp=toc(c_tmp);fprintf(', completed in %5.2fs \n', c_tmp);
    
    % Compute "Ray"  and "h" 
    fprintf('\tComputing Actuator-Sensor Feedback'); c_tmp=tic;
    for iw = 1:nt
        HGs.Ray(:,:,iw)  = DATA.Ray(:,:,iw);
    end
    c_tmp=toc(c_tmp);fprintf(', completed in %5.2fs \n', c_tmp);
    
    HGs.Gl      = IFFT(HGs.Glhat);
    HGs.Glplus  = IFFT(HGs.Glplushat);
    HGs.Glminus = IFFT(HGs.Glminushat);
    HGs.Hl      = IFFT(HGs.Hlhat);
    HGs.Hlplus  = IFFT(HGs.Hlplushat);
    HGs.Hlminus = IFFT(HGs.Hlminushat);
    
    if_hg = isfield(DATA,'RazdRfzRfydhat');
    if ~if_hg % if hg is not present, compute g and h
        fprintf('\tComputing g from Sensor adj-dir run'); c_tmp=tic;
        for iw = 1:nt
            HGs.Grhat(:,:,iw) =  DATA.RfzRfydhat(:,:,iw);
        end
        c_tmp=toc(c_tmp);fprintf(', completed in %5.2fs \n', c_tmp);
        
        fprintf('\tComputing h from Actuator Impulse response'); c_tmp=tic;
        for iw = 1:nt
            HGs.Hrhat(:,:,iw)         = -DATA.Razhat(:,:,iw)';
        end
        HGs.Hr      = IFFT(HGs.Hrhat);
        HGs.Gr      = IFFT(HGs.Grhat);
        c_tmp=toc(c_tmp);fprintf(', completed in %5.2fs \n', c_tmp);
    else
        fprintf('\tComputing hg term from adj-dir-adj run'); c_tmp=tic;
        for iw = 1:nt
            HGs.HrGrhat(:,:,iw) = -DATA.RazdRfzRfydhat(:,:,iw);
        end
        HGs.HrGr      = IFFT(HGs.HrGrhat);
        c_tmp=toc(c_tmp);fprintf(', completed in %5.2fs \n', c_tmp);
    end
    
    
    fprintf('\tPre inverting H, H_+-, G and G_-+'); c_tmp=tic();
    for iw = 1:nt
        HGs.iHlhat     (:,:,iw) = inv(HGs.Hlhat     (:,:,iw));
        HGs.iGlhat     (:,:,iw) = inv(HGs.Glhat     (:,:,iw));
        HGs.iHlminushat(:,:,iw) = inv(HGs.Hlminushat(:,:,iw));
        HGs.iHlplushat (:,:,iw) = inv(HGs.Hlplushat (:,:,iw));
        HGs.iGlminushat(:,:,iw) = inv(HGs.Glminushat(:,:,iw));
        HGs.iGlplushat (:,:,iw) = inv(HGs.Glplushat (:,:,iw));
    end
    c_tmp=toc(c_tmp);fprintf(', completed in %5.2fs \n', c_tmp);
    
    HGs.t=DATA.t;
    HGs.freq=DATA.freq;
    fprintf('Total Time  %5.2fs \n\n ', toc(clock));

