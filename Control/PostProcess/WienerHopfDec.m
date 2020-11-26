function [Hp,Hm,Fp]=WienerHopfDec(wlist,H,omega0,matrixFreeflag,tol,Hi)
    % Compute Wiener-Hopt Decompostion using Fredholm integral equation
    % (Daniele 2017).
    % wlist : uniformily spaced frequency vector
    % G : (n,n,nw) kernel, with n being the size of the vector and nw the
    % number of dicretized frequencies
    %  Gp, Gm : Plus and minus decompositions
    % omega0 : Fredholm parameter with negative imaginary part .
    % G = Gm * Gp
    
    
    %Swhiches between the construction of the matrix and the matrix-free
    % methods.
    if ~exist('matrixFreeflag')
        matrixFreeflag = true;
    elseif isempty(matrixFreeflag)
        matrixFreeflag = true;
    end
    
    if(wlist(1)==0)
        %matrix construction method
        disp(['Error, first frequency in list is zero!' ... 
               'Probably wrong ordering, using ifftshift to fix it.']); 
        wlist = ifftshift(wlist);
        H = ifftshift(H,3);
        fftshitAtEnd = true;
    else
        fftshitAtEnd = false;
    end

    n  = size(H,1) ;
    nw = length(wlist);
    dw = wlist(2)-wlist(1);
    
    %Define functions to move from vector (suitable for use of linear 
    %   problems solvers) and matrices (used in the rest of the code) 
    %   representations
    toVec   = @(x) reshape(permute(x,[1,3,2]),n*nw,1);
    fromVec = @(x) ipermute(reshape(x,n,nw,1),[1,3,2]);

%% Explicty Matrix Construction Method

    Fp0 = zeros(size(H));

    if ~exist("Hi")
        for i=1:nw
            Hi(:,:,i)=inv(H(:,:,i));
        end
    end

    LHS = @(x) x-1/2i*multiprod(Hi,chilbert(multiprod(H,x)))+1/2i*chilbert(x);

    RHS = zeros(n,n,nw);
    for iw=1:nw
            RHS(:,:,iw) = Hi(:,:,iw)/(wlist(iw) - omega0);
    end

    Fp=zeros(size(RHS));
    for i=1:n
        [Fp_tmp,~]= gmres( @(x) toVec(LHS(fromVec(x))), toVec(RHS(:,i,:)), ...
                            25,tol,25) ;
        Fp(:,i,:) = fromVec(Fp_tmp); 
    end

    Hp = zeros(n,n,nw);
    Hm = zeros(n,n,nw);
    for iw=1:nw
        Hp(:,:,iw) =        inv(Fp(:,:,iw)*(wlist(iw) - omega0));
        Hm(:,:,iw) = H(:,:,iw)*(Fp(:,:,iw)*(wlist(iw) - omega0));
    end

    % Constructing Gp and Gm
    if(fftshitAtEnd)
        disp('Reverting function to previous frequecy ordering ...'); 
        Hp = fftshift(Hp,3);
        Hm = fftshift(Hm,3);
    end

end 

