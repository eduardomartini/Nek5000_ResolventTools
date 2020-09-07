function [Gp,Gm,Fp]=WienerHopfDec(wlist,G,omega0,matrixFreeflag,tol)
    % Compute Wiener-Hopt Decompostion using Fredholm integral equation
    % (Daniele 2017).
    % wlist : uniformily spaced frequency vector
    % G : (n,n,nw) kernel, with n being the size of the vector and nw the
    % number of dicretized frequencies
    %  Gp, Gm : Plus and minus decompositions
    % omega0 : Fredholm parameter with negative imaginary part .
    % G = Gm * Gp

    if ~exist('matrixFreeflag')
        matrixFreeflag = true;
    elseif isempty(matrixFreeflag)
        matrixFreeflag = true;
    end
    
    if(wlist(1)==0)
        %matrix construction method
        disp('Error, first frequency in list is zero! Probably wrong ordering. Trying to fix'); 
        wlist = ifftshift(wlist);
        G = ifftshift(G,3);
        fftshitAtEnd = true;
    else
        fftshitAtEnd = false;
    end

    n  = size(G,1) ;
    nw = length(wlist);
    dw = wlist(2)-wlist(1);
    %define functions to move from vector and tensor representation
    toVec   = @(x) reshape(permute(x,[1,3,2]),n*nw,1);
    fromVec = @(x) ipermute(reshape(x,n,nw,1),[1,3,2]);

    %% Explicty Matrix Construction Method
    if ~matrixFreeflag
        fprintf('Constructing Fredmhol equation Matrices ...'); tic();
        M = zeros(n*nw,n*nw);
        for iw=1:nw
            for jw=1:nw
                i = (1:n)+n*(iw-1);
                j = (1:n)+n*(jw-1);
                if (iw==jw) 
                    if (iw>1 & iw<nw)
                        M(i,j) = dw/(2i*pi)*(G(:,:,iw+1)-G(:,:,jw-1)) / (wlist(iw+1)-wlist(jw-1)); 
                    end
                else
                    M(i,j) = dw/(2i*pi)*(G(:,:, iw )-G(:,:, jw )) / (wlist( iw )-wlist( jw ));
                end            
            end
        end
        for iw=1:nw
            i = (1:n)+n*(iw-1);
            M(i,i) = M(i,i) + G(:,:,iw);
        end
            RHS = zeros(nw*n,n);
        for iw=1:nw
                i = (1:n)+n*(iw-1);
                RHS(i,:) = eye(n)/(wlist(iw) - omega0);
        end
%         fprintf(' %4.2f seconds spent \n',toc()); 

        fprintf('Solving equation...');tic();
        Fp = M\RHS;

        % Constructing Gp and Gm
        Gp = zeros(n,n,iw);
        Gm = zeros(n,n,iw);
        for iw=1:nw
            i = (1:n)+n*(iw-1);
            Gp(:,:,iw) = inv(Fp(i,:)*(wlist(iw) - omega0));
            Gm(:,:,iw) = G(:,:,iw)*(Fp(i,:)*(wlist(iw) - omega0));
        end    
%         fprintf(' %4.2f seconds spent \n',toc());
    
        Fp2=zeros(size(G));
        for i=1:n
            Fp2(:,i,:) = fromVec(Fp(:,i));
        end
        Fp=Fp2;
    else
        %%
%         fprintf('Matrix free solution for Fredmhol equation...\n'); tic();
%         if size(G,3) > 1e3
%             disp(['Computing initial guess']);
%             wlist_red= linspace(min(wlist),max(wlist),500);
%             G_red    = interpTensors(wlist,G,wlist_red);
%             [~,~,Fp0]=WienerHopfDec(wlist_red,G_red,omega0,false);
%             Fp0      = interpTensors(wlist_red,Fp0,wlist);
%         else
            Fp0 = zeros(size(G));
%         end
        
        for i=1:nw
            Gi(:,:,i)=inv(G(:,:,i));
        end
        
        LHS = @(x) multiprod(G,x)-1/2i*chilbert(multiprod(G,x))+multiprod(G/2i,chilbert(x));

        RHS2 = zeros(n,n,nw);
        for iw=1:nw
                RHS2(:,:,iw) = eye(n)/(wlist(iw) - omega0);
        end

%       construct preconditioner
        prec = @(x) toVec(multiprod(Gi,fromVec(x)))  ;
        Fp2=zeros(size(RHS2));
        for i=1:n
            tic
            [Fp2_tmp,~]= gmres( @(x) toVec(LHS(fromVec(x))), toVec(RHS2(:,i,:)), ...
                                25,tol,25,prec,[],toVec(Fp0(:,i,:)) ) ;
            Fp2(:,i,:) = fromVec(Fp2_tmp); 
%             toc()
        end
%         plot(wlist,M*Fp,wlist,squeeze(M*toVec(Fp2)),':',wlist,squeeze(LHS(fromVec(Fp))),'--')
        Gp2 = zeros(n,n,nw);
        Gm2 = zeros(n,n,nw);
        for iw=1:nw
            i = (1:n)+n*(iw-1);
            Gp2(:,:,iw) =        inv(Fp2(:,:,iw)*(wlist(iw) - omega0));
            Gm2(:,:,iw) = G(:,:,iw)*(Fp2(:,:,iw)*(wlist(iw) - omega0));
        end
%        fprintf(' %4.2f seconds spent \n',toc());
%         plot(wlist,squeeze(Gm),wlist,squeeze(Gm2))
        %%
       
       Gm=Gm2;
       Gp=Gp2;
       

    end

    % Constructing Gp and Gm
    if(fftshitAtEnd)
        disp('Reverting function to previous frequecy ordering ...'); 
        Gp = fftshift(Gp,3);
        Gm = fftshift(Gm,3);
    end

end 

