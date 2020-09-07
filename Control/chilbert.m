function  y =  chilbert(x,overSample)
    if ~exist('overSample'); overSample=20;end
    n=length(x);
    y=zeros(size(x));

    if false % tailored method (not working)
        signw= [0;ones(n*overSample/2-1,1);0;-ones(n*overSample/2-1,1)];
        signw = reshape(signw,1,1,n*overSample);
        for i=1:size(x,1)
            for j=1:size(x,2)
                ytmp =ifft(fft(x(i,j,:),n*overSample,3).*(signw),n,3);
                y(i,j,:)=reshape(ytmp,1,1,n);
            end
        end
    else % using matlab function)
    %     disp(['Oversample ' num2str(overSample)]);
        for i=1:size(x,1)
            for j=1:size(x,2)
                xr= squeeze(real(x(i,j,:)));
                xi= squeeze(imag(x(i,j,:)));
                ytmp  = imag(hilbert(xr,n*overSample))+1i*imag(hilbert(xi,n*overSample));
                y(i,j,:) = ytmp(1:n);
            end
        end
end
    