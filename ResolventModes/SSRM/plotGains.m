clear all
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave
    nfreq=glob('gains_SS*0.dat');
else
    nfreq=length(dir('gains_SS*')); 
end

for i = 1:nfreq
    gains(:,:,i) = dlmread(sprintf('gains_SS_%03.0f.dat',i));
end
freq = squeeze(gains(1,1,:));
gains = gains(:,2:end,:);

nIters=size(gains,2);

subplot(211)
    plot(freq,squeeze(gains(:,end,:))','-o');
    xlabel('freq');
    ylabel('$\tilde \sigma_i$');
    
subplot(223)
    iif = 1;
    plot(1:nIters,squeeze(gains(:,:,iif)).','-o')
    xlabel('iter');
    ylabel('$\tilde \sigma_i$');
    title(sprintf('freq=%.4f',freq(iif)))
subplot(224)
    iif = 15;
    plot(1:nIters,squeeze(gains(:,:,iif)).','-o')
    xlabel('iter');
    ylabel('$\tilde \sigma_i$');
    title(sprintf('freq=%.4f',freq(iif)))
