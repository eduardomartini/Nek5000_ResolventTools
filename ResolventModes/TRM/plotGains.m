clear all
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave
    nfreq=glob('gains_krylov*.dat');
else
    nfreq=length(dir('gains_krylov*.dat')); 
end

for i = 1:nfreq
    gains_krylov(:,:,i) = dlmread(sprintf('gains_krylov_%03.0f.dat',i));
end
freq = squeeze(gains_krylov(1,1,:));
gains_krylov = gains_krylov(:,2:end,:);

nIters=size(gains_krylov,2);


gains_powIter=dlmread('gains_powerIter.dat');
freq=gains_powIter(:,1);
gains_powIter=gains_powIter(:,2:end);

NormList  = dlmread('run_freqNorms.dat'); 
    NormList  = NormList(:,2:end);  
NormListf = dlmread('run_freqNorms_filtered.dat');
    NormListf = NormListf(:,2:end);




%%  Figures
f=figure; 
    plot(freq,gains_powIter(:,2:2:end));
    hold on
    a=gca;
    for i=2:2:nIters
        a.ColorOrderIndex = i  ;
        plot(freq,squeeze(gains_krylov(:,i,1)),'o');
    end
    
    leg1 = split(sprintf('Iter n=%4.1f;',(2:2:nIters)/2),';');  
    leg2 = split(sprintf('Iter n=%4.1f;',(2:2:nIters)/2),';');  
    legend(leg1{1:end-1},leg2{1:end-1})
    ylabel('Gains')
    xlabel('Frequency')
    
f=figure; 
subplot(211)
    plot(freq,NormList);
    ylabel('Norm')
    xlabel('Frequency')
    title('Before Filter')
    ylim([0,max(ylim)]);
subplot(212)
    plot(freq,NormListf);
    title('After Filter')
    ylabel('Norm')
    xlabel('Frequency')
    leg = split(sprintf('%4.1f;',0:nIters-1),';');  
    legend(leg)
    ylim([0,max(ylim)]);

