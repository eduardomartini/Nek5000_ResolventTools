
parameters = dlmread('params.input');
dt         = parameters(2)*parameters(3);
cutOff     = parameters(4);
df         = parameters(5);

minNpoint = 1/(cutOff*dt);

N  = round(1/(dt*df));
fs = 1/dt;  % Sampling Frequency

%create filter
F_A = dlmread(ampListFile);
F       = F_A(:,1);
amps    = F_A(:,2);


if F(end)<fs/2
    x = [F(1:end-1).',linspace(F(end),fs/2,20)];
else
    x = F(:).';
end

if x(end)>fs/2
      disp('Frequency vector has frequencies above the nyquist frequency!!! ');
      disp(['F_nyqst =', num2str(fs/2), ', max = ', num2str(max(F))]);
      disp('Aborting')
      exit
end

y = x*0; 
y(1:length(amps)) = amps; 
y(length(amps):end) = amps(end); 

penalty = (tanh(tan(pi*(x-(cutOff+x(end))/2)/(cutOff-x(end))))+1)*0.5;
penalty(end)=0;
penalty(x<cutOff)=1;

hfilt = fir2 (N, x/(fs/2),penalty./y ); %create filter

dlmwrite('FIR_Coefs.txt',hfilt.');

