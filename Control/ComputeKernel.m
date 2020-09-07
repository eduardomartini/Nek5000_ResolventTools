%% Script to compute optimal control using Wiener-Hopf formalism, based on a DNS readings
% All frequency dependent matrices and vectors are stored with the
% frequency being the last index.
clear all

%% Control Law Parameters
% shape indexes for sensors, actuators and targets. Leave targets empy for
% full rank
iy = [ 2  ]; ny = length(iy);
ia = [ 28 ]; na = length(ia);
iz = [  5  ]; nz = length(iz);
%Sensor noise and actuation penalty : positive values for absolute value,
%negative to relative values
N=-1e-2;
P=1e-2;
% Discretization used for control Law. Input Tmax and dt. 
Tmax = 800;
dt   = 0.005;

% Other flags
actuatorDataFrom='adjoint';% adjoint : target and sensor adjoint runs, direct : actuator direct runs
computeTargetsCSD = false;
computeGains      = false;

%% 

initiClock=tic();

computeGains=true;

% Set up some functions to help plotting
sortX = @(x) ifftshift( x-(x>x(end)/2)*x(end)); % Makes a [0,x] array into a [-x/2,x/2].
sortY = @(x) ifftshift( x,2);
pfun = @(x) reshape(x,size(x,1)*size(x,2),size(x,3)); % Function for plot Tensors
getFreqVec = @(t) [((0:length(t)-1) - (0:length(t)-1 > (length(t)-1)/2 )*length(t))/(t(end)-t(1)) ];
%% Defines indexes for sensors, actuators and targets (corresponding to shape index in the Nek code

% target CSDS can only be computed if adjoint-direct target data is used.
computeTargetsCSD = computeTargetsCSD &  strcmp(actuatorDataFrom,'adjoint');
% Gains can only be computed if direct-adjoint actuator data is used.
computeGains      = computeGains  & ~strcmp(actuatorDataFrom,'adjoint');

folder = './'; %folder where the results are stored.

nt   = round(Tmax/dt); 
Tmax = nt*dt;
t    = (0:nt-1)*dt;
ts   = t - nt/2*dt;
df   = 1/t(end);
freq = getFreqVec (t);
freqS= fftshift(freq);

%% Loads Data

%reads data from adj-dir-adj run for sensors, and adj-dir run for targets.
if ~isempty(iz)
    [DATA]= readData1(ts,folder,iy,ia,iz,computeTargetsCSD,actuatorDataFrom);
else
    [DATA]= readData2(ts,folder,iy,ia);
end

%%
if ~isempty(iz)
    figure;
        subplot(221)
            plot(DATA.freq,abs(pfun(DATA.RfyRfydhat)),'-b',DATA.freq,abs(pfun(DATA.RfzRfydhat)),'-r')
            xlim([-1,1]*1)
            xlabel('$\omega$');ylabel('CSDs');
        subplot(222)
            plot(DATA.t,(pfun(DATA.RfyRfyd)),'-b',DATA.t,(pfun(DATA.RfzRfyd)),'-r')
            xlim([-10,100])
            xlabel('$t$');ylabel('Sensor/Target');
        subplot(223)
            plot(DATA.freq,abs(pfun(DATA.Rayhat)),'-b',DATA.freq,abs(pfun(DATA.Razhat)),'-r')
            xlim([-1,1]*1)
            xlabel('$\omega$');ylabel('Sensor/Target');
        subplot(224)
            plot(DATA.t,(pfun(DATA.Ray)),'-b',DATA.t,(pfun(DATA.Raz)),'-r')
            xlim([-10,100])
            xlabel('$t$');ylabel('Sensor/Target');
else
    figure;
        subplot(221)
            yyaxis left 
                plot(DATA.freq,abs(pfun(DATA.RfyRfydhat)),'b-',DATA.freq,abs(pfun(DATA.RazdRazhat)),'b:')
            yyaxis right 
                plot(DATA.freq,abs(pfun(DATA.Rayhat)));
            xlim([-1,1]*1)
            legend('$CSD_{yy}$','$R_{ad}^\dagger R_{ad}$','Act Imp.Resp');
            xlabel('$\omega$');ylabel('');
        subplot(222)
            yyaxis left 
                plot(DATA.t,(pfun(DATA.RfyRfyd)),'b-',DATA.t,(pfun(DATA.RazdRaz)),'b:')
            yyaxis right 
                plot(DATA.t,(pfun(DATA.Ray)));
            legend('$CSD_{yy}$','$R_{ad}^\dagger R_{ad}$','Act Imp.Resp');
            xlim([-1,1]*20)
            xlabel('$t$');ylabel('');
        subplot(223)
            plot(DATA.freq,abs(pfun(DATA.RazdRfzRfydhat)),'b')
            xlim([-1,1]*1)
            xlabel('$\omega$');ylabel('$R_{az}^\dagger R_{fz}R_{fy}^\dagger$');
        subplot(224)
            plot(DATA.t,(pfun(DATA.RazdRfzRfyd)),'b')
            xlim([-1,1]*20)
            xlabel('$t$');ylabel('$R_{az}^\dagger R_{fz}R_{fy}^\dagger$');
end
%% Get matrices
HGs = getHGs(DATA,-1e-2,-1e-2,-1i,1e-6); 
%%
figure;
    subplot(211)
        plot(HGs.freq,abs(pfun(HGs.Hhat)),'k',HGs.freq,abs(pfun(HGs.Hminushat)),'b',HGs.freq,abs(pfun (HGs.Hplushat)),'r')
        xlim([-1,1]*2)
        xlabel('$\omega$');ylabel('$|\hat{H}_\pm|$');
    subplot(212)
        yyaxis left
        plot(HGs.t,pfun(HGs.Hminus),'-',HGs.t,pfun(HGs.H),'k-')
        ylabel('$|\hat{H}_-|$')
        yyaxis right
        plot(HGs.t,pfun (HGs.Hplus),'-')
        ylabel('$|\hat{H}_+|$')
        xlim([-1,1]*100)
        xlabel('t');
figure
    subplot(211)
        plot(HGs.freq,abs(pfun(HGs.Ghat)),'k',HGs.freq,abs(pfun(HGs.Gminushat)),'b',HGs.freq,abs(pfun (HGs.Gplushat)),'r:')
        xlim([-1,1]*1)
        xlabel('$\omega$');ylabel('$|\hat{G}_\pm|$');
    subplot(212)
        yyaxis left
        plot(HGs.t,(pfun(HGs.Gminus)),'-',HGs.t,(pfun(HGs.G)),'k-')
        ylabel('$|\hat{G}_-|$');
        yyaxis right
        plot(HGs.t,(pfun (HGs.Gplus)),'-')
        ylabel('$|\hat{G}_+|$');
        xlim([-1,1]*10)
        xlabel('t');
%% Estimation Kernels
if isfield(HGs,'ghat')
    Tu= getEstimationKernels(HGs);
    % Plot Estimation Kernels
    figure;
        subplot(211)
            plot(Tu.freq,abs(pfun(Tu.nchat)),'-b'); hold on;
            plot(Tu.freq,abs(pfun(Tu.tnchat)),':r');
    %         plot(Tu.freq,abs(pfun(Tu.chat)),'--k');
            xlim([-1,1]*2)
            xlabel('$\omega$');ylabel('$|H|$');
        subplot(212)
            plot(Tu.t,(pfun(Tu.nc)),'-b'); hold on;
            plot(Tu.t,(pfun(Tu.tnc)),':r');
    %         plot(Tu.t,(pfun(Tu.c)),'--k');
    %         xlim([-1,1]*2)
            xlabel('$t$');ylabel('$|H|$');
end
%% Compute Control Kernels
Gamma = getControlKernels(HGs);
%%
figure;
    subplot(211)
        plot(Gamma.freq*2*pi,abs(pfun(Gamma.nchat)),'b',Gamma.freq*2*pi,abs(pfun (Gamma.tnchat)),'r:',Gamma.freq*2*pi,abs(pfun (Gamma.chat)),'k--')
        xlim([-1,1]*2)
        title('Control Kernels');
        xlabel('$\omega$');ylabel('$|\hat{\Gamma}_u|$');
    subplot(212)
        plot(Gamma.t,(pfun(Gamma.nc)),'b',Gamma.t,(pfun(Gamma.tnc)),'r:',Gamma.t,(pfun(Gamma.c)),'k--')
        xlim([-10,50])
        xlabel('t');ylabel('$\Gamma$');
        title('Control Kernels');
        
fig_ker=figure;
fig_ker.Position(3:4)=[1400,400];
        plot(Gamma.t,(pfun(Gamma.c)),'k-',Gamma.t,(pfun(Gamma.nc)),'--r',Gamma.t,(pfun(Gamma.tnc)),'g:'); 
        hold on
%         plot(Gamma.t,(pfun(Gamma.cfb)),'k--',Gamma.t,(pfun(Gamma.tncfb)),'g:')
        xlim([-10,50])
        xlabel('$\tau$');ylabel('$\Gamma$');
        legend('Causal','Non-Causal','TNC')
        
fig_kerfb=figure;
fig_kerfb.Position(3:4)=[1400,400];
        plot(Gamma.t,(pfun(Gamma.cfb)),'k-'); 
        hold on
%         plot(Gamma.t,(pfun(Gamma.cfb)),'k--',Gamma.t,(pfun(Gamma.tncfb)),'g:')
        xlim([-10,50])
        xlabel('$\tau$');ylabel("$\Gamma'$");
        legend('Causal','Non-Causal','TNC')
        
        
        
filename=['Kernel_iy=' strrep(num2str(iy),' ','_') '_ia=' strrep(num2str(ia),' ','_') '_iz=' strrep(num2str(iz),' ','_') ];
saveas(fig_ker,[filename '.fig']);
saveas(fig_ker,[filename '.eps'],'epsc');

        
%% Estimate performance
if computeTargetsCSD && ~isempty(iz)
    Czz_nc=zeros(nz,nz,nt);
    Czz_tnc=zeros(nz,nz,nt);
    Czz_c=zeros(nz,nz,nt);
    psd_un=zeros(nt,1);
    psd_nc=zeros(nt,1);
    psd_tnc=zeros(nt,1);
    psd_c=zeros(nt,1);

    for i=1:nt
        Cz1z1 = DATA.RfzRfzdhat(:,:,i);
        H=HGs.Hhat(:,:,i);
        h=HGs.hhat(:,:,i);
        G=HGs.Ghat(:,:,i);
        g=HGs.ghat(:,:,i);

        psd_un(i)=trace(Cz1z1);
        gamma = Gamma.nchat(:,:,i);
        Czz_nc(:,:,i) = Cz1z1 + h'*gamma*G*gamma'*h-h'*gamma*g'-g*gamma'*h;
        psd_nc(i)=trace(Czz_nc(:,:,i));

        gamma = Gamma.tnchat(:,:,i);
        Czz_tnc(:,:,i) = Cz1z1 + h'*gamma*G*gamma'*h-h'*gamma*g'-g*gamma'*h;
        psd_tnc(i)=trace(Czz_tnc(:,:,i));

        gamma = Gamma.chat(:,:,i);
        Czz_c(:,:,i) = Cz1z1 + h'*gamma*G*gamma'*h-h'*gamma*g'-g*gamma'*h;
        psd_c(i)=trace(Czz_c(:,:,i));
    end
    
    figrms=figure
        plot(2*pi*freqS,psd_un /(2*pi),'k',...
             2*pi*freqS,psd_c  /(2*pi),'-b', ...
             2*pi*freqS,psd_nc /(2*pi),':r', ...
             2*pi*freqS,psd_tnc/(2*pi),'--g');
         legend('Uncontroled','Causal','Non-Causal','TNC','Location','Southwest');
         xlabel('$\omega$')
         ylabel('PSD')
         xlim([-1,1]*2)

    filename=['CSD_iy=' strrep(num2str(iy),' ','_') '_ia=' strrep(num2str(ia),' ','_') '_iz=' strrep(num2str(iz),' ','_') ];
    saveas(figrms,[filename '.fig']);
    saveas(figrms,[filename '.eps'],'epsc');

    rms_un =trapz(2*pi*freqS,psd_un /(2*pi));
    rms_c  =trapz(2*pi*freqS,psd_c  /(2*pi));
    rms_nc =trapz(2*pi*freqS,psd_nc /(2*pi));
    rms_tnc=trapz(2*pi*freqS,psd_tnc/(2*pi));

    report= [ ...
    sprintf('---------- Performances ------\n' ) ... 
    sprintf('%9s\t%9s\t%9s\t%9s\n','uncon','causal','non-causal', 'tnc'), ...
    sprintf('%9f\t%9f\t%9f\t%9f\n',[rms_un,rms_c,rms_nc,rms_tnc]), ...
    sprintf('%9f\t%9f\t%9f\t%9f\n',[rms_un,rms_c,rms_nc,rms_tnc]/rms_un*100)];


    filename=['RMS_iy=' strrep(num2str(iy),' ','_') '_ia=' strrep(num2str(ia),' ','_') '_iz=' strrep(num2str(iz),' ','_') ];
     fid = fopen([filename '.txt'],'w');
    fprintf(fid ,'%s',report)
    fprintf('%s\n',report);

end
%% Export Control Law
dt_ext=0.05;
t_ext= 0:dt_ext:50;
CL    = interpTensors(Gamma.t,Gamma.c,t_ext);

formatArray = @(x) strrep(num2str(x),'  ', '_');

suffix = ['iy_' formatArray(iy) '_ia_' formatArray(ia) '_iz_' formatArray(iz) ];
fname_kerInput=['ControlKernel_' suffix   '.input'];

dlmwrite(fname_kerInput,[ny na length(t_ext) dt_ext],' ');    

for i = 1:na
    fname_kerData = sprintf('ControlKernel_%s.data%02.0f',suffix,i);
    dlmwrite(fname_kerInput,fname_kerData,'delimiter','','-append');
    CLi= reshape(CL(i,:,:),ny,length(t_ext));
    dlmwrite(fname_kerData,real(CLi.'),'delimiter',' ','precision','%.6f');
end
save(['Kernels' suffix '.mat'], 'Gamma','HGs','iy','ia','iz');

%%
disp(['Total time: ' num2str(toc(initiClock)) 's']);


figure
plot(Gamma.freq*2*pi,abs((pfun(Gamma.nchat))),'k', ...
     Gamma.freq*2*pi,abs((pfun(Gamma.chat-Gamma.nchat))), ...
     Gamma.freq*2*pi,abs((pfun(Gamma.tnchat-Gamma.nchat)))); 
 legend('$|\hat \Gamma_{c}|$','$|\hat \Gamma_{nc}-\hat \Gamma_{c}|$','$|\hat \Gamma_{tnc}-\hat \Gamma_{c}|$','Location','Best');
 xlim([-1,1]*2)
 xlabel('$\omega$');
 grid on;


