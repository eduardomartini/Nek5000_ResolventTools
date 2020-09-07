%% Script to compute optimal control using Wiener-Hopf formalism, based on a DNS readings
% All frequency dependent matrices and vectors are stored with the
% frequency being the last index.

set(0,'defaulttextInterpreter'  ,'latex');
set(0,'defaultlegendInterpreter','latex'); 
set(0,'defaultAxesFontsize'     ,14) ;
set(0,'defaultLegendFontsize'   ,16) ;
set(0,'defaultLineLineWidth'    ,3 ) ;
set(0,'DefaultLineMarkerSize'   ,8 ) ;

leg={} ;
initiClock=tic();
%
computeGains=true;
% Set up some functions to help plotting
sortX = @(x) ifftshift( x-(x>x(end)/2)*x(end)); % Makes a [0,x] array into a [-x/2,x/2].
sortY = @(x) ifftshift( x,2);
pfun = @(x) reshape(x,size(x,1)*size(x,2),size(x,3)); % Function for plot Tensors
getFreqVec = @(t) [((0:length(t)-1) - (0:length(t)-1 > (length(t)-1)/2 )*length(t))/(t(end)-t(1)) ];
%% Defines indexes for sensors, actuators and targets (corresponding to shape index in the Nek code
iy = [ 2  ]; ny = length(iy);
ia = [ 32 ]; na = length(ia);
iz = [ 5  ]; nz = length(iz);

actuatorDataFrom='adjoint';% adjoint : target and sensor adjoint runs, direct : actuator direct runs
computeTargetsCSD = false;

folder = './'; %folder where the results are stored.

% Discretization used for control Law. Input Tmax and dt. 
Tmax = 800;
for dt   = [0.5,0.1,0.05,0.02,0.01];    
nt   = round(Tmax/dt); 
Tmax = nt*dt;
t    = (0:nt-1)*dt;
ts   = t - nt/2*dt;
df   = 1/t(end);
freq = getFreqVec (t);
freqS= fftshift(freq);

[DATA,~]= readData(ts,folder,iy,ia,iz,actuatorDataFrom,true);

Wh_facTime=tic();

HGs = getHGs(DATA,1e-4,1e-3,-1i,1e-6); 
Gamma = getControlKernels(HGs);
WHtime = toc(Wh_facTime);

figure(1);
    subplot(321)
        plot(HGs.t,pfun(HGs.Hminus),'-') ;
        hold on;
        ylabel('$H_-$')
        xlim([-1,1]*100)
        xlabel('t');
    subplot(322)
        plot(HGs.t,pfun (HGs.Hplus),'-')
        hold on;
        ylabel('$H_+$')
        xlim([-1,1]*100)
        xlabel('t');
    subplot(323)
        plot(HGs.t,(pfun(HGs.Gminus)),'-')
        hold on;
        ylabel('$G_-$');
        xlim([-1,1]*10)
        xlabel('t');
    subplot(324)
        plot(HGs.t,(pfun (HGs.Gplus)),'-')
        hold on;
        ylabel('${G}_+$');
        xlim([-1,1]*10)
        xlabel('t');
    subplot(313)        
        plot(Gamma.t,pfun(Gamma.c));
        hold on;
        ylabel('$\Gamma_c$');
        xlim([-1,10])
        xlabel('t');
        leg{end+1} =  sprintf('dt = %4.3f, time= %4.1fs',dt,WHtime) ;
        legend(leg);
        disp(WHtime)
end