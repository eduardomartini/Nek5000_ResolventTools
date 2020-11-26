set(0,'defaultAxesFontsize'     ,16) ;
set(0,'defaultLegendFontsize'   ,14) ;
iy=2;
iz=5;

folder='./rnd1';
% data= dlmread(sprintf('%s/ProjShapes_01_rnd.dat',folder),'',1,0);
% t=data(:,1); y=data(:,iy+1);z=data(:,iz+1);
% save('data_small.mat','t','y','z');

% data= dlmread(sprintf('%s/ProjShapes_01.dat',folder),'',1,0);
load('data_large.mat');
% t=data(:,1); y=data(:,iy+1);z=data(:,iz+1);
% save('data_large.mat','t','y','z');

% load('data_small.mat');
% DATA =  load('data_large.mat');

% clear data;
% figure
%     plot(t(1:10000),y(1:10000));
%%
dt_data = t(2,1)-t(1,1);
nt_data = length(t);
%%
clear ryy1 ryz1 rzz1
% iendList = round(150/dt):round(100/dt):nt;
iendList = [ round(logspace(log10(round(100/dt_data)),log10(nt_data),5)) ];

for i=1:length(iendList)
    disp(i)
    istart  = round(50/dt_data);
    iend    = iendList(i);

    [ryy1(:,i),lagsyy]     = xcorr(y(istart:iend),y(istart:iend),round( 20/dt_data));
    [ryz1(:,i),lagsyz]     = xcorr(z(istart:iend),y(istart:iend),round(100/dt_data));
    [rzz1(:,i),lagszz]     = xcorr(z(istart:iend),z(istart:iend),round( 20/dt_data));

    %scale correlation w.r.t. rms values
    scale = max(DATA.RfyRfyd)/max(ryy1(:,i)); %.6373;%rms(data(end/5:end,iy+1)).^2;

    msy = 1 %.6373;%rms(data(end/5:end,iy+1)).^2;
    msz = 1 %19.333;%rms(data(end/5:end,iz+1)).^2;
    
    
    ryz1(:,i)= ryz1(:,i)*scale; 
    ryy1(:,i)= ryy1(:,i)*scale; 
    rzz1(:,i)= rzz1(:,i)*scale;
    
%     ryz1(:,i)= ryz1(:,i)*sqrt(msy*msz/(max(ryy1(:,i))*max(rzz1(:,i)))); 
%     ryy1(:,i)= ryy1(:,i)*    (msy    /(max(ryy1(:,i)))               ); 
%     rzz1(:,i)= rzz1(:,i)*    (    msz/               (max(rzz1(:,i)))); 
end
t_empy  = lagsyy*dt_data;
t_empz = lagsyz*dt_data;

%% smoothing
maskyy = (1-tanh( (abs(t_empy)-20)*1 ))/2;
maskyz = (1-tanh( (abs(t_empz-20)-40)*1 ))/2;
plot(t_empz,maskyz)

for i=1:length(iendList)
    ryy1(:,i) = maskyy'.*ryy1(:,i);
    ryz1(:,i) = maskyz'.*ryz1(:,i);
end

%% Target 

Gemp =  interp1(t_empy,ryy1,DATA.t,'cubic',0);
gemp =  interp1(t_empz,ryz1,DATA.t,'cubic',0);

G = pfun(DATA.RfyRfyd);
g = pfun(DATA.RfzRfyd);

errG = Gemp-G.';
errGNorm=sqrt(diag(errG'*errG))/norm(G);

errg = (gemp-g.');
errgNorm=sqrt(diag(errg'*errg))/norm(g);

%% Control Kernels comparison

dt=DATA.t(2)-DATA.t(1);
nt=length(DATA.t);
df=DATA.freq(2)-DATA.freq(1);

 FFT = @(x) ifftshift(ifft( fftshift(    x   ,3),[],3)   ,3)*dt*nt; % fft defined for conventions used;
IFFT = @(x) ifftshift( fft( fftshift(    x   ,3),[],3)   ,3)*df; % fft defined for conventions used;
DATA_emp={};
HGs_emp={};
Gamma_emp={};

nList = size(ryy1,2);


for i= 1:nList
    DATA_emp{i}=DATA;
    DATA_emp{i}.RfyRfyd(:) = Gemp(:,i);
    DATA_emp{i}.RfzRfyd(:) = gemp(:,i);
    
    DATA_emp{i}.RfyRfydhat(:) = FFT(DATA_emp{i}.RfyRfyd);
    DATA_emp{i}.RfzRfydhat(:) = FFT(DATA_emp{i}.RfzRfyd);

    HGs_emp{i} = getHGs(DATA_emp{i},P,N,-1i,1e-6); 
    Gamma_emp{i} = getControlKernels(HGs_emp{i});
end
HGs   = getHGs(DATA,P,N,-1i,1e-6); 
Gamma = getControlKernels(HGs);


%%
subplotloc=getSubPlotLocations([],[0.075,0.025,0.01,0.1],[0.01,0.01],1,4);
f=figure;
f.Position(3:4)=[1200,600];

cmap = hot(nList*2);
% cmap = (hot(5*2));
cmap = flip(cmap(1:end/2,:),1);
leg={};

fs1 = subplot('Position',subplotloc(1,:));
fs2 = subplot('Position',subplotloc(2,:));
fs3 = subplot('Position',subplotloc(3,:));
fs4 = subplot('Position',subplotloc(4,:));


for i= 1:nList(end)
    leg{end+1}=sprintf('$\\Delta T=%3.1e$',iendList(i)*dt_data);
        subplot(fs1);
                plot(Gamma.t,(pfun(Gamma_emp{i}.c)),'color',cmap(i,:)); 
                hold on
                xlim([-10,20])
                xlabel('$\tau$');ylabel('$\Gamma$');
        subplot(fs2);
              plot(DATA.t,pfun(HGs_emp{i}.G),'color',cmap(i,:)); hold on;
                hold on
                xlim([-10,10])
                xlabel('$\tau$');ylabel('$\mathbf G$');
        subplot(fs3);
              plot(DATA.t,pfun(HGs_emp{i}.g),'color',cmap(i,:)); hold on;
                hold on
                xlim([10,40])
                xlabel('$\tau$');ylabel('$\mathbf g$');
        K_emp(:,i) = pfun(Gamma_emp{i}.c);
        
end

leg{end+1}='Analitical';
    subplot(fs1);
        plot(Gamma.t,(pfun(Gamma.c)),':b');
        legend(leg,'Location','east');
        xlim([-.5,15]);
        ylim([-.5,1.5]);
        
        grid on;
    subplot(fs2);
        plot(Gamma.t,pfun(HGs.G),':b');
        xlim([-10,10]);
        ylim([-.29,.75]);
        grid on;
    subplot(fs3);
        plot(Gamma.t,pfun(HGs.g),':b');
        grid on;
        ylim([-1,1]*2.4);

% plot errors  

K = pfun(Gamma.c);
errK = K_emp-K.';
errKNorm=sqrt(diag(errK'*errK))/norm(K);

subplot(fs4);
    T=iendList*dt_data-50;
    loglog(T,(errGNorm),'-o', ...
           T,(errgNorm),'-o', ...
           T,errKNorm,'-o', ...
           T,10*T.^-.5,'k:'       );
       legend('$\mathbf G$','$\mathbf g$','$\Gamma$','$T^{-1/2}$');
    xlim(iendList([1,end])*dt_data-50)
    xlabel('data length');
    ylabel('error norm ');
    grid
    
fs1.XTickLabel="";
fs2.XTickLabel="";
fs3.XTickLabel="";
fs1.FontSize=18;
fs2.FontSize=18;
fs3.FontSize=18;
fs4.FontSize=18;
fs4.YTick= 10.^[-3:0];
    
saveas(gcf,'Emp_Conv.fig');
saveas(gcf,'Emp_Conv.eps','epsc');

