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
fSize=[1000,150];
aSize=[.075,.3,.9,.65];

fig_ker=figure; fig_ker.Position(3:4)=fSize;
fig_G  =figure; fig_G.Position(3:4)=fSize;
fig_g  =figure; fig_g.Position(3:4)=fSize;
% fig_deb=figure; fig_deb.Position(3:4)=[1400,400];
hold on; 

cmap = hot(nList*2);
% cmap = (hot(5*2));
cmap = flip(cmap(1:end/2,:),1);
leg={};


for i= 1:nList(end)
    leg{end+1}=sprintf('$\\Delta T=%3.1e$',iendList(i)*dt_data);

%     figure(fig_deb)
%         subplot(421)
%             plot(Gamma.t,(pfun(HGs_emp{i}.G)),'color',cmap(i,:)); 
%             hold on
%             xlim([-10,10])
%             xlabel('$\tau$');ylabel('$G,H$');
%         subplot(422)
%             plot(Gamma.t,(pfun(HGs_emp{i}.g)),'color',cmap(i,:)); 
%             hold on
%             xlim([-10,60])
%             xlabel('$\tau$');ylabel('$g$');
%         subplot(423)
%             plot(Gamma.t,(pfun(HGs_emp{i}.Gminus)),'color',cmap(i,:)); 
%             hold on
%             xlim([-10,10])
%             ylim([-1,5])
%             xlabel('$\tau$');ylabel('$G_-$');
%         subplot(424)
%             plot(Gamma.t,(pfun(HGs_emp{i}.Gplus)),'color',cmap(i,:)); 
%             hold on
%             xlim([-10,10])
%             ylim([-1,5])
%             xlabel('$\tau$');ylabel('$G_+$');
%         subplot(425)
%             plot(Gamma.t,(pfun(HGs_emp{i}.Hminus)),'color',cmap(i,:)); 
%             hold on
%             xlim([-10,10])
%             ylim([-1,5])
%             xlabel('$\tau$');ylabel('$H_-$');
%         subplot(426)
%             plot(Gamma.t,(pfun(HGs_emp{i}.Hplus)),'color',cmap(i,:)); 
%             hold on
%             xlim([-10,10])
%             ylim([-1,5])
%             xlabel('$\tau$');ylabel('$H_+$');
%         subplot(414)
%             plot(Gamma.t,(pfun(Gamma_emp{i}.c)),'color',cmap(i,:)); 
%             hold on
%             xlim([-10,50])
%             xlabel('$\tau$');ylabel('$\Gamma$');

        figure(fig_ker)
                plot(Gamma.t,(pfun(Gamma_emp{i}.c)),'color',cmap(i,:)); 
                hold on
                xlim([-10,20])
                xlabel('$\tau$');ylabel('$\Gamma$');
        figure(fig_G)
              plot(DATA.t,pfun(HGs_emp{i}.G),'color',cmap(i,:)); hold on;
                hold on
                xlim([-10,10])
                xlabel('$\tau$');ylabel('$G$');
        figure(fig_g)
              plot(DATA.t,pfun(HGs_emp{i}.g),'color',cmap(i,:)); hold on;
                hold on
                xlim([10,40])
                xlabel('$\tau$');ylabel('$g$');
        K_emp(:,i) = pfun(Gamma_emp{i}.c);
        
end

leg{end+1}='Analitical';
%    figure(fig_deb)
%     subplot(421)
%         plot(Gamma.t,(pfun(HGs.G)),'--b'); 
%     subplot(422)
%         plot(Gamma.t,(pfun(HGs.g)),'--b'); 
%     subplot(423)
%         plot(Gamma.t,(pfun(HGs.Gminus)),'--b'); 
%     subplot(424)
%         plot(Gamma.t,(pfun(HGs.Gplus)),'--b'); 
%     subplot(425)
%         plot(Gamma.t,(pfun(HGs.Hminus)),'--b'); 
%     subplot(426)
%         plot(Gamma.t,(pfun(HGs.Hplus)),'--b'); 
%     subplot(414)
%         plot(Gamma.t,(pfun(Gamma.c)),'--b'); 
   figure(fig_ker)
        plot(Gamma.t,(pfun(Gamma.c)),':b');
        legend(leg,'Location','east');
        ylim([-5.5,2.5]);
        a=gca;a.Position=aSize
        xlim([-.5,15]);
        ylim([-.5,1.5]);
        grid on;
        saveas(gcf,'Emp_Gamma_Conv.fig');
        saveas(gcf,'Emp_Gamma_Conv.eps','epsc');
   figure(fig_G)
        plot(Gamma.t,pfun(HGs.G),':b');
        legend(leg,'Location','east');
        a=gca;a.Position=aSize
        xlim([-10,10]);
        ylim([-.25,.75]);
        grid on;
        saveas(gcf,'Emp_Gbig_Conv.fig');
        saveas(gcf,'Emp_Gbig_Conv.eps','epsc');
%         ylim([-5.5,2.5]);
   figure(fig_g)
        plot(Gamma.t,pfun(HGs.g),':b');
        legend(leg,'Location','east');
        a=gca;a.Position=aSize;
        grid on;
        saveas(gcf,'Emp_gsmall_Conv.fig');
        saveas(gcf,'Emp_gsmall_Conv.eps','epsc');
        %         ylim([-5.5,2.5]);

% plot errors        

K = pfun(Gamma.c);
errK = K_emp-K.';
errKNorm=sqrt(diag(errK'*errK))/norm(K);

f=figure; f.Position(3:4)=fSize;
    T=iendList*dt_data-50;
    loglog(T,(errGNorm),'-o', ...
           T,(errgNorm),'-o', ...
           T,errKNorm,'-o', ...
           T,10*T.^-.5,'k:'       );
       legend('$G$','$g$','$\Gamma$','$T^{-1/2}$');
    xlim(iendList([1,end])*dt_data-50)
    xlabel('data length');
    ylabel('error norm ');
    a=gca;a.Position=aSize;
    grid
    saveas(gcf,'Emp_error_Conv.fig');
    saveas(gcf,'Emp_error_Conv.eps','epsc');

