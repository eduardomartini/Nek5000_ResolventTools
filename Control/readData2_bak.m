function [DATA,RAWDATA]=readData2(ts,folder,iy,ia,iz)
    % Read data for control of full rank forces and targets 
    % READS DATA FROM ADJ-DIR-ADJ iterations for sensors, and ADJ-DIR run for
    % targets.
    
    fprintf('Reading raw data...');clock=tic();
    ny = length(iy);
    na = length(ia);
    nz = length(iz);
    nt = length(ts);
    dt = ts(2)-ts(1);
    df = 1/(nt*dt);

     FFT = @(x) ifftshift(ifft( fftshift(    x   ,3),[],3)   ,3)*dt*nt; % fft defined for conventions used;
    IFFT = @(x) ifftshift( fft( fftshift(    x   ,3),[],3)   ,3)*df; % fft defined for conventions used;
    
    %reads raw estimation data
    for i=1:ny
        data_temp = dlmread(sprintf('%s/adj-dir-adj_%1.0f/ProjShapes_01_adj-dir.dat',folder,iy(i)),'',1,0);
        for j=1:ny; RfyRfyd    {i,j} = data_temp(:,[1,1+iy(j)]);end
        for j=1:nz; RfzRfyd    {i,j} = data_temp(:,[1,1+iz(j)]);end
        
        data_temp = dlmread(sprintf('%s/adj-dir-adj_%1.0f/ProjShapes_01_adj-dir-adj.dat',folder,iy(i)),'',1,0);
        for j=1:na; RazdRfzRfyd{i,j} = data_temp(:,[1,1+ia(j)]);end
    end
    
    for i=1:ny
    end

    %reads raw actuation data
    for i=1:ny
        data_temp    = dlmread(sprintf('%s/adj-dir-adj_%1.0f/ProjShapes_01_adj.dat',folder,iy(i)),'',1,0);
        for j=1:na
            Ray{i,j} = data_temp(:,[1,1+ia(j)]);
        end
    end
    
    for i=1:nz
        data_temp    = dlmread(sprintf('%s/adj-dir_%1.0f/ProjShapes_01_adj.dat',folder,iz(i)),'',1,0);
        for j=1:na
            %note adjoint run runs backward in time (which needs to be corrected), and provides  Razd.
            % By not correcting the time we directly obtain Rad, instead of
            % Radz.
            Raz{i,j} = data_temp(:,[1,1+ia(j)]);
        end
        data_temp    = dlmread(sprintf('%s/adj-dir_%1.0f/ProjShapes_01_adj-dir.dat',folder,iz(i)),'',1,0);
        for j=1:na
            RazRazd{i,j} = data_temp(:,[1,1+ia(j)]);
        end
    end

    %reads raw target data, if CSDs are desired
    for i=1:nz   
        data_temp = dlmread(sprintf('%s/adj-dir_%1.0f/ProjShapes_01_adj-dir.dat',folder,iz(i)),'',1,0);
        for j=1:nz   
            RfzRfzd{i,j} = data_temp(:,[1,1+iz(j)]);
        end
    end

    %interpolate estimation data
    RfyRfyd_interp    =zeros(ny,ny,nt);
    RfzRfyd_interp    =zeros(nz,ny,nt);
    RazdRfzRfyd_interp=zeros(na,ny,nt);
    
    for i=1:ny
        for j=1:ny
            t_est = RfyRfyd{i,j}(:,1);
            RfyRfyd_interp(j,i,:)=reshape(interp1(t_est,RfyRfyd{i,j}(:,2),ts,'PCHIP',0).',[],1,nt); 
        end
        for j=1:nz
            t_est = RfzRfyd{i,j}(:,1);
            RfzRfyd_interp(j,i,:)=reshape(interp1(t_est,RfzRfyd{i,j}(:,2),ts,'PCHIP',0).',[],1,nt);     
        end
        for j=1:na
            % Time for adjoint run is corrected.
            t_est = -RazdRfzRfyd{i,j}(:,1);
            RazdRfzRfyd_interp(j,i,:)=reshape(interp1(t_est,RazdRfzRfyd{i,j}(:,2),ts,'PCHIP',0).',[],1,nt);     
        end
    end
    
    %interpolate actuator data    
    Ray_interp   =zeros(ny,na,nt);
    Raz_interp   =zeros(nz,na,nt);
    for j=1:na
        for i=1:nz
            t_act = Raz{i,j}(:,1);
            Raz_interp (i,j,:)=reshape(interp1(t_act,Raz{i,j}(:,2).',ts,'PCHIP',0)     ,[],1,nt); 
        end
        for i=1:ny
            t_act = Ray{i,j}(:,1);
            Ray_interp (i,j,:)=reshape(interp1(t_act,Ray{i,j}(:,2).',ts,'PCHIP',0)     ,[],1,nt); 
        end
    end
    
    %interpolate estimation target
    RfzRfzd_interp =zeros(nz,nz,nt);
    for i=1:nz
        for j=i:nz
            t_est = RfzRfzd{i,j}(:,1);
            RfzRfzd_interp(j,i,:)=reshape(interp1(t_est,RfzRfzd{i,j}(:,2),ts,'PCHIP',0).',[],1,nt);     
        end
    end
    
DATA.RfyRfyd=RfyRfyd_interp;
DATA.RfzRfyd=RfzRfyd_interp;
DATA.Ray=Ray_interp;
DATA.Raz=Raz_interp;
DATA.RfzRfzd=RfzRfzd_interp;
DATA.RazdRfzRfyd=RazdRfzRfyd_interp;

RAWDATA.estY=RfyRfyd;
RAWDATA.estZ=RfzRfyd;
RAWDATA.actY=RfyRfyd;
RAWDATA.actZ=RfzRfyd;



fields=fieldnames(DATA);
for fieldId=1:length(fields)
    field = fields{fieldId};
    DATA.([field  'hat']) = FFT(DATA.(field));
end

DATA.t   =ts;
% getFreqVec = @(t) [((0:length(t)-1) - (0:length(t)-1 > (length(t)-1)/2 )*length(t))/(t(end)-t(1)) ];
getFreqVec = @(t) ((0:nt-1) - ceil(nt/2) )/(nt*dt);
DATA.freq= getFreqVec (ts);

disp([' Done in ' num2str(toc(clock)) 's']);

%%
% subplot(211)
% plot(ts,squeeze(DATA.estY(1,1,:)))
% subplot(212)
% plot(DATA.freq,squeeze(DATA.estYhat(1,1,:)),DATA.freq,squeeze(FFT(DATA.estY(1,1,:))),':r')
% 
% trapz(DATA.freq,squeeze(DATA.estYhat(1,1,:)))
