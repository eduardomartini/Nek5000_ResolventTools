function [DATA]=readData2(ts,folder,iy,ia)
    % Read data for control of full rank forces and targets 
    % READS DATA FROM ADJ-DIR-ADJ iterations for sensors, and ADJ-DIR run for
    % targets.
    
    fprintf('Reading raw data...');clock=tic();
    ny = length(iy);
    na = length(ia);
    nt = length(ts);
    dt = ts(2)-ts(1);
    df = 1/(nt*dt);

     FFT = @(x) ifftshift(ifft( fftshift(    x   ,3),[],3)   ,3)*dt*nt; % fft defined for conventions used;
    IFFT = @(x) ifftshift( fft( fftshift(    x   ,3),[],3)   ,3)*df; % fft defined for conventions used;
    
    %reads raw estimation data
    for i=1:ny
        data_temp = dlmread(sprintf('%s/Run_%1.0f_adj-dir_full/ProjShapes_01.dat',folder,iy(i)),'',1,0);
        for j=1:ny; RfyRfyd    {j,i} = data_temp(:,[1,1+iy(j)]);end
        
        data_temp = dlmread(sprintf('%s/Run_%1.0f_adj-dir-adj_full/ProjShapes_01.dat',folder,iy(i)),'',1,0);
        for j=1:na; RazdRfzRfyd{j,i} = data_temp(:,[1,1+ia(j)]);end
    end

    %reads raw actuation data
    for i=1:na
        data_temp    = dlmread(sprintf('%s/Run_%1.0f_dir_full/ProjShapes_01.dat',folder,ia(i)),'',1,0);
        for j=1:ny
            Ray{j,i} = data_temp(:,[1,1+iy(j)]);
        end
        
        data_temp    = dlmread(sprintf('%s/Run_%1.0f_dir-adj_full/ProjShapes_01.dat',folder,ia(i)),'',1,0);
        for j=1:na
            RazdRaz{j,i} = data_temp(:,[1,1+ia(j)]);
        end
    end

    %interpolate estimation data
    RfyRfyd_interp    =zeros(ny,ny,nt);
    RazdRfzRfyd_interp=zeros(na,ny,nt);
    
    for i=1:ny
        for j=1:ny
            t_est = RfyRfyd{i,j}(:,1);
            RfyRfyd_interp(j,i,:)=reshape(interp1(t_est,RfyRfyd{i,j}(:,2),ts,'PCHIP',0).',[],1,nt); 
        end
        for j=1:na
            % Time for adjoint run is corrected.
            t_est = -RazdRfzRfyd{j,i}(:,1);
            RazdRfzRfyd_interp(j,i,:)=reshape(interp1(t_est,RazdRfzRfyd{j,i}(:,2),ts,'PCHIP',0).',[],1,nt);     
        end
    end
    
    %interpolate actuator data    
    Ray_interp   =zeros(ny,na,nt);
    RazdRaz_interp   =zeros(na,na,nt);
    for j=1:na
        for i=1:ny
            t_act = Ray{i,j}(:,1);
            Ray_interp (i,j,:)=reshape(interp1(t_act,Ray{i,j}(:,2).',ts,'PCHIP',0)     ,[],1,nt); 
        end
        for i=1:na
            t_act = -RazdRaz{i,j}(:,1);
            RazdRaz_interp (i,j,:)=reshape(interp1(t_act,RazdRaz{i,j}(:,2).',ts,'PCHIP',0)     ,[],1,nt); 
        end
    end
    
    
DATA.RfyRfyd    = RfyRfyd_interp;
DATA.RazdRaz    = RazdRaz_interp ;
DATA.Ray        = Ray_interp;
DATA.RazdRfzRfyd= RazdRfzRfyd_interp;


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

