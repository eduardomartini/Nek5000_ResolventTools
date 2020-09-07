function [DATA,RAWDATA]=readData1(ts,folder,iy,ia,iz,actuatorDataFrom,computeTargetsCSD)
    
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
        data_temp = dlmread(sprintf('%s/adjdir%1.0f/ProjShapes_01_dir.dat',folder,iy(i)),'',1,0);
        for j=1:ny; estData_y{i,j} = data_temp(:,[1,1+iy(j)]);end
        for j=1:nz; estData_z{i,j} = data_temp(:,[1,1+iz(j)]);end
    end

    %reads raw actuation data
    if strcmp(actuatorDataFrom,'direct')
        for j=1:na
            data_temp  = dlmread(sprintf('%s/diradj%1.0f/ProjShapes_01_dir.dat',folder,ia(j)),'',1,0);
            for i=1:ny
                actData_y{i,j}  = data_temp(:,[1,1+iy(i)]);
            end
            for j=1:nz
                actData_z{i,j}  = data_temp(:,[1,1+iz(i)]);
            end
        end
    elseif  strcmp(actuatorDataFrom,'adjoint')
        for i=1:ny
            data_temp    = dlmread(sprintf('%s/adjdir%1.0f/ProjShapes_01_adj.dat',folder,iy(i)),'',1,0);
            for j=1:na
                actData_y{i,j} = data_temp(:,[1,1+ia(j)]);
            end
        end
        for i=1:nz
            data_temp    = dlmread(sprintf('%s/adjdir%1.0f/ProjShapes_01_adj.dat',folder,iz(i)),'',1,0);
            for j=1:na
                actData_z{i,j} = data_temp(:,[1,1+ia(j)]);
            end
        end
    end

    %reads raw target data, if CSDs are desired
    if computeTargetsCSD
        for i=1:nz   
            data_temp = dlmread(sprintf('%s/adjdir%1.0f/ProjShapes_01_dir.dat',folder,iz(i)),'',1,0);
            for j=1:nz   
                tarData_z{i,j} = data_temp(:,[1,1+iz(j)]);
            end
        end
    end

    %interpolate estimation data
    estdata_y_interp   =zeros(ny,ny,nt);
    estdata_z_interp   =zeros(nz,ny,nt);
    for i=1:ny
        for j=1:ny
            t_est = estData_y{i,j}(:,1);
            estdata_y_interp(j,i,:)=reshape(interp1(t_est,estData_y{i,j}(:,2),ts,'PCHIP',0).',[],1,nt); 
        end
        for j=1:nz
            t_est = estData_z{i,j}(:,1);
            estdata_z_interp(j,i,:)=reshape(interp1(t_est,estData_z{i,j}(:,2),ts,'PCHIP',0).',[],1,nt);     
        end
    end
    
    %interpolate actuator data    
    actdata_y_interp   =zeros(ny,na,nt);
    actdata_z_interp   =zeros(nz,na,nt);
    for j=1:na
        for i=1:nz
            t_act = actData_z{i,j}(:,1);
            actdata_z_interp (i,j,:)=reshape(interp1(t_act,actData_z{i,j}(:,2).',ts,'PCHIP',0)     ,[],1,nt); 
        end
        for i=1:ny
            t_act = actData_y{i,j}(:,1);
            actdata_y_interp (i,j,:)=reshape(interp1(t_act,actData_y{i,j}(:,2).',ts,'PCHIP',0)     ,[],1,nt); 
        end
    end
    
    %interpolate estimation targer
    tardata_z_interp =zeros(nz,nz,nt);
    if computeTargetsCSD
        for i=1:nz
            for j=i:nz
                t_est = tarData_z{i,j}(:,1);
                tardata_z_interp(j,i,:)=reshape(interp1(t_est,tarData_z{i,j}(:,2),ts,'PCHIP',0).',[],1,nt);     
            end
        end
    end
DATA.estY=estdata_y_interp;
DATA.estZ=estdata_z_interp;
DATA.actY=actdata_y_interp;
DATA.actZ=actdata_z_interp;
if computeTargetsCSD
    DATA.tarZ=tardata_z_interp;
end

RAWDATA.estY=estData_y;
RAWDATA.estZ=estData_z;
RAWDATA.actY=actData_y;
RAWDATA.actZ=actData_z;



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
