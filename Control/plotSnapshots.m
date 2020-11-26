clear time
close all
iplot = 80*2;
file={};
% file={'bf_bfsNL0.f00001'};
file ={     {'UncontrolRun_lincon'   , [],[],[]}  
            {'Run_lincon_ControlKernel_iy_2_ia_30_iz_5' ,[-.5,.5],[.5,.5],[10.5,0]} 
            {'Run_lincon_ControlKernel_iy_2_ia_30_iz_' ,[-.5,.5],[.5,.5],[]} 
            {'Run_lincon_ControlKernel_iy_1_2_3_ia_30_iz_5' ,[-.5,.5;-.5,.75;-.5,.25],[.5,.5],[10.5,0]} 
            {'Run_lincon_ControlKernel_iy_1_2_3_ia_30_iz_' ,[-.5,.5;-.5,.75;-.5,.25],[.5,.5],[]} 
            {'Run_lincon_ControlKernel_iy_1_2_3_ia_29_30_31_iz_4_5_6' ,[-.5,.5;-.5,.75;-.5,.25],[.5,.5;.5,.25;.5,.75],[10.5,0.5;10.5,0;10.5,-.5]} 
            {'Run_lincon_ControlKernel_iy_1_2_3_ia_29_30_31_iz_' ,[-.5,.5;-.5,.75;-.5,.25],[.5,.5;.5,.25;.5,.75],[]} 
           };
        
for ifolder=1:length(file)        
    file{ifolder}{5} = sprintf('%s/i01bfs0.f%05.0f',file{ifolder}{1},iplot);
end
 %%

%     file = {"kalmanGain0.f00001","kalmanGain0.f00002","kalmanGain0.f00003"};


mesh=readnek('bf_bfsNL0.f00001');

subplotloc=getSubPlotLocations([],[0.05,0.025,0.01,0.1],[0.01,0.01],2,ceil(length(file)/2))
f=figure;
f.Position(3:4)=[700,500];
uRange = [-1,1]*30;

for i=1:length(file)
    subplot('Position',subplotloc(i+(i>1),:));
    [data,~,~,time(i),~,~,~,~,~] = readnek(file{i}{5});
    [x,y,u]  = PlotNek(mesh(:,:,1:2),data(:,:,1),[],[],100,50);

    ncolors=128;
    clevels=linspace(uRange(1),uRange(2),ncolors);
    colormap(CustomColorMap(ncolors,[0,0,1],[1,0,0],[1,1,1]));
    contourf(x,y,u,clevels,'linecolor','none');
    caxis(uRange);
    hold on;
    for yy = file{i}{2}' 
%         t = text(yy(1,:),yy(2,:),'$y$');
%         t.FontSize=16;
        plot(yy(1,:),yy(2,:),'ko');
    end
    for yy = file{i}{3}'
%         t = text(yy(1,:),yy(2,:),'$a$');
%         t.FontSize=16;
        plot(yy(1,:),yy(2,:),'kx');
    end
    for yy = file{i}{4}'
%         t = text(yy(1,:),yy(2,:),'$z$');
%         t.FontSize=16;
        plot(yy(1,:),yy(2,:),'k^');
    end
    
    rectangle('Position',[-10,-1,10,1],'facecolor','w');

    t = text(-10,-1,['$(' 'a'+i-1 ')$']);
%     t.VerticalAlignment='top';
%     t.HorizontalAlignment='right';
    t.VerticalAlignment='bottom';
    t.HorizontalAlignment='left';
    t.FontSize=20;

    if mod(i+(i>1),2)==0
        a=gca;
        a.YTickLabel='';
    else
      ylabel('y')
    end
    
    if i<length(file)-1
        a=gca;
        a.XTickLabel='';
        a.YTick=[0,1];
    else
        xlabel('x')
    end
    
end
%%   
f = figure;
f.Position(3:4) = [1000,300];
aSize=[.075,.15,.9,.80];

leg={};
for ifolder=1:length(file)        
    data = dlmread([file{ifolder}{1} '/pertNormXp.txt']);
    plot(data(:,1),data(:,2));hold on;
    leg{ifolder}=['$(' 'a'+ifolder-1 ')$'];
end
xlabel('time');
ylabel('perturbation norm');
grid on;
xlim([400,600]);
ylim([0,max(ylim)]);
a=gca;a.Position=aSize;
        
legend(leg);