
    iplot = 150;
    file={};
    file{1} = sprintf('rnd_old/r01bfs0.f%05.0f',iplot);
    folder='control_2__5__38';
    file{end+1} = sprintf('%s/c01bfs0.f%05.0f',folder,iplot);
    folder='control_1_2_3__5__38';
    file{end+1} = sprintf('%s/c01bfs0.f%05.0f',folder,iplot);
    folder='control_1_2_3__4_5_6__38';
    file{end+1} = sprintf('%s/c01bfs0.f%05.0f',folder,iplot);
    folder='control_1_2_3__4_5_6__37_38_39';
    file{end+1} = sprintf('%s/c01bfs0.f%05.0f',folder,iplot/6);


%     file = {"kalmanGain0.f00001","kalmanGain0.f00002","kalmanGain0.f00003"};


    mesh=readnek('bf_bfsNL0.f00001');

    subplotloc=getSubPlotLocations([],[0.1,0.05,0.02,0.075],[0.01,0.01],1,5)
    f=figure;
    f.Position(3:4)=[600,700];
    for i=1:length(file)
        subplot('Position',subplotloc(i,:));
        data = readnek(file{i});
        [x,y,u] = PlotNek(mesh(:,:,1:2),data(:,:,1),[],[],100,50);

        ncolors=64;
        clevels=linspace(-1,1,ncolors)*0.2;%max(abs(u(:)));
        colormap(CustomColorMap(ncolors,[0,0,1],[1,0,0],[1,1,1]));
        contourf(x,y,u,clevels,'linecolor','none');
        caxis([-1,1]*0.2);
        hold on;
        if i==2      
            t = text(-0.5,0.5,'$y$');
            t.FontSize=16;
            t = text(10.5,0.0,'$z$');
            t.FontSize=16;
            t = text(1.5,0.5,'$a$');
            t.FontSize=16;
        elseif i==3 ; 
            t = text(-0.5,0.25,'$y$');
            t.FontSize=16;
            t = text(-0.5,0.75,'$y$');
            t.FontSize=16;
            t = text(-0.5,0.5,'$y$');
            t.FontSize=16;
            t = text(10.5,0.0,'$z$');
            t.FontSize=16;
            t = text(1.5,0.5,'$a$');
            t.FontSize=16;
        elseif i==4 ; 
            t = text(-0.5,0.25,'$y$');
            t.FontSize=16;
            t = text(-0.5,0.75,'$y$');
            t.FontSize=16;
            t = text(-0.5,0.5,'$y$');
            t.FontSize=16;
            t = text(1.5,0.5,'$a$');
            t.FontSize=16;
            t = text(10.5,-0.5,'$z$');
            t.FontSize=16;
            t = text(10.5,0.0,'$z$');
            t.FontSize=16;
            t = text(10.5,0.5,'$z$');
            t.FontSize=16;
        elseif i==5 ; 
            t = text(-0.5,0.25,'$y$');
            t.FontSize=16;
            t = text(-0.5,0.75,'$y$');
            t.FontSize=16;
            t = text(-0.5,0.5,'$y$');
            t.FontSize=16;
            t = text(1.5,0.25,'$a$');
            t.FontSize=16;
            t = text(1.5,0.5,'$a$');
            t.FontSize=16;
            t = text(1.5,0.75,'$a$');
            t.FontSize=16;
            t = text(10.5,-0.5,'$z$');
            t.FontSize=16;
            t = text(10.5,0.0,'$z$');
            t.FontSize=16;
            t = text(10.5,0.5,'$z$');
            t.FontSize=16;
        end
        if i<5;
            a=gca;
            a.XTickLabel='';
            a.YTick=[0,1];
        end
        rectangle('Position',[-10,-1,10,1],'facecolor','w');
        ylabel('y')
    end
    xlabel('x')