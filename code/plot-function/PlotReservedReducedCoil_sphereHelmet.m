function [] = PlotReservedReducedCoil_sphereHelmet(resultprefix,label,mapped_coilcordf,x11f,y11f,z11f,rr)
    load(fullfile(resultprefix,'deletelist.mat'));
    resultprefix = fullfile(resultprefix,'ReservedCoilImg');
    mkdir(resultprefix)
    height = 2*rr;
    N_coils = size(mapped_coilcordf,3);
    f2 = figure('Name','Reduce Result');
    set(gcf,'position',[0,100,2000,550])
    for i = 1:length(deletelist)
        for view_index = 1:5
            for j = 1:N_coils
                % attention the transfer of the coordinate
                if(find(deletelist(1:i)==j))
                    subplot(2,5,view_index+5)
                    x1f = mapped_coilcordf(1,:,j);
                    y1f = mapped_coilcordf(2,:,j);
                    z1f = mapped_coilcordf(3,:,j);
                    xc1 = mean(x1f);
                    yc1 = mean(y1f);
                    zc1 = mean(z1f);
                    text(double(xc1),double(yc1),double(zc1),num2str(j),'FontSize',10,'color','w')
                    plot3(x1f,y1f,z1f,'linewidth',1.5), hold on
                    R = rr;   
                    [X,Y,Z] = cylinder(R,100);
                    Z = (Z - 1)*height/2+0.01;
                    test1 = surf(X,Y,Z,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
                    set(test1,'FaceColor',[0 0 0],'FaceAlpha',0.3)
                    hold on
                    R1=rr;
                    [X1,Y1,Z1] = sphere;
                    Z1(Z1<-0.001) = NaN;
                    X1 = X1*R1; Y1 = Y1*R1; Z1 = Z1*R1;
                    test11 = surf(X1,Y1,Z1,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
                    set(test11 ,'FaceColor',[0 0 0],'FaceAlpha',0.3)
                    hold on
                    xlabel x
                    ylabel y
                    zlabel z
                else
                    subplot(2,5,view_index)
                    x1f = mapped_coilcordf(1,:,j);
                    y1f = mapped_coilcordf(2,:,j);
                    z1f = mapped_coilcordf(3,:,j);
                    xc1 = mean(x1f);
                    yc1 = mean(y1f);
                    zc1 = mean(z1f);
                    text(double(xc1),double(yc1),double(zc1),num2str(j),'FontSize',10,'color','w')
                    plot3(x1f,y1f,z1f,'linewidth',1.5), hold on
                    R = rr;   
                    [X,Y,Z] = cylinder(R,100);
                    Z = (Z - 1)*height/2+0.01;
                    test1 = surf(X,Y,Z,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
                    set(test1,'FaceColor',[0 0 0],'FaceAlpha',0.3)
                    hold on
                    R1=rr;
                    [X1,Y1,Z1] = sphere;
                    Z1(Z1<-0.001) = NaN;
                    X1 = X1*R1; Y1 = Y1*R1; Z1 = Z1*R1;
                    test11 = surf(X1,Y1,Z1,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
                    set(test11 ,'FaceColor',[0 0 0],'FaceAlpha',0.3)
                    hold on
                    xlabel x
                    ylabel y
                    zlabel z
                end
            end
            for row = 1:2
                subplot(2,5,(row-1)*5+view_index)
                if(view_index == 1)
                    view(180,0);
                elseif(view_index == 2)
                    view(90,0);
                elseif(view_index == 3)
                    view(0,90);
                elseif(view_index == 4)
                    view(0,0);
                elseif(view_index == 5)
                    view(270,0);
                end
                plot3(x11f,y11f,z11f,'b.','MarkerSize',6),hold on; 
                hold on
            end
        end
        currentnum = N_coils - i;
        t = sgtitle(['Reduced Result, Reserved Coil Number ',num2str(currentnum)]);
        t.FontSize = 15;
        t.FontWeight = 'bold';
        savefig(f2,fullfile(resultprefix,[label,'_reservedCoil_',num2str(currentnum),'.fig']));
        saveas(f2,fullfile(resultprefix,[label,'_reservedCoil_',num2str(currentnum),'.png']));
        clf
    end
end

