function [] = PlotReductionCurve_slab(STD_curve,path,N_coils)
    
    currentnum = N_coils;
    [iteration, ~] = size(STD_curve);
    for i = 1:iteration
        xaxis(i) = currentnum;
        reducedper = round(currentnum/100)+1;
        currentnum = currentnum - reducedper;
    end
    f1 = figure('Name','Slab Shimming Results per Iteration');
    set(gcf,'position',[0,100,2000,400])
    hold on
    t = title('Slab Shimming Results(STD) per Iteration');
    t.FontSize = 15;
    t.FontWeight = 'bold';
    [~,slab_num] = size(STD_curve);
    for i = 1:slab_num
        tempSTD = mean(STD_curve(:,i),2);
        plot(xaxis,tempSTD,'-+','MarkerSize',2,'LineWidth',2),grid on; hold on
%         plot(xaxis,tempSTD,'-+','MarkerSize',2,'Color',[0.4660 0.6740 0.1880],'LineWidth',2),grid on; hold on
    end
    legend('Slab1','Slab2','Slab3','Slab4','Slab5','Slab6')
%     legend('Slab1','Slab2','Slab6')
    xlabel('Reserved Coil Number','FontSize',15,'FontWeight','bold')
    ylabel('B0shim/B0','FontSize',15,'FontWeight','bold')
    saveas(f1,fullfile(path,'6_SlabShimmingResultsperIteration.png'))
%     close all
end

