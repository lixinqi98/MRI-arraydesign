function [] = PlotReductionCurve(mean_curve,STD_curve,path,N_coils)
    f1 = figure('Name','Shimming Results(Mean) per Iteration');
    set(gcf,'position',[0,100,2000,400])
    hold on
    t = title('Shimming Results(Mean) per Iteration');
    t.FontSize = 15;
    t.FontWeight = 'bold';
    iteration = length(mean_curve);
    currentnum = N_coils;
    for i = 1:iteration
        xaxis(i) = currentnum;
        reducedper = round(currentnum/100)+1;
        currentnum = currentnum - reducedper;
    end
    p1 = plot(xaxis,mean_curve,'-+','MarkerSize',2,'Color',[0.4660 0.6740 0.1880],'LineWidth',2);grid on;
    xlabel('Reserved Coil Number','FontSize',15,'FontWeight','bold')
    ylabel('B0shimMean/B0Mean','FontSize',15,'FontWeight','bold')
    saveas(f1,fullfile(path,'ShimmingResults(Mean)perIteration.png'))

    f2 = figure('Name','Shimming Results(STD) Curve');
    set(gcf,'position',[0,100,2000,400])
    hold on
    t = title('Shimming Results(STD) Curve');
    t.FontSize = 15;
    t.FontWeight = 'bold';
    iteration = length(STD_curve);
    p2 = plot(xaxis,STD_curve,'-+','MarkerSize',2,'Color',[0.8500 0.3250 0.0980],'LineWidth',2);grid on;
    xlabel('Reserved Coil Number','FontSize',15,'FontWeight','bold')
    ylabel('B0shimSTD/B0STD','FontSize',15,'FontWeight','bold')
    saveas(f2,fullfile(path,'ShimmingResults(STD)perIteration.png'))

    f3 = figure('Name','Shimming Results per Iteration');
    set(gcf,'position',[0,100,2000,400])
    hold on
    t = title('Shimming Results(STD&Mean) per Iteration');
    t.FontSize = 15;
    t.FontWeight = 'bold';
    y1 = mean_curve;
    y2 = STD_curve;
    yyaxis left
    plot(xaxis,y1,'-+','MarkerSize',2,'Color',[0.4660 0.6740 0.1880],'LineWidth',2),grid on;
    yyaxis right
    plot(xaxis,y2,'-*','MarkerSize',2,'Color',[0.8500 0.3250 0.0980],'LineWidth',2),grid on;
    legend('Mean','STD')
    xlabel('Reserved Coil Number','FontSize',15,'FontWeight','bold')
    ylabel('B0shim/B0','FontSize',15,'FontWeight','bold')
    saveas(f3,fullfile(path,'ShimmingResultsperIteration.png'))
    close all
end

