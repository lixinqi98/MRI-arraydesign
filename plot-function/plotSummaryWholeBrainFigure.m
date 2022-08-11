% plot summary figure of whole brain
figure,
% x = ['Slab1', 'Slab2', 'Slab3', 'Slab4', 'Slab5', 'Slab6'];
x = 1:6;
% 99to32 v4
y1 = [0.799873889	0.759798537	0.684684599	0.652003367	0.632181211	0.680873966];
err1 = [0.032456526	0.056384219	0.061201197	0.116922819	0.135418769	0.079076094];
% 
% errorbar(x,y1,err1);

% 99to53 v4
y2 = [0.771702054	0.731952463	0.646318426	0.595106397	0.604954694	0.672258952];
err2 = [0.032580514	0.048938059	0.061034282	0.086287373	0.135282239	0.095473588];

% SH7
y3 = [0.782308647	0.774670349	0.647190874	0.605276975	0.579505269	0.660759795]; 
err3 = [0.078683939	0.054554622	0.073824606	0.108661382	0.121440077	0.130108294];

% SH6
y4 = [0.797176741	0.801911193	0.665587212	0.640535811	0.625303309	0.679551559];
err4 = [0.08742835	0.058349209	0.074635502	0.127396088	0.160070121	0.126997189];

% SH5
y5 = [0.83627985	0.828358119	0.690800545	0.678065686	0.657549905	0.673307971];
err5 = [0.092370171	0.059809371	0.074713655	0.118316657	0.181164374	0.125270701];

% SH4
y6 = [0.9121318	0.906393514	0.717191831	0.699897988	0.689757505	0.790995175];
err6 = [0.109581383	0.077262523	0.08699477	0.121860777	0.1697267	0.140983151];

% SH3
y7 = [0.951064997	0.979403102	0.785427867	0.748536123	0.74151899	0.85095227];
err7 = [0.117013858	0.085291075	0.089673813	0.102276712	0.133705499	0.099100783];

% % SH2
% y8 = [0.7932 0.8651 0.8145 0.7961 0.8118 0.8375];
% err8 = [0.1053 0.0655 0.1039 0.1203 0.1042 0.0759];

% iPRES 32
y9 = [0.920210697	0.910984501	0.780070703	0.795758842	0.896331423	0.870666477];
err9 = [0.060766913	0.10519652	0.070103435	0.074528874	0.130137448	0.109672771];

% y1 = y1./y8;
% y2 = y2./y8;
% y3 = y3./y8;
% y4 = y4./y8;
% y5 = y5./y8;
% y6 = y6./y8;
% y7 = y7./y8;
% y9 = y9./y8;
% y1 = y1./y6;
% y2 = y2./y6;
% y3 = y3./y6;
% y4 = y4./y6;
% y5 = y5./y6;
% 
% y9 = y9./y6;
% y6 = y6./y6;
% 
% err1 = err1./y8;
% err2 = err2./y8;
% err3 = err3./y8;
% err4 = err4./y8;
% err5 = err5./y8;
% err6 = err6./y8;
% err7 = err7./y8;
% err9 = err9./y8;

% e1 = errorbar(x, y1, err1, '-o', 'MarkerSize', 5, 'MarkerEdgeColor','red', 'MarkerFaceColor','red'); hold on
% 
% e7 = errorbar(x, y7, err7, '.-'); hold on % SH3
% e6 = errorbar(x+0.1, y6, err6, '.-'); hold on % SH4
% e5 = errorbar(x+0.2, y5, err5, '.-'); hold on % SH5
% e9 = errorbar(x+0.3, y9, err9, '.-'); hold on % iPRES 32
% e1 = errorbar(x+0.4, y1, err1, '.-'); hold on % 99to32
% e4 = errorbar(x+0.5, y4, err4, '.-'); hold on % SH6
% e3 = errorbar(x+0.6, y3, err3, '.-'); hold on % SH7
% e2 = errorbar(x+0.7, y2, err2, '.-'); hold on % 99to53

err = [0 0 0 0 0 0];
% e7 = plot(x, y7, '.-', 'LineWidth', 3, 'Color', [0.6 0.6 0.6]); hold on % SH3
% e6 = plot(x, y6, '.-', 'LineWidth', 3, 'Color', [0.5 0.5 0.5]); hold on % SH4
% e5 = plot(x, y5, '.-', 'LineWidth', 3, 'Color', [0.4 0.4 0.4]); hold on % SH5
% e9 = plot(x, y9, '--.', 'LineWidth', 1, 'Color', 'b'); hold on % iPRES 32
% e1 = plot(x, y1, '--.', 'LineWidth', 1, 'Color', 'g'); hold on % 99to32
% e4 = plot(x, y4, '.-', 'LineWidth', 3, 'Color', [0.3 0.3 0.3]); hold on % SH6
% e3 = plot(x, y3, '.-', 'LineWidth', 3, 'Color', [0.2 0.2 0.2]); hold on % SH7
% e2 = plot(x, y2, '--.', 'LineWidth', 1, 'Color', 'r'); hold on % 99to53
e7 = plot(x, y7, '.-', 'LineWidth', 3, 'Color', 'b'); hold on % SH3
e6 = plot(x, y6, '.-', 'LineWidth', 3); hold on % SH4
e5 = plot(x, y5, '.-', 'LineWidth', 3, 'Color', 'g'); hold on % SH5
e9 = plot(x, y9, '--.', 'LineWidth', 1, 'Color', 'b'); hold on % iPRES 32
e1 = plot(x, y1, '--.', 'LineWidth', 1, 'Color', 'g'); hold on % 99to32
e4 = plot(x, y4, '.-', 'LineWidth', 3); hold on % SH6
e3 = plot(x, y3, '.-', 'LineWidth', 3, 'Color', 'r'); hold on % SH7
e2 = plot(x, y2, '--.', 'LineWidth', 1, 'Color', 'r'); hold on % 99to53
xlim([0 7])
legend('SH3', 'SH4', 'SH5', 'iPRES32', 'UNIC 99to32', 'SH6', 'SH7', 'UNIC99to53', 'Location', 'northeastoutside')
% legend('SH3', 'SH4', 'SH5', 'SH6', 'SH7', 'Location', 'northeastoutside')
ylim([0.55 1])