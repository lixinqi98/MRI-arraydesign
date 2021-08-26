% use errorbar to plot the slab shimming results

figure,
% x = ['Slab1', 'Slab2', 'Slab3', 'Slab4', 'Slab5', 'Slab6'];
x = 1:6;
% 99to32 v4
y1 = [0.6550
0.6374
0.4430
0.4308
0.4697
0.5569];
% y1 = [0.6689 0.6143 0.5033 0.5078 0.5523 0.6350];
% err1 = [0.0717, 0.0551, 0.0620, 0.0774, 0.0734, 0.0910];
% % 
% % errorbar(x,y1,err1);
% 
% % 99to53 v4
y2 = [0.6041
0.5809
0.4136
0.3877
0.4403
0.5325];
% y2 = [0.6227 0.5552 0.4589 0.4593 0.5121 0.6073];
% err2 = [0.0737 0.0587 0.0627 0.0735 0.0825 0.0917];
% 
% % SH7
y3 = [0.6078
0.6142
0.4135
0.3752
0.4561
0.5363];
% y3 = [0.6283 0.5840 0.4650 0.4421 0.5273 0.6106]; 
% err3 = [0.0688 0.0687 0.0818 0.0751 0.0826 0.0914];
% 
% % SH6
y4 = [0.6307
0.6351
0.4248
0.3860
0.4698
0.5488];
% y4 = [0.6437 0.6107 0.4828 0.4595 0.5455 0.6221];
% err4 = [0.0697 0.0721 0.0883 0.0749  0.0850 0.0908];
% 
% % SH5
y5 = [0.6522
0.6852
0.4608
0.4187
0.4835
0.5718];
% y5 = [0.6677 0.6567 0.5344 0.5014 0.5654 0.6414];
% err5 = [0.0706 0.0748 0.0931 0.0802 0.0875 0.0938];
% 
% % SH4
y6 = [0.6921
0.7225
0.5259
0.4772
0.5252
0.6000];
% y6 = [0.7034 0.6990 0.6017 0.5725 0.6137 0.6683];
% err6 = [0.0805 0.0615 0.0940 0.0906 0.0913 0.0982];
% 
% % SH3
y7 = [0.7414
0.8533
0.6566
0.5774
0.5942
0.6426];
% y7 = [0.7509 0.8233 0.7158 0.6558 0.6701 0.7106];
% err7 = [0.0953 0.0741 0.1083 0.1012 0.0962 0.1088];
% 
% % SH2
% y8 = [0.7932 0.8651 0.8145 0.7961 0.8118 0.8375];
% err8 = [0.1053 0.0655 0.1039 0.1203 0.1042 0.0759];
% 
% % iPRES 32
y9 = [0.6575
0.6711
0.4939
0.5294
0.6173
0.6575];
% y9 = [0.6646 0.6539 0.5626 0.5628 0.6114 0.6715];
% err9 = [0.0724 0.0580 0.0806 0.0874 0.0805 0.0986];

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

err = [0 0 0 0 0 0];
% e7 = plot(x, y7, '.-', 'LineWidth', 3, 'Color', [0.6 0.6 0.6]); hold on % SH3
% e6 = plot(x, y6, '.-', 'LineWidth', 3, 'Color', [0.5 0.5 0.5]); hold on % SH4
% e5 = plot(x, y5, '.-', 'LineWidth', 3, 'Color', [0.4 0.4 0.4]); hold on % SH5
% e9 = plot(x, y9, '--.', 'LineWidth', 1, 'Color', 'b'); hold on % iPRES 32
% e1 = plot(x, y1, '--.', 'LineWidth', 1, 'Color', 'g'); hold on % 99to32
% e4 = plot(x, y4, '.-', 'LineWidth', 3, 'Color', [0.3 0.3 0.3]); hold on % SH6
% e3 = plot(x, y3, '.-', 'LineWidth', 3, 'Color', [0.2 0.2 0.2]); hold on % SH7
% e2 = plot(x, y2, '--.', 'LineWidth', 1, 'Color', 'r'); hold on % 99to53
e7 = plot(x, y7, '.-', 'LineWidth', 3); hold on % SH3
e6 = plot(x, y6, '.-', 'LineWidth', 3, 'Color', 'b'); hold on % SH4
e5 = plot(x, y5, '.-', 'LineWidth', 3, 'Color', 'g'); hold on % SH5
e9 = plot(x, y9, '--.', 'LineWidth', 1, 'Color', 'b'); hold on % iPRES 32
e1 = plot(x, y1, '--.', 'LineWidth', 1, 'Color', 'g'); hold on % 99to32
e4 = plot(x, y4, '.-', 'LineWidth', 3); hold on % SH6
e3 = plot(x, y3, '.-', 'LineWidth', 3, 'Color', 'r'); hold on % SH7
e2 = plot(x, y2, '--.', 'LineWidth', 1, 'Color', 'r'); hold on % 99to53
xlim([0 7])
legend('SH3', 'SH4', 'SH5', 'iPRES32', 'UNIC 99to32', 'SH6', 'SH7', 'UNIC99to53', 'Location', 'northeastoutside')
% legend('SH3', 'SH4', 'SH5', 'SH6', 'SH7', 'Location', 'northeastoutside')
ylim([0.35 0.88])