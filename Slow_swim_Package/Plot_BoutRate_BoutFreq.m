function [] = Plot_BoutRate_BoutFreq(output)
%% subplot pre_Stimulus
Fish_WT=output(find([output.FishGeno]==2));

Fish_Homo=output(find([output.FishGeno]==0));


h1=figure(1); 

subplot(2,2,1)
title('BoutFrequency(Hz)');hold on;
%boxplot_2gp([Fish_WT.BoutFrequency],[Fish_Homo.BoutFrequency]);hold on;
boxplot_mean([Fish_WT.BoutFrequency],[Fish_Homo.BoutFrequency]);hold on;
%ylim([0.1 0.5]);
hold off;
subplot(2,2,2)
title('Bout Rate(Hz)');hold on;
%boxplot_2gp([Fish_WT.BoutRate],[Fish_Homo.BoutRate]);hold on;
boxplot_mean([Fish_WT.BoutRate],[Fish_Homo.BoutRate]);hold on;
%ylim([0.1 0.5]);
hold off;

%saveas(h1,['BoutFrequency and BoutRate.fig'])
%saveas(h1,['BoutFrequency and BoutRate.epsc'])
end

