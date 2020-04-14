function [FWD_Fish_WT,FWD_Fish_Homo,RT_Fish_WT,RT_Fish_Homo] = Plot_FWD_RT(FWD_output,RT_output)

FWD_Fish_WT=FWD_output(find([FWD_output.FishGeno]==2));
FWD_Fish_Homo=FWD_output(find([FWD_output.FishGeno]==0));

RT_Fish_WT=RT_output(find([RT_output.FishGeno]==2));
RT_Fish_Homo=RT_output(find([RT_output.FishGeno]==0));
%% Forward Swim

subplot(2,6,1)
title('Bout Distance (mm)');hold on;
boxplot_mean([FWD_Fish_WT.BoutDistance],[FWD_Fish_Homo.BoutDistance]);hold on;
hold off;

subplot(2,6,2)
title('Bout Duration (sec)');hold on; 
boxplot_mean([FWD_Fish_WT.BoutDuration],[FWD_Fish_Homo.BoutDuration]);hold on;
hold off;

subplot(2,6,3)
title('Bout Speed (mm/sec)');hold on;
boxplot_mean([FWD_Fish_WT.Speed],[FWD_Fish_Homo.Speed]);hold on;
%ylim([0 5]);
%ylim([0 8]);
hold off;

subplot(2,6,4)
title('# Of Oscillations'); hold on;
boxplot_mean([FWD_Fish_WT.NumberOfOscillations],[FWD_Fish_Homo.NumberOfOscillations]);hold on;
%ylim([-2 5]);
hold off;

subplot(2,6,5)
title('TBF');hold on;  
boxplot_mean([FWD_Fish_WT.TBF],[FWD_Fish_Homo.TBF]);hold on;
hold off;

subplot(2,6,6)
title('Median Bend Amplitude (degree)');hold on;  
boxplot_mean([FWD_Fish_WT.MedianBendAmp],[FWD_Fish_Homo.MedianBendAmp]);hold on;
ylim([1.5 5.5]);
hold off;

% Routine Turn

subplot(2,6,7)
title('Bout Distance (mm)');hold on;
boxplot_mean([RT_Fish_WT.BoutDistance],[RT_Fish_Homo.BoutDistance]);hold on;
hold off;

subplot(2,6,8)
title('Bout Duration (sec)');hold on; 
boxplot_mean([RT_Fish_WT.BoutDuration],[RT_Fish_Homo.BoutDuration]);hold on;
%ylim([0 2]);
hold off;

subplot(2,6,9)
title('Bout Speed (mm/sec)');hold on;
boxplot_mean([RT_Fish_WT.Speed],[RT_Fish_Homo.Speed]);hold on;
%ylim([0 5]);
%ylim([0 8]);
hold off;

subplot(2,6,10)
title('# Of Oscillations'); hold on;
boxplot_mean([RT_Fish_WT.NumberOfOscillations],[RT_Fish_Homo.NumberOfOscillations]);hold on;
%ylim([-2 5]);
hold off;

subplot(2,6,11)
title('TBF');hold on;  
boxplot_mean([RT_Fish_WT.TBF],[RT_Fish_Homo.TBF]);hold on;
hold off;

subplot(2,6,12)
title('Median Bend Amplitude (degree)');hold on;  
boxplot_mean([RT_Fish_WT.MedianBendAmp],[RT_Fish_Homo.MedianBendAmp]);hold on;
%ylim([1.5 5.5]);
hold off;

% saveas(h2,'5min_GoodSwimmers_All_thres25_6Parameters.fig')
% saveas(h2,'5min_GoodSwimmers_All_thres25_6Parameters.epsc')
end

