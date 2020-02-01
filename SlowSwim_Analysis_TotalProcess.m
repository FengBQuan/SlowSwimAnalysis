% Clean workspace and load data
close all
clear variables
clc
% load("matlab_workspace_SST_4manip_SlowSwim.mat")
% 
% %clean dataset
% save('SST_4manip_SlowSwim.mat', 'datasetPerBout', 'datasetPerFish')
% clear
load('SST_4manip_SlowSwim.mat');
%% Extract DataSet PreEscape(5min)

DatasetPreEscape= ExtractPreEscape(datasetPerBout);

%% Selecte good swimmers

[DatasetPreEscape_GoodSwimmers,GoodSwimmers]= SelecteGoodSwimmers(DatasetPreEscape,datasetPerFish);

%% CutOff FWD and RT dataset


[DatasetPreEscape_GoodSwimmers_FWD,DatasetPreEscape_GoodSwimmers_RT]  = CutOff_FWD_RT(DatasetPreEscape_GoodSwimmers);
  
% Plot Bout Proportion FWD vs RT 
%(WT_FWD,WT_RT,Homo_FWD,Homo_RT are percentage)
h1=figure(1);
[WT_FWD,WT_RT,Homo_FWD,Homo_RT]=ProportionPie(DatasetPreEscape_GoodSwimmers_FWD,DatasetPreEscape_GoodSwimmers_RT);

%% Bout Frequency (1/IBI), BoutRate (nBout/TotalDuration);

output=ParametersCalculations(DatasetPreEscape_GoodSwimmers);
output = rmfield(output, {'BoutDuration', 'NumberOfOscillations', 'BoutDistance', 'Speed', 'TBF', 'MedianBendAmp', 'MaxBendAmp'});

% Generation of text Table file for statistic analysis by Francois-Xavier
Table_AllBehavior=struct2table(output);
writetable(Table_AllBehavior);

% plot BoutRate and BoutFrequency for all Behavior
h2=figure(2); 
Plot_BoutRate_BoutFreq(output);

%% Forward Swim: BoutDistance, BoutDuration, Speed, NumOfOsc,TBF,MaxBend,MedianBend.
clear output

output=ParametersCalculations(DatasetPreEscape_GoodSwimmers_FWD);
output = rmfield(output, {'BoutRate', 'BoutFrequency'});
FWD_output=output;

% Generation of text Table file for statistic analysis by Francois-Xavier
Table_FWD=struct2table(FWD_output);
writetable(Table_FWD);

%% Routine Turn: BoutDistance, BoutDuration, Speed, NumOfOsc,TBF,MaxBend,MedianBend.
clear output

%datasetPerBout=datasetPerBout_RT1;
output=ParametersCalculations(DatasetPreEscape_GoodSwimmers_RT);
output = rmfield(output, {'BoutRate', 'BoutFrequency'});
RT_output=output;

% Generation of text Table file for statistic analysis 
Table_RT=struct2table(RT_output);
writetable(Table_RT);

clear output
%% Plot FWD vs RT

h3=figure(3); 
[FWD_Fish_WT,FWD_Fish_Homo,RT_Fish_WT,RT_Fish_Homo] =Plot_FWD_RT(FWD_output,RT_output);
% saveas(h2,'5min_GoodSwimmers_All_thres25_6Parameters.fig')
% saveas(h2,'5min_GoodSwimmers_All_thres25_6Parameters.epsc')

save('SST_4manip_SlowSwim_TotalProcess_thre25_moreThan0bouts_wo204.mat');
%% Stat preview

% % FWD Stat preview
[h,p,ci,stats] = ttest2([FWD_Fish_WT.BoutDistance],[FWD_Fish_Homo.BoutDistance])
[h,p,ci,stats] = ttest2([FWD_Fish_WT.BoutDuration],[FWD_Fish_Homo.BoutDuration])
[h,p,ci,stats] = ttest2([FWD_Fish_WT.Speed],[FWD_Fish_Homo.Speed])
[h,p,ci,stats] = ttest2([FWD_Fish_WT.NumberOfOscillations],[FWD_Fish_Homo.NumberOfOscillations])
[h,p,ci,stats] = ttest2([FWD_Fish_WT.TBF],[FWD_Fish_Homo.TBF])
[h,p,ci,stats] = ttest2([FWD_Fish_WT.MedianBendAmp],[FWD_Fish_Homo.MedianBendAmp])

% 
% % RT Stat preview
% 
[h,p,ci,stats] = ttest2([RT_Fish_WT.BoutDistance],[RT_Fish_Homo.BoutDistance])
[h,p,ci,stats] = ttest2([RT_Fish_WT.BoutDuration],[RT_Fish_Homo.BoutDuration])
[h,p,ci,stats] = ttest2([RT_Fish_WT.Speed],[RT_Fish_Homo.Speed])
[h,p,ci,stats] = ttest2([RT_Fish_WT.NumberOfOscillations],[RT_Fish_Homo.NumberOfOscillations])
[h,p,ci,stats] = ttest2([RT_Fish_WT.TBF],[RT_Fish_Homo.TBF])
[h,p,ci,stats] = ttest2([RT_Fish_WT.MedianBendAmp],[RT_Fish_Homo.MedianBendAmp])



