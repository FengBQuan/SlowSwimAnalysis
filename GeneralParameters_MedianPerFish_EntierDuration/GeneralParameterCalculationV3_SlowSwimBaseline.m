%% Clean workspace and load data
close all
clear variables
clc
load("DataSetPreEscape_FWDvsRT_thres25.mat")
%% create output stucture
    
output = struct( 'NTrial', [], 'Fish_ID', [], 'FishGeno',[], ...
        'BoutRate', [], 'BoutDuration',[],'NumberOfOscillations',[],'BoutDistance',[],'Speed',[],'TBF',[],'MedianBendAmp',[],'MaxBendAmp',[],'BoutFrequency',[]);


%% define dataset

datasetPerBout=datasetPerBout_FWD1;
%datasetPerBout=DatasetPreEscape_FWD1;
%datasetPerBout=DatasetPreEscape_RT1;
%% find fish number and genotypes
Fish = unique([datasetPerFish(:).Condition]);
NumberFish=length(Fish);

FishGeno=([datasetPerFish.Genotype]);
Fish_ID = ([datasetPerFish.Condition]);
NTrial = ([datasetPerFish.NTrial]); % Ntiral= N of Clutch

%% calculate for each fish parameters extracted from the structure

for i=1:NumberFish;

    %Calculate IBI and BoutFrequency=1/IBI by using index without 1st one
    
    index{Fish(i)}= find(~([datasetPerBout(:).Condition]-Fish(i)));
   
    if numel(index{Fish(i)})>0;
        index{Fish(i)}(1)=[];
    else if numel(index{Fish(i)})==0;
            index{Fish(i)}=[];
        end
    end
    
    IBI{Fish(i)} = [datasetPerBout(index{Fish(i)}).InstantaneousIBI];
    medianIBI{Fish(i)}=median(IBI{Fish(i)},'omitnan');
    medianBoutFrequency{Fish(i)}=1/medianIBI{Fish(i)};
    
    % Calculate other parameters by using all index
    allindex{Fish(i)}= find(~([datasetPerBout(:).Condition]-Fish(i)));
     
    % if no tracking data, allindex =[], parameters=nan.
    if length(allindex{Fish(i)})==0;
        BoutDuration{Fish(i)}=nan;
        BoutDistance{Fish(i)}=nan;
        Speed{Fish(i)}=nan;
        NumberOfOscillations{Fish(i)}=nan;
        TBF{Fish(i)}=nan;
        MaxBendAmp{Fish(i)}=nan;
        MedianBendAmp{Fish(i)}=nan;
        BoutRate{Fish(i)}=nan;
    else
        
        for h=1:length(allindex{Fish(i)});
            
            display(['currently processing fish ' num2str(i)])
            display(['currently processing bout number ' num2str(h)])
            
            Bend_Amplitude{Fish(i)}{h} = 57.2958*datasetPerBout(allindex{Fish(i)}(h)).Bend_Amplitude;
            MaxBendAmp{Fish(i)}(h)=max(abs(57.2958*datasetPerBout(allindex{Fish(i)}(h)).Bend_Amplitude));
            MedianBendAmp{Fish(i)}(h)=median(abs(57.2958*datasetPerBout(allindex{Fish(i)}(h)).Bend_Amplitude));
            
        end
 
    %calculate parameters for all swims 
    BoutDuration{Fish(i)}=[datasetPerBout(allindex{Fish(i)}).BoutDuration];
    BoutDistance{Fish(i)}=[datasetPerBout(allindex{Fish(i)}).TotalDistance];
    Speed{Fish(i)}=[datasetPerBout(allindex{Fish(i)}).Speed];
    NumberOfOscillations{Fish(i)}=[datasetPerBout(allindex{Fish(i)}).NumberOfOscillations];
    TBF{Fish(i)}=[datasetPerBout(allindex{Fish(i)}).NumberOfOscillations]/[datasetPerBout(allindex{Fish(i)}).BoutDuration];
    BoutRate{Fish(i)}=numel([datasetPerBout(allindex{Fish(i)}).BoutDuration])/300; %300 sec =5min total recording duration
    
    end
    
    
    % outputTable
    output(i).NTrial=NTrial(i); 
    output(i).Fish_ID=Fish_ID(i);
    output(i).FishGeno=FishGeno(i); 
    
    output(i).BoutRate=BoutRate{Fish(i)};% nBouts/DurationTotale
    output(i).BoutDuration=median(BoutDuration{Fish(i)},'omitnan');
    output(i).BoutDistance=median(BoutDistance{Fish(i)},'omitnan');
    output(i).Speed=median(Speed{Fish(i)},'omitnan');
    output(i).NumberOfOscillations=mean(NumberOfOscillations{Fish(i)},'omitnan');
    output(i).TBF=median(TBF{Fish(i)},'omitnan');
    output(i).MedianBendAmp=median(MedianBendAmp{Fish(i)},'omitnan');
    output(i).MaxBendAmp=median(MaxBendAmp{Fish(i)},'omitnan');
    output(i).BoutFrequency=medianBoutFrequency{Fish(i)};%1/IBI
    
    i= i+1;

end;

%% selecte fish who swims more than 20 bouts

output=output(find([output.BoutRate]>(20/300))); 


%% Generation dataset for Fish_WT, Fish_Homo

Fish_WT=output(find([output.FishGeno]==2));

Fish_Homo=output(find([output.FishGeno]==0));

%% subplot pre_Stimulus
h1=figure(1); 

subplot(2,4,1)
title('BoutFrequency(Hz)');hold on;
boxplot_2gp([Fish_WT.BoutFrequency],[Fish_Homo.BoutFrequency]);hold on;
MeanErrorbar([Fish_WT.BoutFrequency],[Fish_Homo.BoutFrequency]);hold on;
%ylim([0.1 0.5]);
hold off;

subplot(2,4,2)
title('Bout Distance (mm)');hold on;
boxplot_2gp([Fish_WT.BoutDistance],[Fish_Homo.BoutDistance]);hold on;
MeanErrorbar([Fish_WT.BoutDistance],[Fish_Homo.BoutDistance]);hold on;
%ylim([0 7]);
%ylim([0 10]);
hold off;

subplot(2,4,3)
title('Bout Duration (sec)');hold on; 
boxplot_2gp([Fish_WT.BoutDuration],[Fish_Homo.BoutDuration]);hold on;
MeanErrorbar([Fish_WT.BoutDuration],[Fish_Homo.BoutDuration]);hold on;
%ylim([0 2]);
hold off;

subplot(2,4,4)
title('Bout Speed (mm/sec)');hold on;
boxplot_2gp([Fish_WT.Speed],[Fish_Homo.Speed]);hold on;
MeanErrorbar([Fish_WT.Speed],[Fish_Homo.Speed]);hold on;
%ylim([0 5]);
%ylim([0 8]);
hold off;

subplot(2,4,5)
title('Bout Rate(Hz)');hold on;
boxplot_2gp([Fish_WT.BoutRate],[Fish_Homo.BoutRate]);hold on;
MeanErrorbar([Fish_WT.BoutRate],[Fish_Homo.BoutRate]);hold on;
%ylim([0.1 0.5]);
hold off;

subplot(2,4,6)
title('# Of Oscillations'); hold on;
boxplot_2gp([Fish_WT.NumberOfOscillations],[Fish_Homo.NumberOfOscillations]);hold on;
MeanErrorbar([Fish_WT.NumberOfOscillations],[Fish_Homo.NumberOfOscillations]);hold on;
%ylim([-2 5]);
hold off;

subplot(2,4,7)
title('TBF');hold on;  
boxplot_2gp([Fish_WT.TBF],[Fish_Homo.TBF]);hold on;
MeanErrorbar([Fish_WT.TBF],[Fish_Homo.TBF]);hold on;
hold off;

subplot(2,4,8)
title('Median Bend Amplitude (degree)');hold on;  
boxplot_2gp([Fish_WT.MedianBendAmp],[Fish_Homo.MedianBendAmp]);hold on;
MeanErrorbar([Fish_WT.MedianBendAmp],[Fish_Homo.MedianBendAmp]);hold on;
ylim([1 6]);
hold off;

% subplot(2,3,6)
% title('Max Bend Amplitude (degree)');hold on;  
% boxplot_2gp([Fish_WT.MaxBendAmp],[Fish_Homo.MaxBendAmp]);hold on;
% MeanErrorbar([Fish_WT.MaxBendAmp],[Fish_Homo.MaxBendAmp]);hold on;
% %ylim([0 40]);
% hold off;


%saveas(h1,['5min_GoodSwimmers_FWD_thres25.fig'])
%saveas(h2,['5min_AllSwimmers_RT.epsc'])

%% Stat

[h,p,ci,stats] = ttest2([Fish_WT.BoutFrequency],[Fish_Homo.BoutFrequency])
[h,p,ci,stats] = ttest2([Fish_WT.BoutDistance],[Fish_Homo.BoutDistance])
[h,p,ci,stats] = ttest2([Fish_WT.BoutDuration],[Fish_Homo.BoutDuration])
[h,p,ci,stats] = ttest2([Fish_WT.Speed],[Fish_Homo.Speed])
[h,p,ci,stats] = ttest2([Fish_WT.NumberOfOscillations],[Fish_Homo.NumberOfOscillations])
[h,p,ci,stats] = ttest2([Fish_WT.TBF],[Fish_Homo.TBF])
[h,p,ci,stats] = ttest2([Fish_WT.MedianBendAmp],[Fish_Homo.MedianBendAmp])
[h,p,ci,stats] = ttest2([Fish_WT.BoutRate],[Fish_Homo.BoutRate])
%% Generation of text Table file for statistic analysis by Francois-Xavier

% Table=struct2table(output);
% writetable(Table);
