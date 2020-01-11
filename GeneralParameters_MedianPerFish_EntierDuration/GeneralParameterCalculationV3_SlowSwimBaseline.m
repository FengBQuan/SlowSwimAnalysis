%% Clean workspace and load data
close all
clear variables
clc
load("DataSetPreEscape_FWDvsRT_thres30.mat")
%% create output stucture
    
output = struct( 'NTrial', [], 'Fish_ID', [], 'FishGeno',[], ...
        'BoutFrequency', [], 'BoutDuration',[],'NumberOfOscillations',[],'BoutDistance',[],'Speed',[],'TBF',[],'MaxBendAmp',[],'MedianBendAmp',[]);


%% Selecte good swimers
GoodSwimers=find(~([datasetPerFish(:).MeanIBI]> 3.3));
datasetGoodFish=datasetPerFish(GoodSwimers);

%% define dataset

datasetPerBout=DatasetPreEscape_FWD1;

%% find fish number and genotypes
Fish = unique([datasetGoodFish(:).Condition]);
NumberFish=length(Fish);

% Define genotype
Fish_temp=Fish;
FishGeno=([datasetGoodFish.Genotype]);
Fish_ID = ([datasetGoodFish.Condition]);
NTrial = ([datasetGoodFish.NTrial]);

Fish_G2=Fish_temp(find(FishGeno( find( Fish_temp ) )==2));
Fish_G1=Fish_temp(find(FishGeno( find( Fish_temp ) )==1));
Fish_G0=Fish_temp(find(FishGeno( find( Fish_temp ) )==0));

%% calculate for each fish parameters extracted from the structure

for i=1:NumberFish;

    %index for IBI
    
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
    
    % All index
    allindex{Fish(i)}= find(~([datasetPerBout(:).Condition]-Fish(i)));

    % now calculate parameters for each bout 
    if length(allindex{Fish(i)})==0;
        BoutDuration{Fish(i)}=nan;
        BoutDistance{Fish(i)}=nan;
        Speed{Fish(i)}=nan;
        NumberOfOscillations{Fish(i)}=nan;
        TBF{Fish(i)}=nan;
        MaxBendAmp{Fish(i)}=nan;
        MedianBendAmp{Fish(i)}=nan;
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
    end
    
    
    % outputTable
    output(i).NTrial=NTrial(i); 
    output(i).Fish_ID=Fish_ID(i);
    output(i).FishGeno=FishGeno(i); 
    
    output(i).BoutFrequency=medianBoutFrequency{Fish(i)};
    output(i).BoutDuration=median(BoutDuration{Fish(i)},'omitnan');
    output(i).BoutDistance=median(BoutDistance{Fish(i)},'omitnan');
    output(i).Speed=median(Speed{Fish(i)},'omitnan');
    output(i).NumberOfOscillations=mean(NumberOfOscillations{Fish(i)},'omitnan');
    output(i).TBF=median(TBF{Fish(i)},'omitnan');
    output(i).MaxBendAmp=median(MaxBendAmp{Fish(i)},'omitnan');
    output(i).MedianBendAmp=median(MedianBendAmp{Fish(i)},'omitnan');
    
    
    i= i+1;

end;

%% Generation of text Table file for statistic analysis by Francois-Xavier

% Table=struct2table(output);
% write(table);

%% Generation dataset for Fish_WT, Fish_Homo

Fish_WT=output(find([output.FishGeno]==2));

Fish_Homo=output(find([output.FishGeno]==0));

%% subplot pre_Stimulus
h1=figure(1); 

title('BoutFrequency(Hz)');hold on;
%boxplot_2gp([Fish_WT.BoutFrequency],[Fish_Homo.BoutFrequency]);hold on;
MeanErrorbar([Fish_WT.BoutFrequency],[Fish_Homo.BoutFrequency]);hold on;
%ylim([0.1 0.5]);
hold off;

saveas(h1,['BoutFrequency(Hz)_mean.fig'])
saveas(h1,['BoutFrequency(Hz)_mean.epsc'])

%%
h2=figure(2); 
subplot(2,3,1)
title('Bout Distance (mm)');hold on;
%boxplot_2gp([Fish_WT.BoutDistance],[Fish_Homo.BoutDistance]);hold on;
MeanErrorbar([Fish_WT.BoutDistance],[Fish_Homo.BoutDistance]);hold on;
%ylim([0 7]);
%ylim([0 10]);
hold off;

subplot(2,3,2)
title('Bout Duration (sec)');hold on; 
%boxplot_2gp([Fish_WT.BoutDuration],[Fish_Homo.BoutDuration]);hold on;
MeanErrorbar([Fish_WT.BoutDuration],[Fish_Homo.BoutDuration]);hold on;
%ylim([0 2]);
hold off;

subplot(2,3,3)
title('Bout Speed (mm/sec)');hold on;
%boxplot_2gp([Fish_WT.Speed],[Fish_Homo.Speed]);hold on;
MeanErrorbar([Fish_WT.Speed],[Fish_Homo.Speed]);hold on;
%ylim([0 5]);
%ylim([0 8]);
hold off;

subplot(2,3,4)
title('# Of Oscillations'); hold on;
%boxplot_2gp([Fish_WT.NumberOfOscillations],[Fish_Homo.NumberOfOscillations]);hold on;
MeanErrorbar([Fish_WT.NumberOfOscillations],[Fish_Homo.NumberOfOscillations]);hold on;
%ylim([-2 5]);
hold off;

subplot(2,3,5)
title('TBF');hold on;  
%boxplot_2gp([Fish_WT.TBF],[Fish_Homo.TBF]);hold on;
MeanErrorbar([Fish_WT.TBF],[Fish_Homo.TBF]);hold on;
hold off;

% subplot(2,3,6)
% title('Max Bend Amplitude (degree)');hold on;  
% boxplot_2gp([Fish_WT.MaxBendAmp],[Fish_Homo.MaxBendAmp]);hold on;
% MeanErrorbar([Fish_WT.MaxBendAmp],[Fish_Homo.MaxBendAmp]);hold on;
% %ylim([0 40]);
% hold off;

subplot(2,3,6)
title('Median Bend Amplitude (degree)');hold on;  
%boxplot_2gp([Fish_WT.MedianBendAmp],[Fish_Homo.MedianBendAmp]);hold on;
MeanErrorbar([Fish_WT.MedianBendAmp],[Fish_Homo.MedianBendAmp]);hold on;
%ylim([0 40]);
hold off;


saveas(h2,['ParametersOverview_FWD_mean.fig'])
saveas(h2,['ParametersOverview_FWD_mean.epsc'])

