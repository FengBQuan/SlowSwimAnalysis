%% Clean workspace and load data
close all
clear variables
clc
load("DataSetPreEscape_FWDvsRT_thres25_new.mat")
%% create output stucture
    
output = struct( 'NTrial', [], 'Fish_ID', [], 'FishGeno',[], ...
         'BoutDuration',[],'NumberOfOscillations',[],'BoutDistance',[],'Speed',[],'TBF',[],'MedianBendAmp',[],'MaxBendAmp',[]);
    
%% define dataset

%for thres 25
%datasetPerBout=datasetPerBout_FWD1;
datasetPerBout=datasetPerBout_RT1;
%datasetPerBout=DatasetPreEscape;

%% set fish number and fish ID
Fish = unique([datasetPerBout(:).Condition]);
NumberFish=length(Fish);


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


    %Calculate IBI and BoutFrequency=1/IBI by using index without 1st one
    
    allindex{Fish(i)}= find(~([datasetPerBout(:).Condition]-Fish(i)));
  
    
    FishGeno{Fish(i)}=unique([datasetPerBout(allindex{Fish(i)}).Genotype]);
    Fish_ID{Fish(i)} =unique([datasetPerBout(allindex{Fish(i)}).Condition]);
    NTrial{Fish(i)} =unique([datasetPerBout(allindex{Fish(i)}).NTrial]); % Ntiral= N of Clutch

    % if no tracking data, allindex =[], parameters=nan.
    if length(allindex{Fish(i)})==0 ;
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
    output(i).NTrial=NTrial{Fish(i)}; 
    output(i).Fish_ID=Fish_ID{Fish(i)};
    output(i).FishGeno=FishGeno{Fish(i)}; 
    
    %output(i).BoutRate=BoutRate{Fish(i)};% nBouts/DurationTotale
    output(i).BoutDuration=median(BoutDuration{Fish(i)},'omitnan');
    output(i).BoutDistance=median(BoutDistance{Fish(i)},'omitnan');
    output(i).Speed=median(Speed{Fish(i)},'omitnan');
    output(i).NumberOfOscillations=mean(NumberOfOscillations{Fish(i)},'omitnan');
    output(i).TBF=median(TBF{Fish(i)},'omitnan');
    output(i).MedianBendAmp=median(MedianBendAmp{Fish(i)},'omitnan');
    output(i).MaxBendAmp=median(MaxBendAmp{Fish(i)},'omitnan');
   %output(i).BoutFrequency=medianBoutFrequency{Fish(i)};%1/IBI
    
    i= i+1;

end;


%% Generation dataset for Fish_WT, Fish_Homo

Fish_WT=output(find([output.FishGeno]==2));

Fish_Homo=output(find([output.FishGeno]==0));

%save('5min_FWD_thres25','FWD_Fish_WT','FWD_Fish_Homo');

%% Generation of text Table file for statistic analysis by Francois-Xavier

Table_RT=struct2table(output);
writetable(Table_RT);

%save('5min_GoodSwimmers_All_thres25');