%% Clean workspace and load data
close all
clear variables
clc
%load("matlab_workspace_SST_4manip_SlowSwim.mat")
%% load the dataset file
InitialFolder='/Users/bong-iquan/MATLAB/sst1.1_mutant/Slow_swim';
 
% Figure out how many trials to analyse
cd(InitialFolder)
list=dir('*_workspace_SST_4manip_SlowSwim*');
%list2=list([list.isdir]==0);%list all folders but no file (make sure I have only folders) 1==folder; 0==file;
NFolder=size(list,1);

%% Initiate variables 
% collecte BoutDuration for all trials 
Compt_NTrial=0;

BD_AllTrials_pre=[];BD_AllTrials_post=[];Osc_AllTrials_pre=[];Osc_AllTrials_post=[];DisB_AllTrials_pre=[];DisB_AllTrials_post=[];
Speed_AllTrials_pre=[];Speed_AllTrials_post=[];TBF_AllTrials_pre=[];TBF_AllTrials_post=[];IBI_AllTrials_pre=[];IBI_AllTrials_post=[];
nBout_AllTrials_pre=[];nBout_AllTrials_post=[];Amplitude_AllTrials_pre=[];Amplitude_AllTrials_post=[];

Median_Amplitude=[];

for ii=1:NFolder
    
    % Selection manually!
       %[fileName,pathName] = uigetfile('.mat');
    
%Selection automatique dans l ordre alphabetique
    
    fileName =list(ii).name ;
    pathName=strcat(list(ii).folder,'/');
    load(strcat(pathName,fileName));
    
output = struct( 'NTrial', [], 'Fish_ID', [], 'FishGeno',[], ...
        'mIBI_pre', [], 'mBoutDuration_pre',[],'mNumberOfOscillations_pre',[],'mBoutDistance_pre',[],'mSpeed_pre',[],'mnBout_pre', [], 'mTBF_pre',[],'mAmplitude_pre',[],...
    'mIBI_post', [], 'mBoutDuration_post',[],'mNumberOfOscillations_post',[],'mBoutDistance_post',[],'mSpeed_post',[], 'mnBout_post',[], 'mTBF_post',[],'mAmplitude_post',[]);
%     'IBI_pre', [], 'TimeBout_pre', [],'BoutDuration_pre',[],'NumberOfOscillations_pre',[],'BoutDistance_pre',[],'Speed_pre',[],'nBout_pre', [], 'TBF_pre',[],...
%     'IBI_post', [], 'TimeBout_post', [],'BoutDuration_post',[],'NumberOfOscillations_post',[],'BoutDistance_post',[],'Speed_post',[], 'nBout_post',[], 'TBF_post',[]);

close all;
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
% find fish number and genotypes

Fish = unique([datasetPerFish(:).Condition]);
NumberFish=length(Fish);
Genotypes = unique([datasetPerBout(:).Genotype]);

% find EscapeWindow

EscapeWindow = [unique([datasetPerBout(:).EscapeWindow1]) unique([datasetPerBout(:).EscapeWindow2])];

%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
% calculate for each fish parameters extracted from the structure

for i=1:NumberFish;
    i
    % Define genotype
Fish_temp=Fish; % to use datasetPerFish
FishGeno=([datasetPerFish.Genotype]);
Fish_ID = ([datasetPerFish.Condition]);
NTrial = ([datasetPerFish.NTrial]);

Fish_G2=Fish_temp(find(FishGeno( find( Fish_temp ) )==2));
Fish_G1=Fish_temp(find(FishGeno( find( Fish_temp ) )==1));
Fish_G0=Fish_temp(find(FishGeno( find( Fish_temp ) )==0));
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%

    % index for IBI
    
    index{Fish(i)}= find(~([datasetPerBout(:).Condition]-Fish(i)));
    index{Fish(i)}(1)=[];
    
    BoutsIBI_Pre_Escape= index{Fish(i)}(find([datasetPerBout(index{Fish(i)}).BoutStart]< EscapeWindow(1)));
    BoutsIBI_Post_Escape= index{Fish(i)}(find([datasetPerBout(index{Fish(i)}).BoutStart]> EscapeWindow(2)));
    
    % all index
    
    %geno_index{Fish(i)} = unique([datasetPerBout([index{Fish(i)}]).Genotype]);
    allindex{Fish(i)}= find(~([datasetPerBout(:).Condition]-Fish(i)));
    
    Bouts_Pre_Escape= allindex{Fish(i)}(find([datasetPerBout(allindex{Fish(i)}).BoutStart]< EscapeWindow(1)));
    Bouts_Post_Escape= allindex{Fish(i)}(find([datasetPerBout(allindex{Fish(i)}).BoutStart]>EscapeWindow(2)));
    
    datasetPreEscape=datasetPerBout(Bouts_Pre_Escape);
    datasetPostEscape=datasetPerBout(Bouts_Post_Escape);
    
    idx_PreEscape{Fish(i)}= find(~([datasetPreEscape(:).Condition]-Fish(i)));
    idx_PostEscape{Fish(i)}= find(~([datasetPostEscape(:).Condition]-Fish(i)));
    
    % now calculate parameters for each bout 
    
%     for h=1:length(allindex{Fish(i)});
%         
%         display(['currently processing fish ' num2str(i)])
%         display(['currently processing bout number ' num2str(h)])
%         
%         % need to convert in degrees... otherwise : TailAngle{Fish(i)}{h}=[datasetPerBout(allindex{Fish(i)}(h)).TailAngle_smoothed]';
%         
%         TailAngle{Fish(i)}{h}=57.2958*[datasetPerBout(allindex{Fish(i)}(h)).TailAngle_smoothed]';
%         
%         Bend_Timing{Fish(i)}{h}=[datasetPerBout(allindex{Fish(i)}(h)).Bend_Timing];
%         
%         Bend_Amplitude{Fish(i)}{h} = 57.2958*datasetPerBout(allindex{Fish(i)}(h)).Bend_Amplitude;
%         
%         AmpMedianPerBout{Fish(i)}{h}=median(abs(Bend_Amplitude{Fish(i)}{h}));
%      
%         TBF{Fish(i)}{h} = [datasetPerBout(allindex{Fish(i)}(h)).InstantaneousTBF];
%       
%         TBFmedianPerBout{Fish(i)}{h}=median(TBF{Fish(i)}{h});
%         
%     end;
    
        %calculate Pre_Escape
    IBI_pre{Fish(i)} = [datasetPerBout(BoutsIBI_Pre_Escape).InstantaneousIBI];
    TimeBout_pre{Fish(i)} = [datasetPerBout(BoutsIBI_Pre_Escape).BoutStart];
    
    BoutDuration_pre{Fish(i)}=[datasetPerBout(Bouts_Pre_Escape).BoutDuration];
    NumberOfOscillations_pre{Fish(i)}=[datasetPerBout(Bouts_Pre_Escape).NumberOfOscillations];
    BoutDistance_pre{Fish(i)}=[datasetPerBout(Bouts_Pre_Escape).TotalDistance];
    Speed_pre{Fish(i)}=[datasetPerBout(Bouts_Pre_Escape).Speed];
    nBout_pre{Fish(i)}= length(Bouts_Pre_Escape);
    
    if length(idx_PreEscape{Fish(i)})==0;
        Pre_TailAngle{Fish(i)}=0;
        Pre_Bend_Amplitude{Fish(i)}= 0;
        Pre_AmpMedianPerBout{Fish(i)}=0;
        Pre_TBF{Fish(i)}= 0;
        Pre_TBFmedianPerBout{Fish(i)}= 0;
        
    else
        for h=1:length(idx_PreEscape{Fish(i)});
            
            display(['currently processing fish ' num2str(i)])
            display(['currently processing bout number ' num2str(h)])
            
            % need to convert in degrees... otherwise : TailAngle{Fish(i)}{h}=[datasetPerBout(allindex{Fish(i)}(h)).TailAngle_smoothed]';
            
            
            Pre_TailAngle{Fish(i)}{h}=57.2958*[datasetPreEscape(idx_PreEscape{Fish(i)}(h)).TailAngle_smoothed]';
            Pre_Bend_Timing{Fish(i)}{h}=[datasetPreEscape(idx_PreEscape{Fish(i)}(h)).Bend_Timing];
            Pre_Bend_Amplitude{Fish(i)}{h} = 57.2958*datasetPreEscape(idx_PreEscape{Fish(i)}(h)).Bend_Amplitude;
            Pre_AmpMedianPerBout{Fish(i)}(h)=median(abs(Pre_Bend_Amplitude{Fish(i)}{h}));
            Pre_TBF{Fish(i)}{h} = [datasetPreEscape(idx_PreEscape{Fish(i)}(h)).InstantaneousTBF];
            Pre_TBFmedianPerBout{Fish(i)}(h)=median(Pre_TBF{Fish(i)}{h});
        end
        
    end
    
    
         %calculate Post_Escape    
    IBI_post{Fish(i)} = [datasetPerBout(BoutsIBI_Post_Escape).InstantaneousIBI];
    TimeBout_post{Fish(i)} = [datasetPerBout(BoutsIBI_Post_Escape).BoutStart];
    
    BoutDuration_post{Fish(i)}=[datasetPerBout(Bouts_Post_Escape).BoutDuration];
    NumberOfOscillations_post{Fish(i)}=[datasetPerBout(Bouts_Post_Escape).NumberOfOscillations];
    BoutDistance_post{Fish(i)}=[datasetPerBout(Bouts_Post_Escape).TotalDistance];
    Speed_post{Fish(i)}=[datasetPerBout(Bouts_Post_Escape).Speed];
    nBout_post{Fish(i)}=length(Bouts_Post_Escape);
    
    if length(idx_PostEscape{Fish(i)})==0;
        Post_TailAngle{Fish(i)}=0;
        Post_Bend_Amplitude{Fish(i)} = 0;
        Post_AmpMedianPerBout{Fish(i)}=0;
        Post_TBF{Fish(i)}{h} = 0;
        Post_TBFmedianPerBout{Fish(i)}= 0;
        
    else
        for h=1:length(idx_PostEscape{Fish(i)});
            
            display(['currently processing fish ' num2str(i)])
            display(['currently processing bout number ' num2str(h)])
            
            % need to convert in degrees... otherwise : TailAngle{Fish(i)}{h}=[datasetPerBout(allindex{Fish(i)}(h)).TailAngle_smoothed]';
            
            Post_TailAngle{Fish(i)}{h}=57.2958*[datasetPostEscape(idx_PostEscape{Fish(i)}(h)).TailAngle_smoothed]';
            Post_Bend_Timing{Fish(i)}{h}=[datasetPostEscape(idx_PostEscape{Fish(i)}(h)).Bend_Timing];
            Post_Bend_Amplitude{Fish(i)}{h} = 57.2958*datasetPostEscape(idx_PostEscape{Fish(i)}(h)).Bend_Amplitude;
            Post_AmpMedianPerBout{Fish(i)}(h)=median(abs(Post_Bend_Amplitude{Fish(i)}{h}));
            Post_TBF{Fish(i)}{h} = [datasetPostEscape(idx_PostEscape{Fish(i)}(h)).InstantaneousTBF];
            Post_TBFmedianPerBout{Fish(i)}(h)=median(Post_TBF{Fish(i)}{h});
        end
        
    end
           
    % calculate parameters for all swims     
%       IBI{Fish(i)} = [datasetPerBout(index{Fish(i)}).InstantaneousIBI];
%       TimeBoutIBI{Fish(i)} = [datasetPerBout(index{Fish(i)}).BoutStart];
%       
%       TimeBout{Fish(i)} =[datasetPerBout(allindex{Fish(i)}).BoutStart];
%       nBout{Fish(i)}=(nBout_pre{Fish(i)})+(nBout_post{Fish(i)});
      
%     BoutDuration{Fish(i)}=[datasetPerBout(allindex{Fish(i)}).BoutDuration];
%     Distance{Fish(i)}=[datasetPerBout(allindex{Fish(i)}).TotalDistance];
%     Speed{Fish(i)}=[datasetPerBout(allindex{Fish(i)}).Speed];
%     NumberOfOscillations{Fish(i)}=[datasetPerBout(allindex{Fish(i)}).NumberOfOscillations];
 
    output(i).NTrial=NTrial(i); 
    output(i).Fish_ID=Fish_ID(i);
    output(i).FishGeno=FishGeno(i); 
    
    output(i).mIBI_pre=log(median(IBI_pre{Fish(i)},'omitnan'));
    output(i).mBoutDuration_pre=median(BoutDuration_pre{Fish(i)},'omitnan');
    output(i).mNumberOfOscillations_pre=median(NumberOfOscillations_pre{Fish(i)},'omitnan');
    output(i).mSpeed_pre=median(Speed_pre{Fish(i)},'omitnan');
    output(i).mBoutDistance_pre=median(BoutDistance_pre{Fish(i)},'omitnan');
    output(i).mnBout_pre=median(nBout_pre{Fish(i)},'omitnan');
    output(i).mTBF_pre=median(Pre_TBFmedianPerBout{Fish(i)},'omitnan');
    output(i).mAmplitude_pre=median(Pre_AmpMedianPerBout{Fish(i)},'omitnan');
    
    output(i).mIBI_post=log(median(IBI_post{Fish(i)},'omitnan'));
    output(i).mBoutDuration_post=median(BoutDuration_post{Fish(i)},'omitnan');
    output(i).mNumberOfOscillations_post=median(NumberOfOscillations_post{Fish(i)},'omitnan');
    output(i).mSpeed_post=median(Speed_post{Fish(i)},'omitnan');
    output(i).mBoutDistance_post=median(BoutDistance_post{Fish(i)},'omitnan');
    output(i).mnBout_post=median(nBout_post{Fish(i)},'omitnan');
    output(i).mTBF_post=median(Post_TBFmedianPerBout{Fish(i)},'omitnan');
    output(i).mAmplitude_post=median(Post_AmpMedianPerBout{Fish(i)},'omitnan');
    
 %output works only for one folder
    
    i= i+1;

end;



%% Collecte data from all Trials

BD_AllTrials_pre{1,1+Compt_NTrial}=BoutDuration_pre;
BD_AllTrials_post{1,1+Compt_NTrial}=BoutDuration_post;

Osc_AllTrials_pre{1,1+Compt_NTrial}=NumberOfOscillations_pre;
Osc_AllTrials_post{1,1+Compt_NTrial}=NumberOfOscillations_post;

DisB_AllTrials_pre{1,1+Compt_NTrial}=BoutDistance_pre;
DisB_AllTrials_post{1,1+Compt_NTrial}=BoutDistance_post;

Speed_AllTrials_pre{1,1+Compt_NTrial}=Speed_pre;
Speed_AllTrials_post{1,1+Compt_NTrial}=Speed_post;

IBI_AllTrials_pre{1,1+Compt_NTrial}=IBI_pre;
IBI_AllTrials_post{1,1+Compt_NTrial}=IBI_post;

nBout_AllTrials_pre{1,1+Compt_NTrial}=nBout_pre;
nBout_AllTrials_post{1,1+Compt_NTrial}=nBout_post;

TBF_AllTrials_pre{1,1+Compt_NTrial}=Pre_TBFmedianPerBout; 
TBF_AllTrials_post{1,1+Compt_NTrial}=Post_TBFmedianPerBout;

Amplitude_AllTrials_pre{1,1+Compt_NTrial}=Pre_AmpMedianPerBout;
Amplitude_AllTrials_post{1,1+Compt_NTrial}=Post_AmpMedianPerBout;

Compt_NTrial=Compt_NTrial+1;
 
%% mean Calculation

for l=1:length(Fish_G2);
    l
    % Calculate medianBoutDuration
    G2medianBoutDuration_pre{ii}(l)=median(BD_AllTrials_pre{ii}{(Fish_G2(l))});
    G2medianBoutDuration_post{ii}(l)=median(BD_AllTrials_post{ii}{(Fish_G2(l))});    
   
    % Calculate medianNumber of Oscillation
    G2medianNumOfOsc_pre{ii}(l)=median(Osc_AllTrials_pre{ii}{(Fish_G2(l))});
    G2medianNumOfOsc_post{ii}(l)=median(Osc_AllTrials_post{ii}{(Fish_G2(l))}); 
     % Calculate meanNumber of Oscillation
    G2meanNumOfOsc_pre{ii}(l)=mean(Osc_AllTrials_pre{ii}{(Fish_G2(l))});
    G2meanNumOfOsc_post{ii}(l)=mean(Osc_AllTrials_post{ii}{(Fish_G2(l))});
    
    % Calculate medianBoutDistance 
    G2medianBoutDistance_pre{ii}(l)=median(DisB_AllTrials_pre{ii}{(Fish_G2(l))});
    G2medianBoutDistance_post{ii}(l)=median(DisB_AllTrials_post{ii}{(Fish_G2(l))});
    
    % Calculate medianSpeed
    G2medianSpeed_pre{ii}(l)=median(Speed_AllTrials_pre{ii}{(Fish_G2(l))});
    G2medianSpeed_post{ii}(l)=median(Speed_AllTrials_post{ii}{(Fish_G2(l))});
    
    % Calculate median TBF 
    G2medianTBF_pre{ii}(l)=median(TBF_AllTrials_pre{ii}{(Fish_G2(l))});
    G2medianTBF_post{ii}(l)=median(TBF_AllTrials_post{ii}{(Fish_G2(l))}); 
        % Calculate mean TBF 
    G2meanTBF_pre{ii}(l)=mean(TBF_AllTrials_pre{ii}{(Fish_G2(l))});
    G2meanTBF_post{ii}(l)=mean(TBF_AllTrials_post{ii}{(Fish_G2(l))}); 
    
    % Calculate median IBI
    G2medianIBI_pre{ii}(l)=log(median(IBI_AllTrials_pre{ii}{(Fish_G2(l))}));
    G2medianIBI_post{ii}(l)=log(median(IBI_AllTrials_post{ii}{(Fish_G2(l))}));
     
    % Calculate median nBout
    G2median_nBout_pre{ii}(l)=nBout_AllTrials_pre{ii}{(Fish_G2(l))};
    G2median_nBout_post{ii}(l)=nBout_AllTrials_post{ii}{(Fish_G2(l))};
    
     % Calculate median Amplitude 
    G2medianAmpmitude_pre{ii}(l)=median(Amplitude_AllTrials_pre{ii}{(Fish_G2(l))});
    G2medianAmpmitude_post{ii}(l)=median(Amplitude_AllTrials_post{ii}{(Fish_G2(l))}); 
 
end;
 
 
%Heterozygote
for l=1:length(Fish_G1)
    l
    % Calculate medianBoutDuration
    G1medianBoutDuration_pre{ii}(l)=median(BD_AllTrials_pre{ii}{(Fish_G1(l))});
    G1medianBoutDuration_post{ii}(l)=median(BD_AllTrials_post{ii}{(Fish_G1(l))});
    
    % Calculate medianNumber of Oscillation
    G1medianNumOfOsc_pre{ii}(l)=median(Osc_AllTrials_pre{ii}{(Fish_G1(l))});
    G1medianNumOfOsc_post{ii}(l)=median(Osc_AllTrials_post{ii}{(Fish_G1(l))});
    % Calculate meanNumber of Oscillation
    G1meanNumOfOsc_pre{ii}(l)=mean(Osc_AllTrials_pre{ii}{(Fish_G1(l))});
    G1meanNumOfOsc_post{ii}(l)=mean(Osc_AllTrials_post{ii}{(Fish_G1(l))});
    
    % Calculate medianBoutDistance
    
    G1medianBoutDistance_pre{ii}(l)=median(DisB_AllTrials_pre{ii}{(Fish_G1(l))});
    G1medianBoutDistance_post{ii}(l)=median(DisB_AllTrials_post{ii}{(Fish_G1(l))});
    
    % Calculate medianSpeed
    G1medianSpeed_pre{ii}(l)=median(Speed_AllTrials_pre{ii}{(Fish_G1(l))});
    G1medianSpeed_post{ii}(l)=median(Speed_AllTrials_post{ii}{(Fish_G1(l))});
 
    % Calculate median IBI
    G1medianIBI_pre{ii}(l)=log(median(IBI_AllTrials_pre{ii}{(Fish_G1(l))}));
    G1medianIBI_post{ii}(l)=log(median(IBI_AllTrials_post{ii}{(Fish_G1(l))}));
    
    % Calculate median TBF
    G1medianTBF_pre{ii}(l)=median(TBF_AllTrials_pre{ii}{(Fish_G1(l))});
    G1medianTBF_post{ii}(l)=median(TBF_AllTrials_post{ii}{(Fish_G1(l))});
            % Calculate mean TBF 
    G1meanTBF_pre{ii}(l)=mean(TBF_AllTrials_pre{ii}{(Fish_G1(l))});
    G1meanTBF_post{ii}(l)=mean(TBF_AllTrials_post{ii}{(Fish_G1(l))}); 

    % Calculate median nBout
    G1median_nBout_pre{ii}(l)=nBout_AllTrials_pre{ii}{(Fish_G1(l))};
    G1median_nBout_post{ii}(l)=nBout_AllTrials_post{ii}{(Fish_G1(l))};
    
         % Calculate median Amplitude 
    G1medianAmpmitude_pre{ii}(l)=median(Amplitude_AllTrials_pre{ii}{(Fish_G1(l))});
    G1medianAmpmitude_post{ii}(l)=median(Amplitude_AllTrials_post{ii}{(Fish_G1(l))}); 
 
   
end;
 
 
%Homozygote
for l=1:length(Fish_G0)
    l
    % Calculate medianBoutDuration
    G0medianBoutDuration_pre{ii}(l)=median(BD_AllTrials_pre{ii}{(Fish_G0(l))});
    G0medianBoutDuration_post{ii}(l)=median(BD_AllTrials_post{ii}{(Fish_G0(l))});   
   
    % Calculate medianNumber of Oscillation
    G0medianNumOfOsc_pre{ii}(l)=median(Osc_AllTrials_pre{ii}{(Fish_G0(l))});
    G0medianNumOfOsc_post{ii}(l)=median(Osc_AllTrials_post{ii}{(Fish_G0(l))});
         % Calculate meanNumber of Oscillation
    G0meanNumOfOsc_pre{ii}(l)=mean(Osc_AllTrials_pre{ii}{(Fish_G0(l))});
    G0meanNumOfOsc_post{ii}(l)=mean(Osc_AllTrials_post{ii}{(Fish_G0(l))});
    
    % Calculate medianBoutDistance   
    G0medianBoutDistance_pre{ii}(l)=median(DisB_AllTrials_pre{ii}{(Fish_G0(l))});
    G0medianBoutDistance_post{ii}(l)=median(DisB_AllTrials_post{ii}{(Fish_G0(l))});
    
    % Calculate medianSpeed
    G0medianSpeed_pre{ii}(l)=median(Speed_AllTrials_pre{ii}{(Fish_G0(l))});
    G0medianSpeed_post{ii}(l)=median(Speed_AllTrials_post{ii}{(Fish_G0(l))});
    
    % Calculate median TBF
    G0medianTBF_pre{ii}(l)=median(TBF_AllTrials_pre{ii}{(Fish_G0(l))});
    G0medianTBF_post{ii}(l)=median(TBF_AllTrials_post{ii}{(Fish_G0(l))}); 
            % Calculate mean TBF 
    G0meanTBF_pre{ii}(l)=mean(TBF_AllTrials_pre{ii}{(Fish_G0(l))});
    G0meanTBF_post{ii}(l)=mean(TBF_AllTrials_post{ii}{(Fish_G0(l))}); 
    
    % Calculate median IBI
    G0medianIBI_pre{ii}(l)=log(median(IBI_AllTrials_pre{ii}{(Fish_G0(l))}));
    G0medianIBI_post{ii}(l)=log(median(IBI_AllTrials_post{ii}{(Fish_G0(l))}));
     
    % Calculate median nBout
    G0median_nBout_pre{ii}(l)=nBout_AllTrials_pre{ii}{(Fish_G0(l))};
    G0median_nBout_post{ii}(l)=nBout_AllTrials_post{ii}{(Fish_G0(l))};
    
    % Calculate median Amplitude 
    G0medianAmpmitude_pre{ii}(l)=median(Amplitude_AllTrials_pre{ii}{(Fish_G0(l))});
    G0medianAmpmitude_post{ii}(l)=median(Amplitude_AllTrials_post{ii}{(Fish_G0(l))}); 
 
end;

ii=ii+1;

end;

%% subplot pre_Stimulus
h1=figure(1); 
title('Pre_Stimulus');hold on;

subplot(2,4,1)
title('Bout Duration (sec)');hold on;
%boxplot_2groupes(G2medianBoutDuration_pre,G2medianBoutDuration_post);hold on;
boxplot_2groupes(G2medianBoutDuration_pre,G0medianBoutDuration_pre);hold on;
ylim([0.1 0.6]);hold off;

subplot(2,4,2)
title('Oscillations (n)'); hold on;
%boxplot_2groupes(G2medianNumOfOsc_pre,G2medianNumOfOsc_post);hold on;
boxplot_2groupes(G2meanNumOfOsc_pre,G2meanNumOfOsc_post);hold on;
%boxplot_2groupes(G2medianNumOfOsc_pre,G0medianNumOfOsc_pre);hold on;
ylim([0 7]);hold off;

subplot(2,4,3)
title('Bout Distance (mm)');hold on;
%boxplot_2groupes(G2medianBoutDistance_pre,G2medianBoutDistance_post);hold on;
boxplot_2groupes(G2medianBoutDistance_pre,G0medianBoutDistance_pre);hold on;
ylim([0 2]);hold off;

subplot(2,4,4)
title('Bout Speed (mm/sec)');hold on;
%boxplot_2groupes(G2medianSpeed_pre,G2medianSpeed_post);hold on;
boxplot_2groupes(G2medianSpeed_pre,G0medianSpeed_pre);hold on;
ylim([0 5]);
hold off;

subplot(2,4,5)
title('IBI (log)');hold on;  
%boxplot_2groupes(G2medianIBI_pre,G2medianIBI_post); hold on;
boxplot_2groupes(G2medianIBI_pre,G0medianIBI_pre); hold on;
ylim([-2 5]);
hold off;

subplot(2,4,6)
title('nBout');hold on;  
%boxplot_2groupes(G2median_nBout_pre,G2median_nBout_post); hold on;
boxplot_2groupes(G2median_nBout_pre,G0median_nBout_pre); hold on;
hold off;

subplot(2,4,7)
title('TBF');hold on;  
%boxplot_2groupes(G2medianTBF_pre,G2medianTBF_post); hold on;
boxplot_2groupes(G2meanTBF_pre,G0meanTBF_pre); hold on;
%boxplot_2groupes(G2medianTBF_pre,G0medianTBF_pre); hold on;
hold off;

subplot(2,4,8)
title('Amplitude in degree');hold on;  
%boxplot_2groupes(G2medianAmpmitude_pre,G2medianAmpmitude_post); hold on;
boxplot_2groupes(G2medianAmpmitude_pre,G0medianAmpmitude_pre); hold on;
hold off;


%% Post_Stimulus subplot 
h2=figure(2); 
title('Post_Stimulus');hold on;

subplot(2,4,1)
title('Bout Duration (sec)');hold on;
%boxplot_2groupes(G2medianBoutDuration_post,G0medianBoutDuration_post);hold on;
boxplot_4groupes(G2medianBoutDuration_pre,G2medianBoutDuration_post,G0medianBoutDuration_pre,G0medianBoutDuration_post);hold on;
ylim([0.1 0.6]);hold off;

subplot(2,4,2)
title('Oscillations (n)'); hold on;
%boxplot_2groupes(G2medianNumOfOsc_post,G0medianNumOfOsc_post);hold on;
%boxplot_4groupes(G2medianNumOfOsc_pre,G2medianNumOfOsc_post,G0medianNumOfOsc_pre,G0medianNumOfOsc_post);hold on;
boxplot_4groupes(G2meanNumOfOsc_pre,G2meanNumOfOsc_post,G0meanNumOfOsc_pre,G0meanNumOfOsc_post);hold on;
ylim([0 7]);hold off;

subplot(2,4,3)
title('Bout Distance (mm)');hold on;
%boxplot_2groupes(G2medianBoutDistance_post,G0medianBoutDistance_post);hold on;
boxplot_4groupes(G2medianBoutDistance_pre,G2medianBoutDistance_post,G0medianBoutDistance_pre,G0medianBoutDistance_post);hold on;
ylim([0 2]);hold off;

subplot(2,4,4)
title('Bout Speed (mm/sec)');hold on;
%boxplot_2groupes(G2medianSpeed_post,G0medianSpeed_post);hold on;
boxplot_4groupes(G2medianSpeed_pre,G2medianSpeed_post,G0medianSpeed_pre,G0medianSpeed_post);hold on;
ylim([0 5]);
hold off;

subplot(2,4,5)
title('IBI (sec)');hold on;   
%boxplot_2groupes(G2medianIBI_post,G0medianIBI_post); hold on;
boxplot_4groupes(G2medianIBI_pre,G2medianIBI_post,G0medianIBI_pre,G0medianIBI_post); hold on;
ylim([-2 5]);
hold off;

subplot(2,4,6)
title('nBout');hold on;  
%boxplot_2groupes(G2median_nBout_pre,G2median_nBout_post); hold on;
boxplot_4groupes(G2median_nBout_pre,G2median_nBout_post,G0median_nBout_pre,G0median_nBout_post); hold on;
hold off;

subplot(2,4,7)
title('TBF');hold on;  
%boxplot_2groupes(G2medianTBF_pre,G2medianTBF_post); hold on;
%boxplot_2groupes(G2medianTBF_pre,G2medianTBF_post,G0medianTBF_pre,G0medianTBF_post); hold on;
boxplot_4groupes(G2meanTBF_pre,G2meanTBF_post,G0meanTBF_pre,G0meanTBF_post); hold on;
hold off;

subplot(2,4,8)
title('Amplitude in degree');hold on;  
%boxplot_2groupes(G2medianAmpmitude_pre,G2medianAmpmitude_post); hold on;
boxplot_4groupes(G2medianAmpmitude_pre,G2medianAmpmitude_post,G0medianAmpmitude_pre,G0medianAmpmitude_post); hold on;
hold off;



