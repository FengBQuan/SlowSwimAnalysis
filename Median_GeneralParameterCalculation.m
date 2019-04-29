 %% Clean workspace and load data
clear
%load('matlab sst_20191101_protocole1.mat')

%% load the dataset file
InitialFolder='/Users/bong-iquan/MATLAB/sst1.1_mutant/Slow_swim';
 
% Figure out how many trials to analyse
cd(InitialFolder)
list=dir('*_workspace_SST_*');
%list2=list([list.isdir]==0);%list all folders but no file (make sure I have only folders) 1==folder; 0==file;
NFolder=size(list,1);

%% Initiate variables 
Compt_NTrial=0;

BD_AllTrials_pre=[];
BD_AllTrials_post=[];

Osc_AllTrials_pre=[];
Osc_AllTrials_post=[];

DisB_AllTrials_pre=[];
DisB_AllTrials_post=[];

Speed_AllTrials_pre=[];
TBF_AllTrials_post=[];

TBF_AllTrials_pre=[];
TBF_AllTrials_post=[];

IBI_AllTrials_pre=[];
IBI_AllTrials_post=[];

nBout_AllTrials_pre=[];
nBout_AllTrials_post=[];

%TBF_AllTrials=[];


for ii=1:NFolder
    
    % Selection manually!
       %[fileName,pathName] = uigetfile('.mat');
    
%Selection automatique dans l ordre alphabetique
    
    fileName =list(ii).name ;
    pathName=strcat(list(ii).folder,'/');
    
    load(strcat(pathName,fileName));
    
output = struct('Fish', [], 'NumberFish',[], 'Genotypes', [], 'geno_index', [], 'G0', [], 'G1', [], ...
    'G2', [], 'Small', [], 'EscapeWindow', [], 'IBI_pre', [], 'IBI_post', [], 'IBI', [], 'Well_ID', [], ...
    'Speed', [], 'index', [], 'allindex', [], 'TailAngle', [], 'Distance', [], 'TBF', [], 'TimeBout_pre', [], ...
    'TimeBout_post', [], 'NumberOfOscillations',[], 'Bend_Amplitude', [], 'Bend_Timing',[],'results',[],...
    'BoutDuration_pre',[], 'BoutDuration_post',[], 'NumberOfOscillations_pre',[], 'NumberOfOscillations_post',[], ...
     'BoutDistance_pre',[],'BoutDistance_post',[],'Speed_pre',[],'Speed_post',[],'nBout_pre', [], 'nBout_post',[],...
     'FishGeno',[],'Fish_G2',[], 'Fish_G1',[],'Fish_G0',[],'TBF_pre',[],'TBF_post',[]);

close all;


if contains(fileName,'20180926')
    EscapeWindow = [105110 105310];
elseif contains(fileName,'20190111')
    EscapeWindow = [30133 30275];
elseif contains(fileName,'20190213')
    EscapeWindow = [30037 30286];
else
    EscapeWindow = [30082 30282];
end;


Fish = unique([datasetPerFish(:).Condition]);
NumberFish=length(Fish);
Genotypes = unique([datasetPerBout(:).Genotype]);


% calculate for each fish parameters extracted from the structure
% datasetPerBout

 
for i=1:NumberFish;
    
    i
    
    index{Fish(i)}= find(~([datasetPerBout(:).Condition]-Fish(i)));
    allindex{Fish(i)}= find(~([datasetPerBout(:).Condition]-Fish(i)));
    
    geno_index{Fish(i)} = unique([datasetPerBout([index{Fish(i)}]).Genotype]);

    %index{Fish(i)}(1)=[];
    
    
    %calculate Pre_Escape
    
    IBI_pre{Fish(i)} = [datasetPerBout(index{Fish(i)}(find([datasetPerBout(index{Fish(i)}).BoutStart]< EscapeWindow(1)))).InstantaneousIBI];
    TimeBout_pre{Fish(i)} = [datasetPerBout(index{Fish(i)}(find([datasetPerBout(index{Fish(i)}).BoutStart]<EscapeWindow(1)))).BoutStart];
    
    BoutDuration_pre{Fish(i)}=[datasetPerBout(allindex{Fish(i)}(find([datasetPerBout(allindex{Fish(i)}).BoutStart]< EscapeWindow(1)))).BoutDuration];
    NumberOfOscillations_pre{Fish(i)}=[datasetPerBout(allindex{Fish(i)}(find([datasetPerBout(allindex{Fish(i)}).BoutStart]< EscapeWindow(1)))).NumberOfOscillations];
    BoutDistance_pre{Fish(i)}=[datasetPerBout(allindex{Fish(i)}(find([datasetPerBout(allindex{Fish(i)}).BoutStart]< EscapeWindow(1)))).TotalDistance];
    Speed_pre{Fish(i)}=[datasetPerBout(allindex{Fish(i)}(find([datasetPerBout(allindex{Fish(i)}).BoutStart]< EscapeWindow(1)))).Speed];
    
    nBout_pre{Fish(i)}= length(allindex{Fish(i)}(find([datasetPerBout(allindex{Fish(i)}).BoutStart]< EscapeWindow(1))));
     
     %calculate Post_Escape
    IBI_post{Fish(i)}= [datasetPerBout(index{Fish(i)}(find([datasetPerBout(index{Fish(i)}).BoutStart]> EscapeWindow(2)))).InstantaneousIBI];
    TimeBout_post{Fish(i)} = [datasetPerBout(index{Fish(i)}(find([datasetPerBout(index{Fish(i)}).BoutStart]>EscapeWindow(2)))).BoutStart];
    
    BoutDuration_post{Fish(i)}=[datasetPerBout(allindex{Fish(i)}(find([datasetPerBout(allindex{Fish(i)}).BoutStart]>EscapeWindow(2)))).BoutDuration];
    NumberOfOscillations_post{Fish(i)}=[datasetPerBout(allindex{Fish(i)}(find([datasetPerBout(allindex{Fish(i)}).BoutStart]>EscapeWindow(2)))).NumberOfOscillations];
    BoutDistance_post{Fish(i)}=[datasetPerBout(allindex{Fish(i)}(find([datasetPerBout(allindex{Fish(i)}).BoutStart]>EscapeWindow(2)))).TotalDistance];
    Speed_post{Fish(i)}=[datasetPerBout(allindex{Fish(i)}(find([datasetPerBout(allindex{Fish(i)}).BoutStart]>EscapeWindow(2)))).Speed];
    
    nBout_post{Fish(i)}=length( allindex{Fish(i)}( find( [datasetPerBout(allindex{Fish(i)}).BoutStart]>EscapeWindow(2) ) ) );
    
    
        %calculate parameters for all swims
    IBI{Fish(i)} = [datasetPerBout(index{Fish(i)}).InstantaneousIBI];
    TimeBoutIBI{Fish(i)} = [datasetPerBout(index{Fish(i)}).BoutStart];

    TimeBout{Fish(i)} =[datasetPerBout(allindex{Fish(i)}).BoutStart];
    BoutDuration{Fish(i)}=[datasetPerBout(allindex{Fish(i)}).BoutDuration];
    Distance{Fish(i)}=[datasetPerBout(allindex{Fish(i)}).TotalDistance];
    Speed{Fish(i)}=[datasetPerBout(allindex{Fish(i)}).Speed];
    NumberOfOscillations{Fish(i)}=[datasetPerBout(allindex{Fish(i)}).NumberOfOscillations];
    nBout{Fish(i)}=(nBout_pre{Fish(i)})+(nBout_post{Fish(i)});
    
    %Find genotype position
    G0=find(~[geno_index{:}]);
    G1=find([geno_index{:}]==1);
    G2=find([geno_index{:}]==2);

    % now calculate parameters for each bout 
    
    for h=1:length([allindex{Fish(i)}]);
        
        display(['currently processing fish ' num2str(i)])
        display(['currently processing bout number ' num2str(h)])
        
        % need to convert in degrees... otherwise : TailAngle{Fish(i)}{h}=[datasetPerBout(allindex{Fish(i)}(h)).TailAngle_smoothed]';
        
        TailAngle{Fish(i)}{h}=57.2958*[datasetPerBout(allindex{Fish(i)}(h)).TailAngle_smoothed];
        
        Bend_Timing{Fish(i)}{h}=[datasetPerBout(allindex{Fish(i)}(h)).Bend_Timing];
        
        Bend_Amplitude{Fish(i)}{h} = 57.2958*datasetPerBout(allindex{Fish(i)}(h)).Bend_Amplitude;
     
        TBF{Fish(i)}{h} = [datasetPerBout(allindex{Fish(i)}(h)).InstantaneousTBF];
   TBF_pre{Fish(i)}{h}=[datasetPerBout(allindex{Fish(i)}(find([datasetPerBout(allindex{Fish(i)}(h)).BoutStart]< EscapeWindow(1)))).InstantaneousTBF];
   TBF_post{Fish(i)}{h}=[datasetPerBout(allindex{Fish(i)}(find([datasetPerBout(allindex{Fish(i)}(h)).BoutStart]>EscapeWindow(2)))).InstantaneousTBF];
        
    end;
end;


%% select good swimers

Fish_temp=Fish;
% for i=1:NumberFish;
%     
%         if isempty([IBI_pre{Fish(i)}]);
%             disp(['empty cell array IBI pre for fish' num2str(Fish(i))]);
%            
%             
%         elseif mean([IBI_pre{Fish(i)}]) > 3.3 ;   
%                 disp(['Bad IBI for fish' num2str(Fish(i))]);
%               
%                 Fish_temp(i)=NaN;
%                 index{Fish(i)}=[];
%                 allindex{Fish(i)}=[];
%                 
%             else disp(['all is good with'  num2str(Fish(i))]);
% 
%             end;
%         end;
% 
% %return
% 
% Fish_temp(isnan(Fish_temp))=[];

FishGeno=([datasetPerFish.Genotype]);

Fish_G2=Fish_temp(find(FishGeno( find( Fish_temp ) )==2));
Fish_G1=Fish_temp(find(FishGeno( find( Fish_temp ) )==1));
Fish_G0=Fish_temp(find(FishGeno( find( Fish_temp ) )==0));

%% collecte data for all trials

% 1)BoutDuration
BD_AllTrials_pre{1,1+Compt_NTrial}=BoutDuration_pre;
BD_AllTrials_post{1,1+Compt_NTrial}=BoutDuration_post;

% 2)Num of Oscillation
Osc_AllTrials_pre{1,1+Compt_NTrial}=NumberOfOscillations_pre;
Osc_AllTrials_post{1,1+Compt_NTrial}=NumberOfOscillations_post;

% 3)BoutDistance
DisB_AllTrials_pre{1,1+Compt_NTrial}=BoutDistance_pre;
DisB_AllTrials_post{1,1+Compt_NTrial}=BoutDistance_post;

% 4)Speed
Speed_AllTrials_pre{1,1+Compt_NTrial}=Speed_pre;
Speed_AllTrials_post{1,1+Compt_NTrial}=Speed_post;

% 5)TBF
TBF_AllTrials_pre{1,1+Compt_NTrial}=TBF_pre;
TBF_AllTrials_post{1,1+Compt_NTrial}=TBF_post;

% 6)IBI
IBI_AllTrials_pre{1,1+Compt_NTrial}=IBI_pre;
IBI_AllTrials_post{1,1+Compt_NTrial}=IBI_post;

% 7)Num of Bout
nBout_AllTrials_pre{1,1+Compt_NTrial}=nBout_pre;
nBout_AllTrials_post{1,1+Compt_NTrial}=nBout_post;

TBF_AllTrials{1,1+Compt_NTrial}=TBF;

Compt_NTrial=Compt_NTrial+1;
% 
% % output data
% 
% output.index=index;
% output.allindex=allindex;
% 
% output.G0=G0;
% output.G1=G1;
% output.G2=G2;
% output.Genotypes=Genotypes;
% output.geno_index=geno_index;
% 
% output.Fish=Fish;
% output.EscapeWindow=EscapeWindow;
% output.NumberFish=NumberFish;
% 
% output.IBI=IBI;
% output.TimeBoutIBI=TimeBoutIBI;
% output.TimeBout=TimeBout;
% output.BoutDuration=BoutDuration; 
% output.TailAngle=TailAngle;
% output.Bend_Timing=Bend_Timing;
% output.Bend_Amplitude=Bend_Amplitude;
% output.Speed=Speed;
% output.NumberOfOscillations=NumberOfOscillations;
% 
% 
% output.IBI_pre=IBI_pre;
% output.TimeBout_pre=TimeBout_pre;
% output.BoutDuration_pre=BoutDuration_pre;
% output.NumberOfOscillations_pre=NumberOfOscillations_pre;
% output.Speed_pre=Speed_pre;
% output.BoutDistance_pre=BoutDistance_pre;
% output.nBout_pre=nBout_pre;
% 
% 
% output.IBI_post=IBI_post;
% output.TimeBout_post=TimeBout_post;
% output.BoutDuration_post=BoutDuration_post;
% output.NumberOfOscillations_post=NumberOfOscillations_post;
% output.Speed_post=Speed_post;
% output.BoutDistance_post=BoutDistance_post;
% output.nBout_post=nBout_post;
% 
% output.FishGeno=FishGeno;
% output.Fish_G0=Fish_G0;
% output.Fish_G1=Fish_G1;
% output.Fish_G2=Fish_G2;

%return
%% mean Calculation

%WT
for l=1:length(Fish_G2);
    l
   
    % BoutDuration
    G2medianBoutDuration_pre{ii}(l)=median(BD_AllTrials_pre{ii}{(Fish_G2(l))},'omitnan');
    G2medianBoutDuration_post{ii}(l)=median(BD_AllTrials_post{ii}{(Fish_G2(l))},'omitnan');
    
   
    % Number of Oscillation
    G2medianNumOfOsc_pre{ii}(l)=median(Osc_AllTrials_pre{ii}{(Fish_G2(l))},'omitnan');
    G2medianNumOfOsc_post{ii}(l)=median(Osc_AllTrials_post{ii}{(Fish_G2(l))},'omitnan');
    
    % BoutDistance
    
    G2medianBoutDistance_pre{ii}(l)=median(DisB_AllTrials_pre{ii}{(Fish_G2(l))},'omitnan');
    G2medianBoutDistance_post{ii}(l)=median(DisB_AllTrials_post{ii}{(Fish_G2(l))},'omitnan');
    
    % Speed
    G2medianSpeed_pre{ii}(l)=median(Speed_AllTrials_pre{ii}{(Fish_G2(l))},'omitnan');
    G2medianSpeed_post{ii}(l)=median(Speed_AllTrials_post{ii}{(Fish_G2(l))},'omitnan');
    
    % TBF
    G2medianTBF_pre{ii}(l)=median(cell2mat(TBF_AllTrials_pre{ii}{(Fish_G2(l))}),'omitnan');
    G2medianTBF_post{ii}(l)=median(cell2mat(TBF_AllTrials_post{ii}{(Fish_G2(l))}),'omitnan');
    
    
    % IBI
    G2medianIBI_pre{ii}(l)=log(median(IBI_AllTrials_pre{ii}{(Fish_G2(l))},'omitnan'));
    G2medianIBI_post{ii}(l)=log(median(IBI_AllTrials_post{ii}{(Fish_G2(l))},'omitnan'));
    
     
    % nBout
    G2median_nBout_pre{ii}(l)=nBout_AllTrials_pre{ii}{(Fish_G2(l))};
    G2median_nBout_post{ii}(l)=nBout_AllTrials_post{ii}{(Fish_G2(l))};
 
end;
 
 
%Heterozygote
for l=1:length(Fish_G1)
    l
    % Calculate medianBoutDuration
    G1medianBoutDuration_pre{ii}(l)=median(BD_AllTrials_pre{ii}{(Fish_G1(l))},'omitnan');
    G1medianBoutDuration_post{ii}(l)=median(BD_AllTrials_post{ii}{(Fish_G1(l))},'omitnan');
    
   
    % Calculate medianNumber of Oscillation
    G1medianNumOfOsc_pre{ii}(l)=median(Osc_AllTrials_pre{ii}{(Fish_G1(l))},'omitnan');
    G1medianNumOfOsc_post{ii}(l)=median(Osc_AllTrials_post{ii}{(Fish_G1(l))},'omitnan');
    
    % Calculate medianBoutDistance
    
    G1medianBoutDistance_pre{ii}(l)=median(DisB_AllTrials_pre{ii}{(Fish_G1(l))},'omitnan');
    G1medianBoutDistance_post{ii}(l)=median(DisB_AllTrials_post{ii}{(Fish_G1(l))},'omitnan');
    
    % Calculate medianSpeed
    G1medianSpeed_pre{ii}(l)=median(Speed_AllTrials_pre{ii}{(Fish_G1(l))},'omitnan');
    G1medianSpeed_post{ii}(l)=median(Speed_AllTrials_post{ii}{(Fish_G1(l))},'omitnan');
 
        % TBF
    G1medianTBF_pre{ii}(l)=median(cell2mat(TBF_AllTrials_pre{ii}{(Fish_G1(l))}),'omitnan');
    G1medianTBF_post{ii}(l)=median(cell2mat(TBF_AllTrials_post{ii}{(Fish_G1(l))}),'omitnan');
    
    % Calculate median IBI
    G1medianIBI_pre{ii}(l)=log(median(IBI_AllTrials_pre{ii}{(Fish_G1(l))},'omitnan'));
    G1medianIBI_post{ii}(l)=log(median(IBI_AllTrials_post{ii}{(Fish_G1(l))},'omitnan'));

     
    % Calculate median nBout
    G1median_nBout_pre{ii}(l)=nBout_AllTrials_pre{ii}{(Fish_G1(l))};
    G1median_nBout_post{ii}(l)=nBout_AllTrials_post{ii}{(Fish_G1(l))};
 
   
end;
 
 
%Homozygote
for l=1:length(Fish_G0)
    l
    % Calculate medianBoutDuration
    G0medianBoutDuration_pre{ii}(l)=median(BD_AllTrials_pre{ii}{(Fish_G0(l))},'omitnan');
    G0medianBoutDuration_post{ii}(l)=median(BD_AllTrials_post{ii}{(Fish_G0(l))},'omitnan');
    
   
    % Calculate medianNumber of Oscillation
    G0medianNumOfOsc_pre{ii}(l)=median(Osc_AllTrials_pre{ii}{(Fish_G0(l))},'omitnan');
    G0medianNumOfOsc_post{ii}(l)=median(Osc_AllTrials_post{ii}{(Fish_G0(l))},'omitnan');
    
    % Calculate medianBoutDistance
    
    G0medianBoutDistance_pre{ii}(l)=median(DisB_AllTrials_pre{ii}{(Fish_G0(l))},'omitnan');
    G0medianBoutDistance_post{ii}(l)=median(DisB_AllTrials_post{ii}{(Fish_G0(l))},'omitnan');
    
    % Calculate medianSpeed
    G0medianSpeed_pre{ii}(l)=median(Speed_AllTrials_pre{ii}{(Fish_G0(l))},'omitnan');
    G0medianSpeed_post{ii}(l)=median(Speed_AllTrials_post{ii}{(Fish_G0(l))},'omitnan');
 
        % Calculate medianTBF
    G0medianTBF_pre{ii}(l)=median(cell2mat(TBF_AllTrials_pre{ii}{(Fish_G0(l))}),'omitnan');
    G0medianTBF_post{ii}(l)=median(cell2mat(TBF_AllTrials_post{ii}{(Fish_G0(l))}),'omitnan');
    
    % Calculate median IBI
    G0medianIBI_pre{ii}(l)=log(median(IBI_AllTrials_pre{ii}{(Fish_G0(l))},'omitnan'));
    G0medianIBI_post{ii}(l)=log(median(IBI_AllTrials_post{ii}{(Fish_G0(l))},'omitnan'));

     
     % Calculate median nBout
    G0median_nBout_pre{ii}(l)=nBout_AllTrials_pre{ii}{(Fish_G0(l))};
    G0median_nBout_post{ii}(l)=nBout_AllTrials_post{ii}{(Fish_G0(l))};
 
% 
end;



end;

%% subplot pre_Stimulus
h11=figure(11); 
title('Pre_Stimulus');hold on;

subplot(1,5,1)
title('Bout Duration (sec)');hold on;
boxplot_2groupes(G2medianBoutDuration_pre,G0medianBoutDuration_pre);hold on;
ylim([0.1 0.6]);hold off;

subplot(1,5,2)
title('Oscillations (n)'); hold on;
boxplot_2groupes(G2medianNumOfOsc_pre,G0medianNumOfOsc_pre);hold on;
ylim([0 7]);hold off;

subplot(1,5,3)
title('Bout Distance (mm)');hold on;
boxplot_2groupes(G2medianBoutDistance_pre,G0medianBoutDistance_pre);hold on;
ylim([0 2]);hold off;

subplot(1,5,4)
title('Bout Speed (mm/sec)');hold on;
boxplot_2groupes(G2medianSpeed_pre,G0medianSpeed_pre);hold on;
ylim([0 5]);
hold off;

subplot(1,5,5)
title('IBI (sec)');hold on;  
boxplot_2groupes(G2medianIBI_pre,G0medianIBI_pre); hold on;
ylim([-2 5]);
hold off;

%% Post_Stimulus subplot 
h11=figure(11); 
title('Post_Stimulus');hold on;

subplot(1,5,1)
title('Bout Duration (sec)');hold on;
boxplot_2groupes(G2medianBoutDuration_post,G0medianBoutDuration_post);hold on;
ylim([0.1 0.6]);hold off;

subplot(1,5,2)
title('Oscillations (n)'); hold on;
boxplot_2groupes(G2medianNumOfOsc_post,G0medianNumOfOsc_post);hold on;
ylim([0 7]);hold off;

subplot(1,5,3)
title('Bout Distance (mm)');hold on;
boxplot_2groupes(G2medianBoutDistance_post,G0medianBoutDistance_post);hold on;
ylim([0 2]);hold off;

subplot(1,5,4)
title('Bout Speed (mm/sec)');hold on;
boxplot_2groupes(G2medianSpeed_post,G0medianSpeed_post);hold on;
ylim([0 5]);
hold off;

subplot(1,5,5)
title('IBI (sec)');hold on;   
boxplot_2groupes(G2medianIBI_post,G0medianIBI_post); hold on;
ylim([-2 5]);
hold off;

%% PLot BoutDuration

h1=figure(1); % Compare bout durations per conditions (genotypes Homo, Het, Wt)
title('Bout Duration (sec)');hold on;
boxplot_4groupes(G2medianBoutDuration_pre,G0medianBoutDuration_pre,G2medianBoutDuration_post,G0medianBoutDuration_post);hold on;
ylim([0.1 0.6]);
hold off;

% [h,p,ci,stats]=ttest2(cell2mat(G2medianBoutDuration_pre),cell2mat(G0medianBoutDuration_pre));
% [h,p,ci,stats]=ttest2(cell2mat(G2medianBoutDuration_post),cell2mat(G0medianBoutDuration_post));

% saveas(h1,['Bout Duration.fig'])
% saveas(h1,['Bout Duration.png'])
%% Number Of Oscillations
h2=figure(2); % Compare Number of Oscillation per conditions (genotypes Homo, Het, Wt)
title('Number Of Oscillations'); hold on;
boxplot_4groupes(G2medianNumOfOsc_pre,G0medianNumOfOsc_pre,G2medianNumOfOsc_post,G0medianNumOfOsc_post);hold on;

ylim([0 7]);
hold off;

% [h,p,ci,stats]=ttest2(cell2mat(G2medianNumOfOsc_pre),cell2mat(G0medianNumOfOsc_pre));
% [h,p,ci,stats]=ttest2(cell2mat(G2medianNumOfOsc_post),cell2mat(G0medianNumOfOsc_post));


% plot_WT_scatter(G2medianNumOfOsc_pre,G2medianNumOfOsc_post);hold on;
% plot_Het_scatter(G1medianNumOfOsc_pre,G1medianNumOfOsc_post);hold on;
% plot_Homo_scatter(G0medianNumOfOsc_pre,G0medianNumOfOsc_post);hold off;
% saveas(h2,['Number of Oscillations.fig'])
% saveas(h2,['Number of Oscillations.png'])  
%% Bout Distance
 
h3=figure(3); % Compare Bout Distance per conditions (genotypes Homo, Het, Wt)
title('Bout Distance (mm)');hold on;
boxplot_4groupes(G2medianBoutDistance_pre,G0medianBoutDistance_pre,G2medianBoutDistance_post,G0medianBoutDistance_post);hold on;

ylim([0 2]);
hold off;

% [h,p,ci,stats]=ttest2(cell2mat(G2medianBoutDistance_pre),cell2mat(G0medianBoutDistance_pre));
% [h,p,ci,stats]=ttest2(cell2mat(G2medianBoutDistance_post),cell2mat(G0medianBoutDistance_post));


 plot_WT_scatter(G2medianBoutDistance_pre,G2medianBoutDistance_post);
% plot_Het_scatter(G1medianBoutDistance_pre,G1medianBoutDistance_post);
 plot_Homo_scatter(G0medianBoutDistance_pre,G0medianBoutDistance_post);

%saveas(h3,['BoutDistance.fig'])
% saveas(h3,['BoutDistance.png'])  
%% Speed 
 
h4=figure(4); % Compare Speed per conditions (genotypes Homo, Het, Wt)
title('Speed (mm/s)');hold on;
boxplot_4groupes(G2medianSpeed_pre,G0medianSpeed_pre,G2medianSpeed_post,G0medianSpeed_post);hold on;
ylim([0 5]);
hold off;

% [h,p,ci,stats]=ttest2(cell2mat(G2medianSpeed_pre),cell2mat(G0medianSpeed_pre));
% [h,p,ci,stats]=ttest2(cell2mat(G2medianSpeed_post),cell2mat(G0medianSpeed_post));


% plot_WT_scatter(G2medianSpeed_pre,G2medianSpeed_post);hold on;
% plot_Het_scatter(G1medianSpeed_pre,G1medianSpeed_post);hold on;
% plot_Homo_scatter(G0medianSpeed_pre,G0medianSpeed_post);hold off;
% saveas(h4,['Speed.fig'])
% saveas(h4,['Speed.png'])  
%% TBF 
 
h5=figure(5); % Compare TBF per conditions (genotypes Homo, Het, Wt)
title('TBF');hold on;
boxplot_4groupes(G2medianTBF_pre,G0medianTBF_pre,G2medianTBF_post,G0medianTBF_post);hold on;
%ylim([0 5]);
hold off;

% [h,p,ci,stats]=ttest2(cell2mat(G2medianTBF_pre),cell2mat(G0medianTBF_pre));
% [h,p,ci,stats]=ttest2(cell2mat(G2medianTBF_post),cell2mat(G0medianTBF_post));


% plot_WT_scatter(G2medianTBF_pre,G2medianTBF_post);hold on;
% plot_Het_scatter(G1medianTBF_pre,G1medianTBF_post);hold on;
% plot_Homo_scatter(G0medianTBF_pre,G0medianTBF_post);hold off; 
%saveas(h5,['TBF.fig'])
%saveas(h5,['TBF.png'])
%% PLot IBI
 
h6=figure(6); % Compare bout durations per conditions (genotypes Homo, Het, Wt)
title('IBI (sec)');   

boxplot_4groupes(G2medianIBI_pre,G0medianIBI_pre,G2medianIBI_post,G0medianIBI_post); hold on;
%ylim([0 5]);
hold off;

% [h,p,ci,stats]=ttest2(cell2mat(G2medianIBI_pre),cell2mat(G0medianIBI_pre));
% [h,p,ci,stats]=ttest2(cell2mat(G2medianIBI_post),cell2mat(G0medianIBI_post));


% plot_WT_scatter(G2medianIBI_pre,G2medianIBI_post);hold on;
% plot_Het_scatter(G1medianIBI_pre,G1medianIBI_post);hold on;
% plot_Homo_scatter(G0medianIBI_pre,G0medianIBI_post);hold off;
% saveas(h6,['IBI.fig'])
% saveas(h6,['IBI.png'])
%% nBout
 
% h7=figure(7); % Compare _nBout per conditions (genotypes Homo, Het, Wt)
% title('nBout'); 
% boxplot_4groupes(G2median_nBout_pre,G0median_nBout_pre,G2median_nBout_post,G0median_nBout_post); hold on;
% hold off;

% plot_WT_scatter(G2median_nBout_pre,G2median_nBout_post);hold on;
% %plot_Het_scatter(G1median_nBout_pre,G1median_nBout_post);hold on;
% plot_Homo_scatter(G0median_nBout_pre,G0median_nBout_post);hold off; 
% saveas(h7,['nBout.fig'])
% saveas(h7,['nBout.png'])  

