 %% Clean workspace and load data
close all
clear variables
clc

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
Speed_AllTrials_post=[];

TBF_AllTrials_pre=[];
TBF_AllTrials_post=[];

IBI_AllTrials_pre=[];
IBI_AllTrials_post=[];

nBout_AllTrials_pre=[];
nBout_AllTrials_post=[];




for ii=1:NFolder
    
    % Selection manually!
       %[fileName,pathName] = uigetfile('.mat');
    
%Selection automatique dans l ordre alphabetique
    
    fileName =list(ii).name ;
    pathName=strcat(list(ii).folder,'/');
    
    load(strcat(pathName,fileName));
    
output = struct( 'NTrial', [], 'Fish_ID', [], 'FishGeno',[], ...
        'mIBI_pre', [], 'mBoutDuration_pre',[],'mNumberOfOscillations_pre',[],'mBoutDistance_pre',[],'mSpeed_pre',[],'mnBout_pre', [], 'mTBF_pre',[],...
    'mIBI_post', [], 'mBoutDuration_post',[],'mNumberOfOscillations_post',[],'mBoutDistance_post',[],'mSpeed_post',[], 'mnBout_post',[], 'mTBF_post',[]);
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

%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
% Define genotype

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
    

    %calculate Pre_Escape
    IBI_pre{Fish(i)} = [datasetPerBout(BoutsIBI_Pre_Escape).InstantaneousIBI];
    TimeBout_pre{Fish(i)} = [datasetPerBout(BoutsIBI_Pre_Escape).BoutStart];
    
    BoutDuration_pre{Fish(i)}=[datasetPerBout(Bouts_Pre_Escape).BoutDuration];
    NumberOfOscillations_pre{Fish(i)}=[datasetPerBout(Bouts_Pre_Escape).NumberOfOscillations];
    BoutDistance_pre{Fish(i)}=[datasetPerBout(Bouts_Pre_Escape).TotalDistance];
    Speed_pre{Fish(i)}=[datasetPerBout(Bouts_Pre_Escape).Speed];
    nBout_pre{Fish(i)}= length(Bouts_Pre_Escape);
    TBF_pre{Fish(i)}=[datasetPerBout(Bouts_Pre_Escape).InstantaneousTBF];


     %calculate Post_Escape
     
    IBI_post{Fish(i)} = [datasetPerBout(BoutsIBI_Post_Escape).InstantaneousIBI];
    TimeBout_post{Fish(i)} = [datasetPerBout(BoutsIBI_Post_Escape).BoutStart];
    
    BoutDuration_post{Fish(i)}=[datasetPerBout(Bouts_Post_Escape).BoutDuration];
    NumberOfOscillations_post{Fish(i)}=[datasetPerBout(Bouts_Post_Escape).NumberOfOscillations];
    BoutDistance_post{Fish(i)}=[datasetPerBout(Bouts_Post_Escape).TotalDistance];
    Speed_post{Fish(i)}=[datasetPerBout(Bouts_Post_Escape).Speed];
    nBout_post{Fish(i)}=length(Bouts_Post_Escape);
    TBF_post{Fish(i)}=[datasetPerBout(Bouts_Post_Escape).InstantaneousTBF];
 
            
    % calculate parameters for all swims
      
      IBI{Fish(i)} = [datasetPerBout(index{Fish(i)}).InstantaneousIBI];
      TimeBoutIBI{Fish(i)} = [datasetPerBout(index{Fish(i)}).BoutStart];
      
      TimeBout{Fish(i)} =[datasetPerBout(allindex{Fish(i)}).BoutStart];
      nBout{Fish(i)}=(nBout_pre{Fish(i)})+(nBout_post{Fish(i)});
      
%     BoutDuration{Fish(i)}=[datasetPerBout(allindex{Fish(i)}).BoutDuration];
%     Distance{Fish(i)}=[datasetPerBout(allindex{Fish(i)}).TotalDistance];
%     Speed{Fish(i)}=[datasetPerBout(allindex{Fish(i)}).Speed];
%     NumberOfOscillations{Fish(i)}=[datasetPerBout(allindex{Fish(i)}).NumberOfOscillations];
   

    % now calculate parameters for each bout 
    
    for h=1:length(idx{Fish(i)});
        
        display(['currently processing fish ' num2str(i)])
        display(['currently processing bout number ' num2str(h)])
        
        % need to convert in degrees... otherwise : TailAngle{Fish(i)}{h}=[datasetPerBout(allindex{Fish(i)}(h)).TailAngle_smoothed]';
        
        TailAngle{Fish(i)}{h}=57.2958*[datasetPerBout(allindex{Fish(i)}(h)).TailAngle_smoothed]';
        
        Bend_Timing{Fish(i)}{h}=[datasetPerBout(allindex{Fish(i)}(h)).Bend_Timing];
        
        Bend_Amplitude{Fish(i)}{h} = 57.2958*datasetPerBout(allindex{Fish(i)}(h)).Bend_Amplitude;
     
        TBF{Fish(i)}{h} = [datasetPerBout(allindex{Fish(i)}(h)).InstantaneousTBF];
          
        
    end;
    
    output(i).NTrial=NTrial(i); 
    output(i).Fish_ID=Fish_ID(i);
    output(i).FishGeno=FishGeno(i); 
    
    output(i).mIBI_pre=mean(IBI_pre{Fish(i)},'omitnan');
    output(i).mBoutDuration_pre=mean(BoutDuration_pre{Fish(i)},'omitnan');
    output(i).mNumberOfOscillations_pre=mean(NumberOfOscillations_pre{Fish(i)},'omitnan');
    output(i).mSpeed_pre=mean(Speed_pre{Fish(i)},'omitnan');
    output(i).mBoutDistance_pre=mean(BoutDistance_pre{Fish(i)},'omitnan');
    output(i).mnBout_pre=mean(nBout_pre{Fish(i)},'omitnan');
    output(i).mTBF_pre=mean(TBF_pre{Fish(i)},'omitnan');
    
    output(i).mIBI_post=mean(IBI_post{Fish(i)},'omitnan');
    output(i).mBoutDuration_post=mean(BoutDuration_post{Fish(i)},'omitnan');
    output(i).mNumberOfOscillations_post=mean(NumberOfOscillations_post{Fish(i)},'omitnan');
    output(i).mSpeed_post=mean(Speed_post{Fish(i)},'omitnan');
    output(i).mBoutDistance_post=mean(BoutDistance_post{Fish(i)},'omitnan');
    output(i).mnBout_post=mean(nBout_post{Fish(i)},'omitnan');
    output(i).mTBF_post=mean(TBF_post{Fish(i)},'omitnan');
    
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

TBF_AllTrials_pre{1,1+Compt_NTrial}=TBF_pre; 
TBF_AllTrials_post{1,1+Compt_NTrial}=TBF_post;



Compt_NTrial=Compt_NTrial+1;
 
%% output data
% output.Fish=Fish;
% output.EscapeWindow=EscapeWindow;
% output.NumberFish=NumberFish;
% output.index=index;
% output.allindex=allindex;
% 
% output.FishGeno=FishGeno;
% output.Fish_G0=Fish_G0;
% output.Fish_G1=Fish_G1;
% output.Fish_G2=Fish_G2;

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
 

% output.IBI=IBI;
% output.TimeBoutIBI=TimeBoutIBI;
% output.TimeBout=TimeBout;
% output.BoutDuration=BoutDuration; 
% output.TailAngle=TailAngle;
% output.Bend_Timing=Bend_Timing;
% output.Bend_Amplitude=Bend_Amplitude;
% output.Speed=Speed;
% output.NumberOfOscillations=NumberOfOscillations;

%return
%% mean Calculation
G2medianBoutDuration_pre{ii}(l)=[];
for l=1:length(Fish_G2);
    l
   
    % Calculate medianBoutDuration
    G2medianBoutDuration_pre{ii}(l)=median(BD_AllTrials_pre{ii}{(Fish_G2(l))});
    G2medianBoutDuration_post{ii}(l)=median(BD_AllTrials_post{ii}{(Fish_G2(l))});
    
   
    % Calculate medianNumber of Oscillation
    G2medianNumOfOsc_pre{ii}(l)=median(Osc_AllTrials_pre{ii}{(Fish_G2(l))});
    G2medianNumOfOsc_post{ii}(l)=median(Osc_AllTrials_post{ii}{(Fish_G2(l))});
    
    % Calculate medianBoutDistance
    
    G2medianBoutDistance_pre{ii}(l)=median(DisB_AllTrials_pre{ii}{(Fish_G2(l))});
    G2medianBoutDistance_post{ii}(l)=median(DisB_AllTrials_post{ii}{(Fish_G2(l))});
    
    % Calculate medianSpeed
    G2medianSpeed_pre{ii}(l)=median(Speed_AllTrials_pre{ii}{(Fish_G2(l))});
    G2medianSpeed_post{ii}(l)=median(Speed_AllTrials_post{ii}{(Fish_G2(l))});
    
    % Calculate median TBF
    
    G2medianTBF_pre{ii}(l)=median(TBF_AllTrials_pre{ii}{(Fish_G2(l))});
    G2medianTBF_post{ii}(l)=median(TBF_AllTrials_post{ii}{(Fish_G2(l))});
    
    
    % Calculate median IBI
    G2medianIBI_pre{ii}(l)=log(median(IBI_AllTrials_pre{ii}{(Fish_G2(l))}));
    G2medianIBI_post{ii}(l)=log(median(IBI_AllTrials_post{ii}{(Fish_G2(l))}));
    
%    % Calculate single IBI
%     G2Single_IBI_pre{ii}(l)= median(IBI_AllTrials_pre{ii}{(Fish_G2(l))}(end)),'omitnan';
%     G2Single_IBI_post{ii}(l)=median(IBI_AllTrials_pre{ii}{(Fish_G2(l))}(1)),'omitnan';
     
    % Calculate median nBout
    G2median_nBout_pre{ii}(l)=nBout_AllTrials_pre{ii}{(Fish_G2(l))};
    G2median_nBout_post{ii}(l)=nBout_AllTrials_post{ii}{(Fish_G2(l))};
 
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
    
 
%    % Calculate single IBI
%     G1Single_IBI_pre{ii}(l)=log(median(IBI_AllTrials_pre{ii}{(Fish_G1(l))}(end)));
%     G1Single_IBI_post{ii}(l)=log(median(IBI_AllTrials_pre{ii}{(Fish_G1(l))}(1)));
    
 
    % Calculate median nBout
    G1median_nBout_pre{ii}(l)=nBout_AllTrials_pre{ii}{(Fish_G1(l))};
    G1median_nBout_post{ii}(l)=nBout_AllTrials_post{ii}{(Fish_G1(l))};
 
   
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
    
    % Calculate medianBoutDistance
    
    G0medianBoutDistance_pre{ii}(l)=median(DisB_AllTrials_pre{ii}{(Fish_G0(l))});
    G0medianBoutDistance_post{ii}(l)=median(DisB_AllTrials_post{ii}{(Fish_G0(l))});
    
    % Calculate medianSpeed
    G0medianSpeed_pre{ii}(l)=median(Speed_AllTrials_pre{ii}{(Fish_G0(l))});
    G0medianSpeed_post{ii}(l)=median(Speed_AllTrials_post{ii}{(Fish_G0(l))});
    
    % Calculate median TBF
    G0medianTBF_pre{ii}(l)=median(TBF_AllTrials_pre{ii}{(Fish_G0(l))});
    G0medianTBF_post{ii}(l)=median(TBF_AllTrials_post{ii}{(Fish_G0(l))});
    
    
    % Calculate median IBI
    G0medianIBI_pre{ii}(l)=log(median(IBI_AllTrials_pre{ii}{(Fish_G0(l))}));
    G0medianIBI_post{ii}(l)=log(median(IBI_AllTrials_post{ii}{(Fish_G0(l))}));
    
%     % Calculate single IBI
%    G0Single_IBI_pre{ii}(l)=log(median(IBI_AllTrials_pre{ii}{(Fish_G0(l))}(end)));
%     G0Single_IBI_post{ii}(l)=log(median(IBI_AllTrials_pre{ii}{(Fish_G0(l))}(1)));
     
%     % Calculate median nBout
    G0median_nBout_pre{ii}(l)=nBout_AllTrials_pre{ii}{(Fish_G0(l))};
    G0median_nBout_post{ii}(l)=nBout_AllTrials_post{ii}{(Fish_G0(l))};
 
end;



ii= ii+1;


end;


%% PLot BoutDuration

% h1=figure(1); % Compare bout durations per conditions (genotypes Homo, Het, Wt)
% title('Bout Duration (sec)');  
% 
% plot_WT_scatter(G2medianBoutDuration_pre,G2medianBoutDuration_post);hold on;
% plot_Het_scatter(G1medianBoutDuration_pre,G1medianBoutDuration_post);hold on;
% plot_Homo_scatter(G0medianBoutDuration_pre,G0medianBoutDuration_post);hold off;
% % saveas(h1,['Bout Duration.fig'])
% % saveas(h1,['Bout Duration.png']) 
% %% Plot Number of Oscillation
%  
% h2=figure(2); % Compare Number of Oscillation per conditions (genotypes Homo, Het, Wt)
% title('Number Of Oscillations'); 
%  
% plot_WT_scatter(G2medianNumOfOsc_pre,G2medianNumOfOsc_post);hold on;
% plot_Het_scatter(G1medianNumOfOsc_pre,G1medianNumOfOsc_post);hold on;
% plot_Homo_scatter(G0medianNumOfOsc_pre,G0medianNumOfOsc_post);hold off;
% % saveas(h2,['Number of Oscillations.fig'])
% % saveas(h2,['Number of Oscillations.png'])  
% %% Bout Distance
%  
% h3=figure(3); % Compare Bout Distance per conditions (genotypes Homo, Het, Wt)
% title('Bout Distance (mm)');  
%  
% plot_WT_scatter(G2medianBoutDistance_pre,G2medianBoutDistance_post);hold on;
% plot_Het_scatter(G1medianBoutDistance_pre,G1medianBoutDistance_post);
% plot_Homo_scatter(G0medianBoutDistance_pre,G0medianBoutDistance_post);hold off;
% %saveas(h3,['BoutDistance.fig'])
% % saveas(h3,['BoutDistance.png'])  
% %% Speed 
%  
% h4=figure(4); % Compare Speed per conditions (genotypes Homo, Het, Wt)
% title('Speed (mm/s)');
%  
% plot_WT_scatter(G2medianSpeed_pre,G2medianSpeed_post);hold on;
% plot_Het_scatter(G1medianSpeed_pre,G1medianSpeed_post);hold on;
% plot_Homo_scatter(G0medianSpeed_pre,G0medianSpeed_post);hold off;
% % saveas(h4,['Speed.fig'])
% % saveas(h4,['Speed.png'])  
% %% TBF 
%  
% h5=figure(5); % Compare TBF per conditions (genotypes Homo, Het, Wt)
% title('TBF'); 
% plot_WT_scatter(G2medianTBF_pre,G2medianTBF_post);hold on;
% plot_Het_scatter(G1medianTBF_pre,G1medianTBF_post);hold on;
% plot_Homo_scatter(G0medianTBF_pre,G0medianTBF_post);hold off; 
% %saveas(h5,['TBF.fig'])
% %saveas(h5,['TBF.png'])
% %% PLot IBI
%  
% h6=figure(6); % Compare bout durations per conditions (genotypes Homo, Het, Wt)
% title('IBI (sec)');    
%  
% plot_WT_scatter(G2medianIBI_pre,G2medianIBI_post);hold on;
% plot_Het_scatter(G1medianIBI_pre,G1medianIBI_post);hold on;
% plot_Homo_scatter(G0medianIBI_pre,G0medianIBI_post);hold off;
% % saveas(h6,['IBI.fig'])
% % saveas(h6,['IBI.png'])
% %% nBout
%  
% h7=figure(7); % Compare _nBout per conditions (genotypes Homo, Het, Wt)
% title('nBout'); 
%  
% plot_WT_scatter(G2median_nBout_pre,G2median_nBout_post);hold on;
% plot_Het_scatter(G1median_nBout_pre,G1median_nBout_post);hold on;
% plot_Homo_scatter(G0median_nBout_pre,G0median_nBout_post);hold off; 
% % saveas(h7,['nBout.fig'])
% % saveas(h7,['nBout.png'])  




