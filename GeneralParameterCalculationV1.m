 %% Clean workspace and load data
 clear 

%% load the dataset file
InitialFolder='/Users/bong-iquan/MATLAB/sst1.1_mutant/Slow_swim';
 
% Figure out how many trials to analyse
cd(InitialFolder)
list=dir('*_workspace_SST_20190220_new2*');
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
     'FishGeno',[],'Fish_G2',[], 'Fish_G1',[],'Fish_G0',[]);

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
    
    for h=1:length([allindex{Fish(i)}]);
        
        display(['currently processing fish ' num2str(i)])
        display(['currently processing bout number ' num2str(h)])
        
        % need to convert in degrees... otherwise : TailAngle{Fish(i)}{h}=[datasetPerBout(allindex{Fish(i)}(h)).TailAngle_smoothed]';
        
        TailAngle{Fish(i)}{h}=57.2958*[datasetPerBout(allindex{Fish(i)}(h)).TailAngle_smoothed]';
        
        Bend_Timing{Fish(i)}{h}=[datasetPerBout(allindex{Fish(i)}(h)).Bend_Timing];
        
        Bend_Amplitude{Fish(i)}{h} = 57.2958*datasetPerBout(allindex{Fish(i)}(h)).Bend_Amplitude;
     
        TBF{Fish(i)}{h} = [datasetPerBout(allindex{Fish(i)}(h)).InstantaneousTBF];
          
        
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

%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
% Define genotype

FishGeno=([datasetPerFish.Genotype]);

Fish_G2=Fish_temp(find(FishGeno( find( Fish_temp ) )==2));
Fish_G1=Fish_temp(find(FishGeno( find( Fish_temp ) )==1));
Fish_G0=Fish_temp(find(FishGeno( find( Fish_temp ) )==0));

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

%WT
for l=1:length(Fish_G2);
    l
   
    % Calculate meanBoutDuration
    G2meanBoutDuration_pre{ii}(l)=mean(BD_AllTrials_pre{ii}{(Fish_G2(l))});
    G2meanBoutDuration_post{ii}(l)=mean(BD_AllTrials_post{ii}{(Fish_G2(l))});
    
   
    % Calculate meanNumber of Oscillation
    G2meanNumOfOsc_pre{ii}(l)=mean(Osc_AllTrials_pre{ii}{(Fish_G2(l))});
    G2meanNumOfOsc_post{ii}(l)=mean(Osc_AllTrials_post{ii}{(Fish_G2(l))});
    
    % Calculate meanBoutDistance
    
    G2meanBoutDistance_pre{ii}(l)=mean(DisB_AllTrials_pre{ii}{(Fish_G2(l))});
    G2meanBoutDistance_post{ii}(l)=mean(DisB_AllTrials_post{ii}{(Fish_G2(l))});
    
    % Calculate meanSpeed
    G2meanSpeed_pre{ii}(l)=mean(Speed_AllTrials_pre{ii}{(Fish_G2(l))});
    G2meanSpeed_post{ii}(l)=mean(Speed_AllTrials_post{ii}{(Fish_G2(l))});
    
    % Calculate mean TBF
    
    G2meanTBF_pre{ii}(l)=mean(TBF_AllTrials_pre{ii}{(Fish_G2(l))});
    G2meanTBF_post{ii}(l)=mean(TBF_AllTrials_post{ii}{(Fish_G2(l))});
    
    
    % Calculate mean IBI
    G2meanIBI_pre{ii}(l)=log(mean(IBI_AllTrials_pre{ii}{(Fish_G2(l))}));
    G2meanIBI_post{ii}(l)=log(mean(IBI_AllTrials_post{ii}{(Fish_G2(l))}));
    
%    % Calculate single IBI
%     G2Single_IBI_pre{ii}(l)= mean(IBI_AllTrials_pre{ii}{(Fish_G2(l))}(end)),'omitnan';
%     G2Single_IBI_post{ii}(l)=mean(IBI_AllTrials_pre{ii}{(Fish_G2(l))}(1)),'omitnan';
     
    % Calculate mean nBout
    G2mean_nBout_pre{ii}(l)=nBout_AllTrials_pre{ii}{(Fish_G2(l))};
    G2mean_nBout_post{ii}(l)=nBout_AllTrials_post{ii}{(Fish_G2(l))};

end;


%Heterozygote
for l=1:length(Fish_G1)
    l
    % Calculate meanBoutDuration
    G1meanBoutDuration_pre{ii}(l)=mean(BD_AllTrials_pre{ii}{(Fish_G1(l))});
    G1meanBoutDuration_post{ii}(l)=mean(BD_AllTrials_post{ii}{(Fish_G1(l))});
    
   
    % Calculate meanNumber of Oscillation
    G1meanNumOfOsc_pre{ii}(l)=mean(Osc_AllTrials_pre{ii}{(Fish_G1(l))});
    G1meanNumOfOsc_post{ii}(l)=mean(Osc_AllTrials_post{ii}{(Fish_G1(l))});
    
    % Calculate meanBoutDistance
    
    G1meanBoutDistance_pre{ii}(l)=mean(DisB_AllTrials_pre{ii}{(Fish_G1(l))});
    G1meanBoutDistance_post{ii}(l)=mean(DisB_AllTrials_post{ii}{(Fish_G1(l))});
    
    % Calculate meanSpeed
    G1meanSpeed_pre{ii}(l)=mean(Speed_AllTrials_pre{ii}{(Fish_G1(l))});
    G1meanSpeed_post{ii}(l)=mean(Speed_AllTrials_post{ii}{(Fish_G1(l))});

    
    % Calculate mean IBI
    G1meanIBI_pre{ii}(l)=log(mean(IBI_AllTrials_pre{ii}{(Fish_G1(l))}));
    G1meanIBI_post{ii}(l)=log(mean(IBI_AllTrials_post{ii}{(Fish_G1(l))}));

    
    % Calculate mean TBF
    G1meanTBF_pre{ii}(l)=mean(TBF_AllTrials_pre{ii}{(Fish_G1(l))});
    G1meanTBF_post{ii}(l)=mean(TBF_AllTrials_post{ii}{(Fish_G1(l))});
    

%    % Calculate single IBI
%     G1Single_IBI_pre{ii}(l)=log(mean(IBI_AllTrials_pre{ii}{(Fish_G1(l))}(end)));
%     G1Single_IBI_post{ii}(l)=log(mean(IBI_AllTrials_pre{ii}{(Fish_G1(l))}(1)));
%     
%     % Calculate mean nBout
    G1mean_nBout_pre{ii}(l)=nBout_AllTrials_pre{ii}{(Fish_G1(l))};
    G1mean_nBout_post{ii}(l)=nBout_AllTrials_post{ii}{(Fish_G1(l))};

   
end;


%Homozygote
for l=1:length(Fish_G0)
    l
    % Calculate meanBoutDuration
    G0meanBoutDuration_pre{ii}(l)=mean(BD_AllTrials_pre{ii}{(Fish_G0(l))});
    G0meanBoutDuration_post{ii}(l)=mean(BD_AllTrials_post{ii}{(Fish_G0(l))});
    
   
    % Calculate meanNumber of Oscillation
    G0meanNumOfOsc_pre{ii}(l)=mean(Osc_AllTrials_pre{ii}{(Fish_G0(l))});
    G0meanNumOfOsc_post{ii}(l)=mean(Osc_AllTrials_post{ii}{(Fish_G0(l))});
    
    % Calculate meanBoutDistance
    
    G0meanBoutDistance_pre{ii}(l)=mean(DisB_AllTrials_pre{ii}{(Fish_G0(l))});
    G0meanBoutDistance_post{ii}(l)=mean(DisB_AllTrials_post{ii}{(Fish_G0(l))});
    
    % Calculate meanSpeed
    G0meanSpeed_pre{ii}(l)=mean(Speed_AllTrials_pre{ii}{(Fish_G0(l))});
    G0meanSpeed_post{ii}(l)=mean(Speed_AllTrials_post{ii}{(Fish_G0(l))});
    
    % Calculate mean TBF
    G0meanTBF_pre{ii}(l)=mean(TBF_AllTrials_pre{ii}{(Fish_G0(l))});
    G0meanTBF_post{ii}(l)=mean(TBF_AllTrials_post{ii}{(Fish_G0(l))});
    
    
    % Calculate mean IBI
    G0meanIBI_pre{ii}(l)=log(mean(IBI_AllTrials_pre{ii}{(Fish_G0(l))}));
    G0meanIBI_post{ii}(l)=log(mean(IBI_AllTrials_post{ii}{(Fish_G0(l))}));
    
%     % Calculate single IBI
%    G0Single_IBI_pre{ii}(l)=log(mean(IBI_AllTrials_pre{ii}{(Fish_G0(l))}(end)));
%     G0Single_IBI_post{ii}(l)=log(mean(IBI_AllTrials_pre{ii}{(Fish_G0(l))}(1)));
     
%     % Calculate mean nBout
    G0mean_nBout_pre{ii}(l)=nBout_AllTrials_pre{ii}{(Fish_G0(l))};
    G0mean_nBout_post{ii}(l)=nBout_AllTrials_post{ii}{(Fish_G0(l))};

end;




end;


%% PLot BoutDuration

h1=figure(1); % Compare bout durations per conditions (genotypes Homo, Het, Wt)
title('Bout Duration (sec)');  

plot_WT_scatter(G2meanBoutDuration_pre,G2meanBoutDuration_post);
hold on;
plot_Het_scatter(G1meanBoutDuration_pre,G1meanBoutDuration_post);
hold on;
plot_Homo_scatter(G0meanBoutDuration_pre,G0meanBoutDuration_post);
hold off;

% saveas(h1,['Bout Duration.fig'])
% saveas(h1,['Bout Duration.png']) 
%% Plot Number of Oscillation

h2=figure(2); % Compare Number of Oscillation per conditions (genotypes Homo, Het, Wt)
title('Number Of Oscillations'); 

plot_WT_scatter(G2meanNumOfOsc_pre,G2meanNumOfOsc_post);
hold on;
plot_Het_scatter(G1meanNumOfOsc_pre,G1meanNumOfOsc_post);
hold on;
plot_Homo_scatter(G0meanNumOfOsc_pre,G0meanNumOfOsc_post);
hold off;

% saveas(h2,['Number of Oscillations.fig'])
% saveas(h2,['Number of Oscillations.png'])  
%% Bout Distance

h3=figure(3); % Compare Bout Distance per conditions (genotypes Homo, Het, Wt)
title('Bout Distance (mm)');  

plot_WT_scatter(G2meanBoutDistance_pre,G2meanBoutDistance_post);
hold on;
plot_Het_scatter(G1meanBoutDistance_pre,G1meanBoutDistance_post);
hold on;
plot_Homo_scatter(G0meanBoutDistance_pre,G0meanBoutDistance_post);
hold off;

%saveas(h3,['BoutDistance.fig'])
% saveas(h3,['BoutDistance.png'])  

%% Speed 

h4=figure(4); % Compare Speed per conditions (genotypes Homo, Het, Wt)
title('Speed (mm/s)');

plot_WT_scatter(G2meanSpeed_pre,G2meanSpeed_post);
hold on;
plot_Het_scatter(G1meanSpeed_pre,G1meanSpeed_post);
hold on;
plot_Homo_scatter(G0meanSpeed_pre,G0meanSpeed_post);
hold off;

% saveas(h4,['Speed.fig'])
% saveas(h4,['Speed.png'])  

%% TBF 
 
h5=figure(5); % Compare TBF per conditions (genotypes Homo, Het, Wt)
title('TBF'); 


plot_WT_scatter(G2meanTBF_pre,G2meanTBF_post);
hold on;
plot_Het_scatter(G1meanTBF_pre,G1meanTBF_post);
hold on;
plot_Homo_scatter(G0meanTBF_pre,G0meanTBF_post);
hold off;
 
%saveas(h5,['TBF.fig'])
%saveas(h5,['TBF.png'])

%% PLot IBI
 
h6=figure(6); % Compare bout durations per conditions (genotypes Homo, Het, Wt)
title('IBI (sec)');    

plot_WT_scatter(G2meanIBI_pre,G2meanIBI_post);
hold on;
plot_Het_scatter(G1meanIBI_pre,G1meanIBI_post);
hold on;
plot_Homo_scatter(G0meanIBI_pre,G0meanIBI_post);
hold off;

% saveas(h6,['IBI.fig'])
% saveas(h6,['IBI.png'])

%% nBout

h7=figure(7); % Compare _nBout per conditions (genotypes Homo, Het, Wt)
title('nBout'); 

plot_WT_scatter(G2mean_nBout_pre,G2mean_nBout_post);
hold on;
plot_Het_scatter(G1mean_nBout_pre,G1mean_nBout_post);
hold on;
plot_Homo_scatter(G0mean_nBout_pre,G0mean_nBout_post);
hold off;
 
% saveas(h7,['nBout.fig'])
% saveas(h7,['nBout.png'])  

%% Single IBI plot
% h4=figure(4); %Plot single IBI before (X) and after (Y) escape
%  
% title('single IBI before (X) and after (Y) escape');    
% 
%   for i=1:NumberFish; 
%         
%         if isempty([IBI_post{Fish(i)}]) | isempty([IBI_pre{Fish(i)}])
%         disp(['empty cell array IBI pre or post for fish' num2str(i)])
%         else 
%         plot(test,test,'k-'); title('single IBI before (X) and after (Y) escape');
%             if [datasetPerBout([index{Fish(i)}]).Genotype]==0;
%                 plot(mean([IBI_pre{Fish(i)}(end)]), mean([IBI_post{Fish(i)}(1)]),'ko-'); hold on;
%             else if [datasetPerBout([index{Fish(i)}]).Genotype]==1;
%                     plot(mean([IBI_pre{Fish(i)}(end)]), mean([IBI_post{Fish(i)}(1)]),'go-'); hold on;
%                 else if [datasetPerBout([index{Fish(i)}]).Genotype]==2;
%                         plot(mean([IBI_pre{Fish(i)}(end)]), mean([IBI_post{Fish(i)}(1)]),'ro-'); hold on;
%                     end;
%                 end;
%             end;
%         end;
%     end;
%     
%   saveas(h4,['single IBI before (X) and after (Y) escape.fig'])
%   saveas(h4,['single IBI before (X) and after (Y) escape.png']) 

%% nBout per min

h5=figure(5);
title ('mean nBout over section time');
       
            nFrames=60634;
            fps=100;
            TimeWindow= 60; %in sec
            Period= nFrames/(fps*TimeWindow); %Period = (210000/(350*60);
% G2

G2nBoutPerMin=[];
G2MeanPerMin=[];
for l=1:length(Fish_G2)
    %fprintf("-- Fish %d --\n",l);
    
    for z= 1:Period;
        %fprintf(" %d min %d\n",z);
        
        %find bout position in the time window
        G2nBoutPerMin{z}{l}= length(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] < (fps*TimeWindow).*z) ));
        
        G2MeanPerMin(z)= mean( cell2mat(G2nBoutPerMin{z}));
        G2_SEM(z)=std(cell2mat(G2nBoutPerMin{z}))/sqrt(length(Fish_G2));
    end;
end;

plot(1:Period, G2MeanPerMin,'bo-');hold on;
errorbar(1:Period, G2MeanPerMin, G2_SEM,'b'); hold on;
xlabel("min");
ylabel('Mean nBoutPerMin');hold on;
%grid();
title(['Mean nBoutPerMin']);hold on;

% % G1
% G1nBoutPerMin=[];
% G1MeanPerMin=[];
% for l=1:length(Fish_G1)
%     %fprintf("-- Fish %d --\n",l);
%     
%     for z= 1:Period;
%         %fprintf(" %d min %d\n",z);
%         
%         %find bout position in the time window
%         G1nBoutPerMin{z}{l}= length(find( ([datasetPerBout( allindex{Fish_G1(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G1(l)} ).BoutStart] < (fps*TimeWindow).*z) ));
%         
%         G1MeanPerMin(z)= mean( cell2mat(G1nBoutPerMin{z}));
%         G1_SEM(z)=std(cell2mat(G1nBoutPerMin{z}))/sqrt(length(Fish_G1));
%     end;
% end;
% 
% plot(1:Period, G1MeanPerMin,'go-');hold on;
% errorbar(1:Period, G1MeanPerMin, G1_SEM,'g'); hold on;
% xlabel("min");
% ylabel('Mean nBoutPerMin');hold on;
% %grid();
% title(['Mean nBoutPerMin']);hold on;

%G0
G0nBoutPerMin=[];
G0MeanPerMin=[];
for l=1:length(Fish_G0)
    %fprintf("-- Fish %d --\n",l);
    
    for z= 1:Period;
        %fprintf(" %d min %d\n",z);
        
        %find bout position in the time window
        G0nBoutPerMin{z}{l}= length(find( ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] < (fps*TimeWindow).*z) ));
        
        G0MeanPerMin(z)= mean( cell2mat(G0nBoutPerMin{z}));
        G0_SEM(z)=std(cell2mat(G0nBoutPerMin{z}))/sqrt(length(Fish_G0));
    end;
end;

plot(1:Period, G0MeanPerMin,'ro-');
errorbar(1:Period, G0MeanPerMin, G0_SEM,'r'); hold on;
xlabel("min");
ylabel('Mean nBoutPerMin');
%grid();
title(['Mean nBoutPerMin']);hold off;
legend('-/-','+/+');

saveas(h5,['Mean nBoutPerMin.fig'])
saveas(h5,['Mean nBoutPerMin.png'])  

%% IBI per 10s

h6=figure(6);
title ('mean nBout over section time');
       
            nFrames=60634;
            fps= unique([datasetPerBout(:).fps]);
            TimeWindow= 60; %in sec
            Period= nFrames/(fps*TimeWindow); %Period = (210000/(350*60);
% G2

G2IBI_time=[];
G2meanIBI_time=[];
G2meanFishIBI_time=[];
for z= 1:Period-1;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G2)
        %fprintf("-- Fish %d --\n",l);
        l
       %while [datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1) & [datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z;
        if z==1;
            idx_TimeWindow= index{Fish_G2(l)}(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2IBI_time{z}{l} = [datasetPerBout(idx_TimeWindow).InstantaneousIBI];
        else
            
         try
            idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2IBI_time{z}{l} = [datasetPerBout(idx_TimeWindow).InstantaneousIBI];
       
            
         catch
            numel(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty cell array all IBIs for time' num2str(z) 'Fish' num2str(l)]);
            G2IBI_time{z}{l}=10;
         end;
          
        end;
 
       %end;
        if isempty(G2IBI_time{z}{l});
            G2IBI_time{z}{l}=10;
        end;
       
         G2meanIBI_time{z}{l}= mean(G2IBI_time{1,z}{1,l},'omitnan');
         G2meanFishIBI_time(z)= mean(cell2mat(G2meanIBI_time{1,z}),'omitnan');
         G2_SEM(z)=std(cell2mat(G2meanIBI_time{1,z}))/sqrt(length(Fish_G2));
    end;
    
end;

%plot
plot(1:(Period-1), G2meanFishIBI_time,'bo-');hold on;
%errorbar(1:(Period-1), G2meanFishIBI_time, G2_SEM,'b'); hold on;
xlabel("sec");
ylabel('Mean IBI/10s');hold on;
%grid();
title(['Mean IBI/10s']);


%G1
G1IBI_time=[];
G1meanIBI_time=[];
G1meanFishIBI_time=[];
for z= 1:Period-1;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G1)
        %fprintf("-- Fish %d --\n",l);
        l
       
        if z==1;
            idx_TimeWindow= index{Fish_G1(l)}(find( ([datasetPerBout( allindex{Fish_G1(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G1(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G1IBI_time{z}{l} = [datasetPerBout(idx_TimeWindow).InstantaneousIBI];
        else
            
         try
            idx_TimeWindow= allindex{Fish_G1(l)}(find( ([datasetPerBout( allindex{Fish_G1(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G1(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G1IBI_time{z}{l} = [datasetPerBout(idx_TimeWindow).InstantaneousIBI];
            
         catch
            numel(find( ([datasetPerBout( allindex{Fish_G1(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G1(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty cell array all IBIs for time' num2str(z) 'Fish' num2str(l)]);
            G1IBI_time{z}{l}=10;
         end;
          
        end;
 
     
        if isempty(G1IBI_time{z}{l});
            G1IBI_time{z}{l}=10;
        end;
        
        G1meanIBI_time{z}{l}= mean(G1IBI_time{1,z}{1,l},'omitnan');
        G1meanFishIBI_time(z)= mean(cell2mat(G1meanIBI_time{1,z}),'omitnan');
        G1_SEM(z)=std(cell2mat(G1meanIBI_time{1,z}))/sqrt(length(Fish_G1));
    end;
end;

plot(1:(Period-1), G1meanFishIBI_time,'go-');hold on;
%errorbar(1:(Period-1), G1meanFishIBI_time, G1_SEM,'g'); hold on;
xlabel("sec");
ylabel('Mean IBI/10s');hold on;
%grid();
title(['Mean IBI/10s']);

%G0
G0IBI_time=[];
G0meanIBI_time=[];
G0meanFishIBI_time=[];
for z= 1:Period-1;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G0)
        %fprintf("-- Fish %d --\n",l);
        l
       
        if z==1;
            idx_TimeWindow= index{Fish_G0(l)}(find( ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G0IBI_time{z}{l} = [datasetPerBout(idx_TimeWindow).InstantaneousIBI];
        else
            
         try
            idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G0IBI_time{z}{l} = [datasetPerBout(idx_TimeWindow).InstantaneousIBI];
            
         catch
            numel(find( ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty cell array all IBIs for time' num2str(z) 'Fish' num2str(l)]);
            G0IBI_time{z}{l}=10;
         end;
          
        end;
 
      
        if isempty(G0IBI_time{z}{l});
            G0IBI_time{z}{l}=10;
        end;
        
       G0meanIBI_time{z}{l}= mean(G0IBI_time{1,z}{1,l},'omitnan');
       G0meanFishIBI_time(z)= mean(cell2mat(G0meanIBI_time{1,z}),'omitnan');
       G0_SEM(z)=std(cell2mat(G0meanIBI_time{1,z}))/sqrt(length(Fish_G0));
    end;
end;
 

plot(1:(Period-1), G0meanFishIBI_time,'ro-');hold on;
%errorbar(1:(Period-1), G0meanFishIBI_time, G0_SEM,'r'); hold on;
xlabel("sec");
ylabel('Mean IBI/10s');hold on;
%grid();
title(['Mean IBI/10s']);hold off;

saveas(h6,['IBIs per 10s.fig'])
saveas(h6,['IBIs per 10s.png'])

