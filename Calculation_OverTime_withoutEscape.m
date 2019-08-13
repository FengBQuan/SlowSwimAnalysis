
%% Clean workspace and load data
close all
clear variables
clc
load('matlab_workspace_SST_4manip_SlowSwim.mat')

% good swimers
    
% set new datasetSlowSwimBouts
%NoSwimFish= find([datasetPerFish(:).MeanIBI]> 10);
GoodSwimers=find(~([datasetPerFish(:).MeanIBI]> 10));
datasetGoodFish=datasetPerFish(GoodSwimers);

Fish = unique([datasetGoodFish(:).Condition]);
NumberFish=length(Fish);


% find EscapeWindow
EscapeWindow = [30000 30350];


EscapeBout=find((([datasetPerBout(:).BoutStart]> EscapeWindow(1)) & ([datasetPerBout(:).BoutStart]< EscapeWindow(2)) ))   ;
SlowSwimBouts=find(~(([datasetPerBout(:).BoutStart]> EscapeWindow(1)) & ([datasetPerBout(:).BoutStart]< EscapeWindow(2)) ) );
datasetSlowSwimBouts= datasetPerBout(SlowSwimBouts);


for i=1:NumberFish;

allindex{Fish(i)}= find(~([datasetSlowSwimBouts(:).Condition]-Fish(i)));

index{Fish(i)}= find(~([datasetSlowSwimBouts(:).Condition]-Fish(i)));

if numel(index{Fish(i)})>0;
    index{Fish(i)}(1)=[];
else if numel(index{Fish(i)})==0;
        index{Fish(i)}=[];
    end
end

IBI_pre{Fish(i)} = [datasetSlowSwimBouts(index{Fish(i)}(find([datasetSlowSwimBouts(index{Fish(i)}).BoutStart]< EscapeWindow(1)))).InstantaneousIBI];
   
%    for h=1:length(allindex{Fish(i)});
%         
%         display(['currently processing fish ' num2str(i)])
%         display(['currently processing bout number ' num2str(h)])
%         
%         TailAngle{Fish(i)}{h}=57.2958*[datasetSlowSwimBouts(allindex{Fish(i)}(h)).TailAngle_smoothed]';
% 
%         Bend_Amplitude{Fish(i)}{h} = 57.2958*[datasetSlowSwimBouts(allindex{Fish(i)}(h)).Bend_Amplitude];
%         
%         Max_Amp{Fish(i)}{h}=max(Bend_Amplitude{Fish(i)}{h});
%         
%         %idx_TurnBout= allindex{Fish(i)}(find( max(abs(57.2958*[datasetPerBout(allindex{Fish(i)}).Bend_Amplitude])) > 30 ));
%       
%     end;

end


Fish_temp=Fish;


FishGeno=([datasetGoodFish.Genotype]);

Fish_G2=Fish_temp(find(FishGeno( find( Fish_temp ) )==2));
Fish_G1=Fish_temp(find(FishGeno( find( Fish_temp ) )==1));
Fish_G0=Fish_temp(find(FishGeno( find( Fish_temp ) )==0));

save('Analysis_SST_4manip_without249frames_20180912_goodFish.mat')

            nFrames= 60385;
            fps= unique([datasetSlowSwimBouts(:).fps]);
            TimeWindow= 120; %in sec
            Period= nFrames/(fps*TimeWindow); %Period = (60034/(100*60);

%% IBI over time


% G2---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
 
G2IBI_time=[];
% G2medianIBI_time=[];
% G2medianFishIBI_time=[];
% G2_IBI_SEM=[];
 
for z= 1:Period-1;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G2)
        %fprintf("-- Fish %d --\n",l);
        l
        
       %while [datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1) & [datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z;
        if z==1;
            idx_TimeWindow= index{Fish_G2(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2_IBI_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).InstantaneousIBI];
        else
            
         try
            idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2_IBI_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).InstantaneousIBI];
  
         catch
            numel(find( ([datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty cell array all IBIs for time' num2str(z) 'Fish' num2str(l)]);
            G2_IBI_time{z}{l}=10;
         end;
          
        end;
 
       %end;
        if isempty(G2_IBI_time{z}{l});
            G2_IBI_time{z}{l}=10;
        end;
       
         G2medianIBI_time{z}{l}= median(G2_IBI_time{1,z}{1,l},'omitnan');
         G2medianFishIBI_time(z)= median(cell2mat(G2medianIBI_time{1,z}),'omitnan');
         
         G2medianFishBoutFreq_time(z)= 1/G2medianFishIBI_time(z);
         G2_IBI_SEM(z)=std(cell2mat(G2medianIBI_time{1,z}))/sqrt(length(Fish_G2));
    end;
    
end;
 
% G0---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
 
G0IBI_time=[];
% G0medianIBI_time=[];
% G0medianFishIBI_time=[];
% G0_IBI_SEM=[];
 
for z= 1:Period-1;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G0)
        %fprintf("-- Fish %d --\n",l);
        l

        
        if z==1;
            idx_TimeWindow= index{Fish_G0(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G0IBI_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).InstantaneousIBI];
        else
            
         try
            idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G0IBI_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).InstantaneousIBI];
            
         catch
            numel(find( ([datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty cell array all IBIs for time' num2str(z) 'Fish' num2str(l)]);
            G0IBI_time{z}{l}=10;
         end;
          
        end;
 
      
        if isempty(G0IBI_time{z}{l});
            G0IBI_time{z}{l}=10;
        end;
        
       G0medianIBI_time{z}{l}= median(G0IBI_time{1,z}{1,l},'omitnan');
       G0medianFishIBI_time(z)= median(cell2mat(G0medianIBI_time{1,z}),'omitnan');
       
       G0medianFishBoutFreq_time(z)= 1/G0medianFishIBI_time(z);
       G0_IBI_SEM(z)=std(cell2mat(G0medianIBI_time{1,z}))/sqrt(length(Fish_G0));
       
    end;
end;
 
%save('BoutFreqErrorbarCheck.mat', 'G0medianIBI_time','G2medianIBI_time')
save('workspace_BoutFrequency.mat')
 
%% BoutDuration over time calculation
 
G2_BoutDuration_time=[];
% G2medianBoutDuration_time=[];
% G2medianFishBoutDuration_time=[];
%  
% G2_Total_BoutDuration_time=[];
% G2medianFishTotalBD_time=[];

for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G2)
        %fprintf("-- Fish %d --\n",l);
        l
       %while [datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1) & [datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z;
%         if z==1;
%             idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
%             G2_BoutDuration_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).BoutDuration];
        %else
            
         try
            idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2_BoutDuration_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).BoutDuration];
       
            
         catch
            numel(find( ([datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty cell array all IBIs for time' num2str(z) 'Fish' num2str(l)]);
            G2_BoutDuration_time{z}{l}=0;
         end;
          
        %end;
 
       %end;
        if isempty(G2_BoutDuration_time{z}{l});
            G2_BoutDuration_time{z}{l}= 0;
        end;
         
         G2medianBoutDuration_time{z}{l}= median(G2_BoutDuration_time{1,z}{1,l},'omitnan'); % median BD for fish l
         G2medianFishBoutDuration_time(z)= median(cell2mat(G2medianBoutDuration_time{1,z}),'omitnan');%median of all BD of all fish in period Z
         G2_SEM(z)=std(cell2mat(G2medianBoutDuration_time{1,z}))/sqrt(length(Fish_G2));
         
         G2_Total_BoutDuration_time{z}{l}= sum(G2_BoutDuration_time{1,z}{1,l},'omitnan'); % sum BD for fish l in period z
         G2medianFishTotalBD_time(z)= median(cell2mat(G2_Total_BoutDuration_time{1,z}),'omitnan');% median of total BD of all fish
         G2_Total_SEM(z)=std(cell2mat(G2_Total_BoutDuration_time{1,z}))/sqrt(length(Fish_G2));
    end;
    
end;
  
% G0---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
G0_BoutDuration_time=[];
% G0medianBoutDuration_time=[];
% G0medianFishBoutDuration_time=[];
%  
% G0_Total_BoutDuration_time=[];
% G0medianFishTotalBD_time=[];
 
for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G0)
        %fprintf("-- Fish %d --\n",l);
        l
       %while [datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1) & [datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z;
%         if z==1;
%             idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
%             G0_BoutDuration_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).BoutDuration];
%         else
            
         try
            idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G0_BoutDuration_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).BoutDuration];
       
            
         catch
            numel(find( ([datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty cell array all IBIs for time' num2str(z) 'Fish' num2str(l)]);
            G0_BoutDuration_time{z}{l}=0;
         end;
          
        %end;
 
       %end;
        if isempty(G0_BoutDuration_time{z}{l});
            G0_BoutDuration_time{z}{l}= 0;
        end;
         
         G0medianBoutDuration_time{z}{l}= median(G0_BoutDuration_time{1,z}{1,l},'omitnan'); % median BD for fish l
         G0medianFishBoutDuration_time(z)= median(cell2mat(G0medianBoutDuration_time{1,z}),'omitnan');%median of all BD of all fish in period Z
         G0_SEM(z)=std(cell2mat(G0medianBoutDuration_time{1,z}))/sqrt(length(Fish_G0));
         
         G0_Total_BoutDuration_time{z}{l}= sum(G0_BoutDuration_time{1,z}{1,l},'omitnan'); % sum BD for fish l in period z
         G0medianFishTotalBD_time(z)= median(cell2mat(G0_Total_BoutDuration_time{1,z}),'omitnan');% median of total BD of all fish
         G0_Total_SEM(z)=std(cell2mat(G0_Total_BoutDuration_time{1,z}))/sqrt(length(Fish_G0));
    end;
    
end;
 
save('workspace_BoutDuration.mat')

%% Distance over time calculations

G2_Distance_time=[];
G2medianDistance_time=[];
G2medianFishDistance_time=[];
 
G2_Total_Distance_time=[];
G2medianFishTotalDistance_time=[];
for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G2)
        %fprintf("-- Fish %d --\n",l);
        l
       %while [datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1) & [datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z;
        if z==1;
            idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2_Distance_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).TotalDistance];
        else
            
         try
            idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2_Distance_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).TotalDistance];
       
            
         catch
            numel(find( ([datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty cell array all IBIs for time' num2str(z) 'Fish' num2str(l)]);
            G2_Distance_time{z}{l}=0;
         end;
          
        end;
 
       %end;
        if isempty(G2_Distance_time{z}{l});
            G2_Distance_time{z}{l}= 0;
        end;
         
         G2medianDistance_time{z}{l}= median(G2_Distance_time{1,z}{1,l},'omitnan'); % median BD for fish l
         G2medianFishDistance_time(z)= median(cell2mat(G2medianDistance_time{1,z}),'omitnan');%median of all BD of all fish in period Z
         G2_SEM(z)=std(cell2mat(G2medianDistance_time{1,z}))/sqrt(length(Fish_G2));
         
         G2_Total_Distance_time{z}{l}= sum(G2_Distance_time{1,z}{1,l},'omitnan'); % sum BD for fish l in period z
         G2medianFishTotalDistance_time(z)= median(cell2mat(G2_Total_Distance_time{1,z}),'omitnan');% median of total BD of all fish
         G2_Total_SEM(z)=std(cell2mat(G2_Total_Distance_time{1,z}))/sqrt(length(Fish_G2));
    end;
    
end;
 
% G0---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
G0_Distance_time=[];
G0medianDistance_time=[];
G0medianFishDistance_time=[];
 
G0_Total_Distance_time=[];
G0medianFishTotalDistance_time=[];
for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G0)
        %fprintf("-- Fish %d --\n",l);
        l
       %while [datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1) & [datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z;
        if z==1;
            idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G0_Distance_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).TotalDistance];
        else
            
         try
            idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G0_Distance_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).TotalDistance];
       
            
         catch
            numel(find( ([datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty cell array all IBIs for time' num2str(z) 'Fish' num2str(l)]);
            G0_Distance_time{z}{l}=0;
         end;
          
        end;
 
       %end;
        if isempty(G0_Distance_time{z}{l});
            G0_Distance_time{z}{l}= 0;
        end;
         
         G0medianDistance_time{z}{l}= median(G0_Distance_time{1,z}{1,l},'omitnan'); % median BD for fish l
         G0medianFishDistance_time(z)= median(cell2mat(G0medianDistance_time{1,z}),'omitnan');%median of all BD of all fish in period Z
         G0_SEM(z)=std(cell2mat(G0medianDistance_time{1,z}))/sqrt(length(Fish_G0));
         
         G0_Total_Distance_time{z}{l}= sum(G0_Distance_time{1,z}{1,l},'omitnan'); % sum BD for fish l in period z
         G0medianFishTotalDistance_time(z)= median(cell2mat(G0_Total_Distance_time{1,z}),'omitnan');% median of total BD of all fish
         G0_Total_SEM(z)=std(cell2mat(G0_Total_Distance_time{1,z}))/sqrt(length(Fish_G0));
    end;
    
end;
 
save('workspace_BoutDistance.mat')

%% nBout per min calculation       

% G2---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
G2nBoutPerMin=[];
G2MedianPerMin=[];
for l=1:length(Fish_G2)
    %fprintf("-- Fish %d --\n",l);
    
    for z= 1:Period;
        %fprintf(" %d min %d\n",z);
        
        %find bout position in the time window
        G2nBoutPerMin{z}{l}= length(find( ([datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] < (fps*TimeWindow).*z) ));
        
        G2MedianPerMin(z)= median( cell2mat(G2nBoutPerMin{z}));
        G2_SEM(z)=std(cell2mat(G2nBoutPerMin{z}))/sqrt(length(Fish_G2));
    end;
end;
 
% G0---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
 
G0nBoutPerMin=[];
G0MedianPerMin=[];
for l=1:length(Fish_G0)
    %fprintf("-- Fish %d --\n",l);
    
    for z= 1:Period;
        %fprintf(" %d min %d\n",z);
        
        %find bout position in the time window
        G0nBoutPerMin{z}{l}= length(find( ([datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] < (fps*TimeWindow).*z) ));
        
        G0MedianPerMin(z)= median( cell2mat(G0nBoutPerMin{z}));
        G0_SEM(z)=std(cell2mat(G0nBoutPerMin{z}))/sqrt(length(Fish_G0));
    end;
end;

 save('workspace_nBout.mat')
%% Amplitude per min       

% G2---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
 
G2_BendAmplitude=[];
G2_MaxAmp_BendAmplitude=[];
G2_medianBout_MaxAmplitude=[];
G2_medianFish_BendAmplitude=[];
 
for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G2)
        %fprintf("-- Fish %d --\n",l);
        l
 
         if numel(find( ([datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty cell array for time' num2str(z) 'Fish' num2str(l)]);
            G2_BendAmplitude{z}{l}{h}=0;
            G2_MaxAmp_BendAmplitude{z}{l}{h}=0;
            G2_medianBout_MaxAmplitude{z}{l}=0;
            G2_medianFish_BendAmplitude{z}=0;
          
         else
         idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            
            for h=1:length(idx_TimeWindow)
                h
            display(['currently processing period ' num2str(z)])    
            display(['currently processing fish ' num2str(l)])
            display(['currently processing bout number ' num2str(h)])
            
            G2_BendAmplitude{z}{l}{h}= 57.2958*[datasetSlowSwimBouts(idx_TimeWindow(h)).Bend_Amplitude]; % all Amplitude in bout h, for fishl, in period z
%             for a= 1:length(G2_BendAmplitude{z}{l}{h})
%             if abs(G2_BendAmplitude{z}{l}{h}(a,1))<2;
%                 G2_BendAmplitude{z}{l}{h}(a,1)=2;
%             end
%             end
            %G2_medianAmp_BendAmplitude{z}{l}{h}=median(abs(G2_BendAmplitude{z}{l}{h}));%median Amplitude of bout h, for fishl, in period z
            G2_MaxAmp_BendAmplitude{z}{l}{h}=max(abs(G2_BendAmplitude{z}{l}{h}));
            end
            
         end
          %G2_medianBout_BendAmplitude{z}{l}= median(cell2mat(G2_medianAmp_BendAmplitude{1,z}{1,l})); % median of all bout for fish l
          G2_medianBout_MaxAmplitude{z}{l}= median(cell2mat(G2_MaxAmp_BendAmplitude{1,z}{1,l}));
          G2_medianFish_BendAmplitude{z}= median(cell2mat(G2_medianBout_MaxAmplitude{1,z}));%median of all fish in period Z
          G2_SEM(z)=std(cell2mat(G2_medianBout_MaxAmplitude{1,z}))/sqrt(length(Fish_G2));
         
        
    end
    
end;
 
% G0---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
G0_BendAmplitude=[];
G0_MaxAmp_BendAmplitude=[];
G0_medianBout_MaxAmplitude=[];
G0_medianFish_BendAmplitude=[];
 
 
for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G0)
        %fprintf("-- Fish %d --\n",l);
        l
 
         if numel(find( ([datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty Bout for time' num2str(z) 'Fish' num2str(l)]);
            G0_BendAmplitude{z}{l}{h}=0;
            G0_MaxAmp_BendAmplitude{z}{l}{h}=0;
            G0_medianBout_MaxAmplitude{z}{l}=0;
            G0_medianFish_BendAmplitude{z}=0;
          
         else
         idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            
            for h=1:length(idx_TimeWindow)
                h
            display(['currently processing period ' num2str(z)])    
            display(['currently processing fish ' num2str(l)])
            display(['currently processing bout number ' num2str(h)])
            
            G0_BendAmplitude{z}{l}{h}= 57.2958*[datasetSlowSwimBouts(idx_TimeWindow(h)).Bend_Amplitude];% all Amplitude in bout h, for fishl, in period z
%             for a= 1:length(G0_BendAmplitude{z}{l}{h})
%             if abs(G0_BendAmplitude{z}{l}{h}(a,1))<2;
%                 G0_BendAmplitude{z}{l}{h}(a,1)= 2;
%             end
%             end
            %G0_medianAmp_BendAmplitude{z}{l}{h}=median(abs(G0_BendAmplitude{z}{l}{h}));
             G0_MaxAmp_BendAmplitude{z}{l}{h}=max(abs(G0_BendAmplitude{z}{l}{h}));%median Amplitude of bout h, for fishl, in period z
            end
 
         end
          G0_medianBout_MaxAmplitude{z}{l}= median(cell2mat(G0_MaxAmp_BendAmplitude{1,z}{1,l})); % median of all bout for fish l
          G0_medianFish_BendAmplitude{z}= median(cell2mat(G0_medianBout_MaxAmplitude{1,z}));%median of all fish in period Z
          G0_SEM(z)=std(cell2mat(G0_medianBout_MaxAmplitude{1,z}))/sqrt(length(Fish_G0));
     
    end
    
end;

save('workspace_MaxiBendAmp.mat')
 
%%  Median NumberOfOscillations over time calculation

% G2---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
 
G2_NumberOfOscillations_time=[];
G2medianNumberOfOscillations_time=[];
G2medianFishNumberOfOscillations_time=[];
 
for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G2)
        %fprintf("-- Fish %d --\n",l);
        l
       %while [datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1) & [datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z;
        if z==1;
            idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2_NumberOfOscillations_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).NumberOfOscillations];
        else
            
         try
            idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2_NumberOfOscillations_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).NumberOfOscillations];
       
            
         catch
            numel(find( ([datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty cell array for time' num2str(z) 'Fish' num2str(l)]);
            G2_NumberOfOscillations_time{z}{l}=0;
         end;
          
        end;
 
       %end;
        if isempty(G2_NumberOfOscillations_time{z}{l});
            G2_NumberOfOscillations_time{z}{l}= 0;
        end;
         
         G2medianNumberOfOscillations_time{z}{l}= median(G2_NumberOfOscillations_time{1,z}{1,l},'omitnan'); % median BD for fish l
         G2medianFishNumberOfOscillations_time(z)= median(cell2mat(G2medianNumberOfOscillations_time{1,z}),'omitnan');%median of all BD of all fish in period Z
         G2_SEM(z)=std(cell2mat(G2medianNumberOfOscillations_time{1,z}))/sqrt(length(Fish_G2));
         
    end;
    
end;
  
% G0---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
G0_NumberOfOscillations_time=[];
G0medianNumberOfOscillations_time=[];
G0medianFishNumberOfOscillations_time=[];
 
for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G0)
        %fprintf("-- Fish %d --\n",l);
        l
       %while [datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1) & [datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z;
        if z==1;
            idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G0_NumberOfOscillations_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).NumberOfOscillations];
        else
            
         try
            idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G0_NumberOfOscillations_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).NumberOfOscillations];
       
            
         catch
            numel(find( ([datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty cell array all IBIs for time' num2str(z) 'Fish' num2str(l)]);
            G0_NumberOfOscillations_time{z}{l}=0;
         end;
          
        end;
 
       %end;
        if isempty(G0_NumberOfOscillations_time{z}{l});
            G0_NumberOfOscillations_time{z}{l}= 0;
        end;
         
         G0medianNumberOfOscillations_time{z}{l}= median(G0_NumberOfOscillations_time{1,z}{1,l},'omitnan'); % median BD for fish l
         G0medianFishNumberOfOscillations_time(z)= median(cell2mat(G0medianNumberOfOscillations_time{1,z}),'omitnan');%median of all BD of all fish in period Z
         G0_SEM(z)=std(cell2mat(G0medianNumberOfOscillations_time{1,z}))/sqrt(length(Fish_G0));
         
    end;
    
end;
 

save('workspace_medianNumOfOsc.mat')

 
%% Mean NumberOfOscillations over time calculation

% G2---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
 
G2_NumberOfOscillations_time=[];
G2meanNumberOfOscillations_time=[];
G2meanFishNumberOfOscillations_time=[];
 
for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G2)
        %fprintf("-- Fish %d --\n",l);
        l
       %while [datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1) & [datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z;
        if z==1;
            idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2_NumberOfOscillations_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).NumberOfOscillations];
        else
            
         try
            idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2_NumberOfOscillations_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).NumberOfOscillations];
       
            
         catch
            numel(find( ([datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty cell array all IBIs for time' num2str(z) 'Fish' num2str(l)]);
            G2_NumberOfOscillations_time{z}{l}=0;
         end;
          
        end;
 
       %end;
        if isempty(G2_NumberOfOscillations_time{z}{l});
            G2_NumberOfOscillations_time{z}{l}= 0;
        end;
         
         G2meanNumberOfOscillations_time{z}{l}= mean(G2_NumberOfOscillations_time{1,z}{1,l},'omitnan'); % mean BD for fish l
         G2meanFishNumberOfOscillations_time(z)= mean(cell2mat(G2meanNumberOfOscillations_time{1,z}),'omitnan');%mean of all BD of all fish in period Z
         G2_SEM(z)=std(cell2mat(G2meanNumberOfOscillations_time{1,z}))/sqrt(length(Fish_G2));
 
    end;
    
end;
  
% G0---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
G0_NumberOfOscillations_time=[];
G0meanNumberOfOscillations_time=[];
G0meanFishNumberOfOscillations_time=[];
 
for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G0)
        %fprintf("-- Fish %d --\n",l);
        l
       %while [datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1) & [datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z;
        if z==1;
            idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G0_NumberOfOscillations_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).NumberOfOscillations];
        else
            
         try
            idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G0_NumberOfOscillations_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).NumberOfOscillations];
       
            
         catch
            numel(find( ([datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty cell array all IBIs for time' num2str(z) 'Fish' num2str(l)]);
            G0_NumberOfOscillations_time{z}{l}=0;
         end;
          
        end;
 
       %end;
        if isempty(G0_NumberOfOscillations_time{z}{l});
            G0_NumberOfOscillations_time{z}{l}= 0;
        end;
         
         G0meanNumberOfOscillations_time{z}{l}= mean(G0_NumberOfOscillations_time{1,z}{1,l},'omitnan'); % mean BD for fish l
         G0meanFishNumberOfOscillations_time(z)= mean(cell2mat(G0meanNumberOfOscillations_time{1,z}),'omitnan');%mean of all BD of all fish in period Z
         G0_SEM(z)=std(cell2mat(G0meanNumberOfOscillations_time{1,z}))/sqrt(length(Fish_G0));
    end;
    
end;

save('workspace_medianNumOfOsc.mat')

 
 
%% TBF(NumberOfOscillations/Duration) over time calculation
      
% G2---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
 
G2_TBF_time=[];
G2medianTBF_time=[];
G2medianFishTBF_time=[];
 
for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G2)
        %fprintf("-- Fish %d --\n",l);
        l
       %while [datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1) & [datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z;
        if z==1;
            idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2_TBF_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).NumberOfOscillations]/[datasetSlowSwimBouts(idx_TimeWindow).BoutDuration];
        else
          
         try
            idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2_TBF_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).NumberOfOscillations]/[datasetSlowSwimBouts(idx_TimeWindow).BoutDuration];
       
            
         catch
            numel(find( ([datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty cell array all IBIs for time' num2str(z) 'Fish' num2str(l)]);
            G2_TBF_time{z}{l}=0;
         end;
          
        end;
 
       %end;
        if isempty(G2_TBF_time{z}{l});
            G2_TBF_time{z}{l}= 0;
        end;
         
         G2medianTBF_time{z}{l}= median(G2_TBF_time{1,z}{1,l},'omitnan'); % median BD for fish l
         G2medianFishTBF_time(z)= median(cell2mat(G2medianTBF_time{1,z}),'omitnan');%median of all BD of all fish in period Z
         G2_SEM(z)=std(cell2mat(G2medianTBF_time{1,z}))/sqrt(length(Fish_G2));
 
    end;
    
end;
  
% G0---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
G0_TBF_time=[];
G0medianTBF_time=[];
G0medianFishTBF_time=[];
 
for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G0)
        %fprintf("-- Fish %d --\n",l);
        l
       %while [datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1) & [datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z;
        if z==1;
            idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G0_TBF_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).NumberOfOscillations]/[datasetSlowSwimBouts(idx_TimeWindow).BoutDuration];
        else
            
         try
            idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G0_TBF_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).NumberOfOscillations]/[datasetSlowSwimBouts(idx_TimeWindow).BoutDuration];
       
            
         catch
            numel(find( ([datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            disp(['empty cell array all IBIs for time' num2str(z) 'Fish' num2str(l)]);
            G0_TBF_time{z}{l}=0;
         end;
          
        end;
 
       %end;
        if isempty(G0_TBF_time{z}{l});
            G0_TBF_time{z}{l}= 0;
        end;
         
         G0medianTBF_time{z}{l}= median(G0_TBF_time{1,z}{1,l},'omitnan'); % median BD for fish l
         G0medianFishTBF_time(z)= median(cell2mat(G0medianTBF_time{1,z}),'omitnan');%median of all BD of all fish in period Z
         G0_SEM(z)=std(cell2mat(G0medianTBF_time{1,z}))/sqrt(length(Fish_G0));
         
    end;
    
end;


save('workspace_TBF.mat')


%% Speed 

% G2---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
 
G2_Speed_time=[];
G2medianSpeed_time=[];
G2medianFishSpeed_time=[];
 
G2_Total_Speed_time=[];
G2medianFishTotalSpeed_time=[];
for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G2)
        %fprintf("-- Fish %d --\n",l);
        l
       %while [datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1) & [datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z;
        if z==1;
            idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2_Speed_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).Speed];
        else
            
         try
            idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2_Speed_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).Speed];
       
            
         catch
            numel(find( ([datasetSlowSwimBouts( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty cell array all IBIs for time' num2str(z) 'Fish' num2str(l)]);
            G2_Speed_time{z}{l}=0;
         end;
          
        end;
 
       %end;
        if isempty(G2_Speed_time{z}{l});
            G2_Speed_time{z}{l}= 0;
        end;
         
         G2medianSpeed_time{z}{l}= median(G2_Speed_time{1,z}{1,l},'omitnan'); % median BD for fish l
         G2medianFishSpeed_time(z)= median(cell2mat(G2medianSpeed_time{1,z}),'omitnan');%median of all BD of all fish in period Z
         G2_SEM(z)=std(cell2mat(G2medianSpeed_time{1,z}))/sqrt(length(Fish_G2));
         
         G2_Total_Speed_time{z}{l}= sum(G2_Speed_time{1,z}{1,l},'omitnan'); % sum BD for fish l in period z
         G2medianFishTotalSpeed_time(z)= median(cell2mat(G2_Total_Speed_time{1,z}),'omitnan');% median of total BD of all fish
         G2_Total_SEM(z)=std(cell2mat(G2_Total_Speed_time{1,z}))/sqrt(length(Fish_G2));
    end;
    
end;
 
% G0---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
G0_Speed_time=[];
G0medianSpeed_time=[];
G0medianFishSpeed_time=[];
 
G0_Total_Speed_time=[];
G0medianFishTotalSpeed_time=[];
for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G0)
        %fprintf("-- Fish %d --\n",l);
        l
       %while [datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1) & [datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z;
        if z==1;
            idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G0_Speed_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).Speed];
        else
            
         try
            idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G0_Speed_time{z}{l} = [datasetSlowSwimBouts(idx_TimeWindow).Speed];
       
            
         catch
            numel(find( ([datasetSlowSwimBouts( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty cell array all IBIs for time' num2str(z) 'Fish' num2str(l)]);
            G0_Speed_time{z}{l}=0;
         end;
          
        end;
 
       %end;
        if isempty(G0_Speed_time{z}{l});
            G0_Speed_time{z}{l}= 0;
        end;
         
         G0medianSpeed_time{z}{l}= median(G0_Speed_time{1,z}{1,l},'omitnan'); % median BD for fish l
         G0medianFishSpeed_time(z)= median(cell2mat(G0medianSpeed_time{1,z}),'omitnan');%median of all BD of all fish in period Z
         G0_SEM(z)=std(cell2mat(G0medianSpeed_time{1,z}))/sqrt(length(Fish_G0));
         
         G0_Total_Speed_time{z}{l}= sum(G0_Speed_time{1,z}{1,l},'omitnan'); % sum BD for fish l in period z
         G0medianFishTotalSpeed_time(z)= median(cell2mat(G0_Total_Speed_time{1,z}),'omitnan');% median of total BD of all fish
         G0_Total_SEM(z)=std(cell2mat(G0_Total_Speed_time{1,z}))/sqrt(length(Fish_G0));
    end;
    
end;
 
save('workspace_Speed.mat')

%% Figure

h1=figure(1)
title ('Parameters over section time');hold on;

subplot(2,4,1)
title ('Bout Frequency over section time');hold on;
 
plot(1:(Period-1), G2medianFishBoutFreq_time,'bo-');hold on;
%errorbar(1:(Period-1), G2medianFishBoutFreq_time, G2_IBI_SEM,'b'); hold on;       
plot(1:(Period-1), G0medianFishBoutFreq_time,'ro-');hold on;
%errorbar(1:(Period-1), G0medianFishBoutFreq_time, G0_IBI_SEM,'r'); hold on;
xlabel("Period per 2 min");
ylabel('Hz');hold on;
%grid();
hold off;


subplot(2,4,2)
title(['BoutDuration over time']);hold on;
 
plot(1:(Period), G2medianFishBoutDuration_time,'bo-');hold on;
errorbar(1:(Period), G2medianFishBoutDuration_time, G2_SEM,'b'); hold on;
 
plot(1:(Period), G0medianFishBoutDuration_time,'ro-');hold on;
errorbar(1:(Period), G0medianFishBoutDuration_time, G0_SEM,'r'); hold on;
xlabel("Period per 2 min");hold on;
ylabel('BoutDuration in sec');hold on;
%grid();
hold off;


subplot(2,4,3)
title(['BoutDistance over time']);hold on;

plot(1:(Period), G2medianFishDistance_time,'bo-');hold on;
errorbar(1:(Period), G2medianFishDistance_time, G2_SEM,'b'); hold on;

plot(1:(Period), G0medianFishDistance_time,'ro-');hold on;
errorbar(1:(Period), G0medianFishDistance_time, G0_SEM,'r'); hold on;
xlabel("Period per 2 min");
ylabel('BoutDistance in mm');hold on;
%grid();
hold off;


subplot(2,4,4)
title ('median nBout over section time');hold on;
 
plot(1:Period, G2MedianPerMin,'bo-');hold on;
errorbar(1:Period, G2MedianPerMin, G2_SEM,'b'); hold on;
 
plot(1:Period, G0MedianPerMin,'ro-');
errorbar(1:Period, G0MedianPerMin, G0_SEM,'r'); hold on;
xlabel("Period per 2 min");
ylabel('nBout');
%legend('-/-','+/+');
%grid();
hold off;



subplot(2,4,5)
title(['Median Bend_Amplitude per min']);hold on;
 
plot(1:(Period), cell2mat(G2_medianFish_BendAmplitude),'bo-');hold on;
errorbar(1:(Period), cell2mat(G2_medianFish_BendAmplitude), G2_SEM,'b'); hold on;
 
plot(1:(Period), cell2mat(G0_medianFish_BendAmplitude),'ro-');hold on;
errorbar(1:(Period), cell2mat(G0_medianFish_BendAmplitude), G0_SEM,'r'); hold on;
xlabel("Period per 2 min");
ylabel('Degree');hold on;
 
hold off;


subplot(2,4,6)
title(['BoutSpeed over time']);hold on;
%G2
plot(1:(Period), G2medianFishSpeed_time,'bo-');hold on;
errorbar(1:(Period), G2medianFishSpeed_time, G2_SEM,'b'); hold on;
%G0
plot(1:(Period), G0medianFishSpeed_time,'ro-');hold on;
errorbar(1:(Period), G0medianFishSpeed_time, G0_SEM,'r'); hold on;
 
xlabel("Period per 2 min");
ylabel('BoutSpeed in mm');hold on;
%grid();
hold off;


subplot(2,4,7)
title([' Mean NumberOfOscillations over time']);hold on;
 
plot(1:(Period), G2meanFishNumberOfOscillations_time,'bo-');hold on;
errorbar(1:(Period), G2meanFishNumberOfOscillations_time, G2_SEM,'b'); hold on;
 
plot(1:(Period), G0meanFishNumberOfOscillations_time,'ro-');hold on;
errorbar(1:(Period), G0meanFishNumberOfOscillations_time, G0_SEM,'r'); hold on;
 
xlabel("Period per 2 min");hold on;
ylabel('NumberOfOscillations in sec');hold on;
%grid();
hold off;

subplot(2,4,8)
title(['TBF(NumberOfOscillations/Duration) over time']);hold on;
 
plot(1:(Period), G2medianFishTBF_time,'bo-');hold on;
errorbar(1:(Period), G2medianFishTBF_time, G2_SEM,'b'); hold on;
 
plot(1:(Period), G0medianFishTBF_time,'ro-');hold on;
errorbar(1:(Period), G0medianFishTBF_time, G0_SEM,'r'); hold on;
 
xlabel("Period per 2 min");hold on;
ylabel('TBF in Hz');hold on;
%grid();
hold off;




saveas(h1,['Parameters per 2min.fig'])
saveas(h1,['Parameters per 2min.epsc'])



