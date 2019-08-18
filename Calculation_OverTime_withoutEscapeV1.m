
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

%% create outputTables


output_Distance = struct( 'Variable', [], 'Clutch', [], 'Fish', [], 'FishGeno',[], ...
    'Period1', [], 'Period2',[],'Period3',[],'Period4',[],'Period5',[]);
output_Duration = struct( 'Variable', [], 'Clutch', [], 'Fish', [], 'FishGeno',[], ...
    'Period1', [], 'Period2',[],'Period3',[],'Period4',[],'Period5',[]);
output_Speed = struct( 'Variable', [], 'Clutch', [], 'Fish', [], 'FishGeno',[], ...
    'Period1', [], 'Period2',[],'Period3',[],'Period4',[],'Period5',[]);
output_NumOfOsc = struct( 'Variable', [], 'Clutch', [], 'Fish', [], 'FishGeno',[], ...
    'Period1', [], 'Period2',[],'Period3',[],'Period4',[],'Period5',[]);
output_TBF = struct( 'Variable', [], 'Clutch', [], 'Fish', [], 'FishGeno',[], ...
    'Period1', [], 'Period2',[],'Period3',[],'Period4',[],'Period5',[]);
output_MaxBendAmp = struct( 'Variable', [], 'Clutch', [], 'Fish', [], 'FishGeno',[], ...
    'Period1', [], 'Period2',[],'Period3',[],'Period4',[],'Period5',[]);




for i=1:NumberFish

allindex{Fish(i)}= find(~([datasetSlowSwimBouts(:).Condition]-Fish(i)));

index{Fish(i)}= find(~([datasetSlowSwimBouts(:).Condition]-Fish(i)));

if numel(index{Fish(i)})>0;
    index{Fish(i)}(1)=[];
else if numel(index{Fish(i)})==0;
        index{Fish(i)}=[];
    end
end

end

Fish_temp=Fish;
FishGeno=([datasetGoodFish.Genotype]);

Fish_G2=Fish_temp(find(FishGeno( find( Fish_temp ) )==2));
Fish_G1=Fish_temp(find(FishGeno( find( Fish_temp ) )==1));
Fish_G0=Fish_temp(find(FishGeno( find( Fish_temp ) )==0));

%% set Period

            nFrames= 60000;
            fps= unique([datasetSlowSwimBouts(:).fps]);
            TimeWindow= 120; %in sec
            Period= nFrames/(fps*TimeWindow); %Period = (60034/(100*60); 

%% Other parameters

Distance_time =[];
BoutDuration_time=[];
Speed_time=[];
NumberOfOscillations_time=[];
TBF_time=[];


medianDistance_time=[];
medianBoutDuration_time=[];
medianSpeed_time=[];
meanNumberOfOscillations_time=[];
medianTBF_time=[];


    for i=1:NumberFish
   
        for z= 1:Period
    
    display([' BoutFrequency currently processing fish ' num2str(i)])
    display([' BoutFrequency currently processing period ' num2str(z)])
    
        if  numel(index{Fish(i)})==0;
            Distance_time{Fish(i)}{z} =nan;
            BoutDuration_time{Fish(i)}{z}=nan;
            Speed_time{Fish(i)}{z}=nan;
            NumberOfOscillations_time{Fish(i)}{z}=nan;
            TBF_time{Fish(i)}{z}=nan;
        else
        
         try
            idx_TimeWindow= allindex{Fish(i)}(find( ([datasetSlowSwimBouts( allindex{Fish(i)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish(i)}).BoutStart] < (fps*TimeWindow).*z) ) );
            
             Distance_time{Fish(i)}{z} = [datasetSlowSwimBouts(idx_TimeWindow).TotalDistance]; 
            BoutDuration_time{Fish(i)}{z} = [datasetSlowSwimBouts(idx_TimeWindow).BoutDuration];
            Speed_time{Fish(i)}{z} = [datasetSlowSwimBouts(idx_TimeWindow).Speed];
            NumberOfOscillations_time{Fish(i)}{z} = [datasetSlowSwimBouts(idx_TimeWindow).NumberOfOscillations];
            TBF_time{Fish(i)}{z}= [datasetSlowSwimBouts(idx_TimeWindow).NumberOfOscillations]/[datasetSlowSwimBouts(idx_TimeWindow).BoutDuration];
  
         catch
            numel(find( ([datasetSlowSwimBouts( allindex{Fish(i)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish(i)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty cell array all parameters for time' num2str(z) 'Fish' num2str(i)]);
            Distance_time{Fish(i)}{z}=0;
            BoutDuration_time{Fish(i)}{z}=0;
            Speed_time{Fish(i)}{z}=0;
            NumberOfOscillations_time{Fish(i)}{z}=0;
            TBF_time{Fish(i)}{z}=0;
            
         end
        end
          

       
       % calculate median of all bouts for each fish
         medianDistance_time{Fish(i)}{z}= median(Distance_time{Fish(i)}{z},'omitnan'); 
         medianBoutDuration_time{Fish(i)}{z}= median(BoutDuration_time{Fish(i)}{z},'omitnan'); 
         medianSpeed_time{Fish(i)}{z}= median(Speed_time{Fish(i)}{z},'omitnan');
         meanNumberOfOscillations_time{Fish(i)}{z}= mean(NumberOfOscillations_time{Fish(i)}{z},'omitnan');
         medianTBF_time{Fish(i)}{z}= median(TBF_time{Fish(i)}{z},'omitnan');
        end
    end
         
    

%% Amplitude per min       


BendAmplitude=[];
MaxAmp_BendAmplitude=[];
medianBout_MaxAmplitude=[];

for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for i=1:NumberFish
        %fprintf("-- Fish %d --\n",l);
        i
        
 
         if numel(find( ([datasetSlowSwimBouts( allindex{Fish(i)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish(i)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty cell array for time' num2str(z) 'Fish' num2str(i)]);
            BendAmplitude{z}{Fish(i)}{h}=0;
            MaxAmp_BendAmplitude{z}{Fish(i)}{h}=0;
            medianBout_MaxAmplitude{z}{Fish(i)}=0;
            medianFish_BendAmplitude(z)=0;
          
         else
         idx_TimeWindow= allindex{Fish(i)}(find( ([datasetSlowSwimBouts( allindex{Fish(i)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish(i)}).BoutStart] < (fps*TimeWindow).*z) ) );
            
            for h=1:length(idx_TimeWindow)
                h
            display(['Amplitude currently processing period ' num2str(z)])    
            display(['Amplitude currently processing fish ' num2str(i)])
            display(['Amplitude currently processing bout number ' num2str(h)])
            
            BendAmplitude{z}{Fish(i)}{h}= 57.2958*[datasetSlowSwimBouts(idx_TimeWindow(h)).Bend_Amplitude]; % all Amplitude in bout h, for fishl, in period z
%             for a= 1:length(G2_BendAmplitude{z}{l}{h})
%             if abs(G2_BendAmplitude{z}{l}{h}(a,1))<2;
%                 G2_BendAmplitude{z}{l}{h}(a,1)=2;
%             end
%             end
            %G2_medianAmp_BendAmplitude{z}{l}{h}=median(abs(G2_BendAmplitude{z}{l}{h}));%median Amplitude of bout h, for fishl, in period z
            MaxAmp_BendAmplitude{z}{Fish(i)}{h}=max(abs(BendAmplitude{z}{Fish(i)}{h}));
            end

            
         end
          %G2_medianBout_BendAmplitude{z}{l}= median(cell2mat(G2_medianAmp_BendAmplitude{1,z}{1,l})); % median of all bout for fish l
          medianBout_MaxAmplitude{z}{Fish(i)}= median(cell2mat(MaxAmp_BendAmplitude{z}{Fish(i)})); % median of all(bouts) MaxBendAmp for each Fish
          %medianFish_BendAmplitude(z)= median(cell2mat(medianBout_MaxAmplitude{1,z}));%median of all fish in period Z; So 5 values
    end
end


 for i=1:NumberFish      
       
         
        % output data
         output_Distance(i).Variable=['SlowSwim_Distance'];
         output_Distance(i).Clutch = datasetGoodFish(i).NTrial;
         output_Distance(i).Fish= Fish(i);
         output_Distance(i).FishGeno = datasetGoodFish(i).Genotype;
         output_Distance(i).Period1= medianDistance_time{Fish(i)}(1);
         output_Distance(i).Period2= medianDistance_time{Fish(i)}(2);
         output_Distance(i).Period3= medianDistance_time{Fish(i)}(3);
         output_Distance(i).Period4= medianDistance_time{Fish(i)}(4);
         output_Distance(i).Period5= medianDistance_time{Fish(i)}(5);
         
         
         output_Duration(i).Variable=['SlowSwim_BoutDuration'];
         output_Duration(i).Clutch = datasetGoodFish(i).NTrial;
         output_Duration(i).Fish= Fish(i);
         output_Duration(i).FishGeno = datasetGoodFish(i).Genotype;
         output_Duration(i).Period1= medianBoutDuration_time{Fish(i)}(1);
         output_Duration(i).Period2= medianBoutDuration_time{Fish(i)}(2);
         output_Duration(i).Period3= medianBoutDuration_time{Fish(i)}(3);
         output_Duration(i).Period4= medianBoutDuration_time{Fish(i)}(4);
         output_Duration(i).Period5= medianBoutDuration_time{Fish(i)}(5);
         
         output_Speed(i).Variable=['SlowSwim_Speed'];
         output_Speed(i).Clutch = datasetGoodFish(i).NTrial;
         output_Speed(i).Fish= Fish(i);
         output_Speed(i).FishGeno = datasetGoodFish(i).Genotype;
         output_Speed(i).Period1= medianSpeed_time{Fish(i)}(1);
         output_Speed(i).Period2= medianSpeed_time{Fish(i)}(2);
         output_Speed(i).Period3= medianSpeed_time{Fish(i)}(3);
         output_Speed(i).Period4= medianSpeed_time{Fish(i)}(4);
         output_Speed(i).Period5= medianSpeed_time{Fish(i)}(5);
         
         
         output_NumOfOsc(i).Variable=['SlowSwim_NumberOfOscillations'];
         output_NumOfOsc(i).Clutch = datasetGoodFish(i).NTrial;
         output_NumOfOsc(i).Fish= Fish(i);
         output_NumOfOsc(i).FishGeno = datasetGoodFish(i).Genotype;
         output_NumOfOsc(i).Period1= meanNumberOfOscillations_time{Fish(i)}(1);
         output_NumOfOsc(i).Period2= meanNumberOfOscillations_time{Fish(i)}(2);
         output_NumOfOsc(i).Period3= meanNumberOfOscillations_time{Fish(i)}(3);
         output_NumOfOsc(i).Period4= meanNumberOfOscillations_time{Fish(i)}(4);
         output_NumOfOsc(i).Period5= meanNumberOfOscillations_time{Fish(i)}(5);
         
         
         
         output_TBF(i).Variable=['SlowSwim_TBF'];
         output_TBF(i).Clutch = datasetGoodFish(i).NTrial;
         output_TBF(i).Fish= Fish(i);
         output_TBF(i).FishGeno = datasetGoodFish(i).Genotype;
         output_TBF(i).Period1= medianTBF_time{Fish(i)}(1);
         output_TBF(i).Period2= medianTBF_time{Fish(i)}(2);
         output_TBF(i).Period3= medianTBF_time{Fish(i)}(3);
         output_TBF(i).Period4= medianTBF_time{Fish(i)}(4);
         output_TBF(i).Period5= medianTBF_time{Fish(i)}(5);
   
      
          output_MaxBendAmp(i).Variable=['SlowSwim_MaxAmp_BendAmplitude'];
          output_MaxBendAmp(i).Clutch = datasetGoodFish(i).NTrial;
          output_MaxBendAmp(i).Fish= Fish(i);
          output_MaxBendAmp(i).FishGeno = datasetGoodFish(i).Genotype;
          output_MaxBendAmp(i).Period1= medianBout_MaxAmplitude{1}{Fish(i)};
          output_MaxBendAmp(i).Period2= medianBout_MaxAmplitude{2}{Fish(i)};
          output_MaxBendAmp(i).Period3= medianBout_MaxAmplitude{3}{Fish(i)};
          output_MaxBendAmp(i).Period4= medianBout_MaxAmplitude{4}{Fish(i)};
          output_MaxBendAmp(i).Period5= medianBout_MaxAmplitude{5}{Fish(i)};
   
    end
    
%% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%

output_Distance = struct2table(output_Distance);
writetable(output_Distance);

output_Duration =struct2table(output_Duration);
writetable(output_Duration);

output_Speed=struct2table(output_Speed);
writetable(output_Speed);

output_NumOfOsc= struct2table(output_NumOfOsc);
writetable(output_NumOfOsc);

output_TBF= struct2table(output_TBF);
writetable(output_TBF);

output_MaxBendAmp= struct2table(output_MaxBendAmp);
writetable(output_MaxBendAmp);

%%  Other parameters G2/G0

for z= 1:Period
    
G2_Distance_time=[];
G2_BoutDuration_time=[];
G2_Speed_time=[];
G2_NumberOfOscillations_time=[];
G2_TBF_time=[];
 
    for l=1:length(Fish_G2)

        %BoutDistance
        G2_Distance_time{z}{l}= Distance_time{Fish_G2(l)}{z};
        G2medianDistance_time{z}{l}= median(G2_Distance_time{1,z}{1,l},'omitnan');
        G2medianFishDistance_time(z)= median(cell2mat(G2medianDistance_time{1,z}),'omitnan');
        G2_Distance_SEM(z)=std(cell2mat(G2medianDistance_time{1,z}))/sqrt(length(Fish_G2));
        
        
        %BoutDuration
        G2_BoutDuration_time{z}{l}= BoutDuration_time{Fish_G2(l)}{z}; 
        G2medianBoutDuration_time{z}{l}= median(G2_BoutDuration_time{1,z}{1,l},'omitnan'); % median BD for fish l
        G2medianFishBoutDuration_time(z)= median(cell2mat(G2medianBoutDuration_time{1,z}),'omitnan');%median of all BD of all fish in period Z
        G2_SEM(z)=std(cell2mat(G2medianBoutDuration_time{1,z}))/sqrt(length(Fish_G2));
        
        %Speed
        G2_Speed_time{z}{l}= Speed_time{Fish_G2(l)}{z}; 
        G2medianSpeed_time{z}{l}= median(G2_Speed_time{1,z}{1,l},'omitnan');
        G2medianFishSpeed_time(z)= median(cell2mat(G2medianSpeed_time{1,z}),'omitnan');
        G2_Speed_SEM(z)=std(cell2mat(G2medianSpeed_time{1,z}))/sqrt(length(Fish_G2));
        
        %NumberOfOscillations
        G2_NumberOfOscillations_time{z}{l}= NumberOfOscillations_time{Fish_G2(l)}{z};
        G2meanNumberOfOscillations_time{z}{l}= mean(G2_NumberOfOscillations_time{1,z}{1,l},'omitnan');
        G2meanFishNumberOfOscillations_time(z)= mean(cell2mat(G2meanNumberOfOscillations_time{1,z}),'omitnan');
        G2_NumberOfOscillations_SEM(z)=std(cell2mat(G2meanNumberOfOscillations_time{1,z}))/sqrt(length(Fish_G2));
        
        %TBF
        G2_TBF_time{z}{l}= TBF_time{Fish_G2(l)}{z};
        G2medianTBF_time{z}{l}= median(G2_TBF_time{1,z}{1,l},'omitnan');
        G2medianFishTBF_time(z)= median(cell2mat(G2medianTBF_time{1,z}),'omitnan');
        G2_TBF_SEM(z)=std(cell2mat(G2medianTBF_time{1,z}))/sqrt(length(Fish_G2));
        
    end
    
    
G0_Distance_time=[];
G0_BoutDuration_time=[];
G0_Speed_time=[];
G0_NumberOfOscillations_time=[];
G0_TBF_time=[];
    
    for l=1:length(Fish_G0) 
        l
        %BoutDistance
        G0_Distance_time{z}{l}= Distance_time{Fish_G0(l)}{z};
        G0medianDistance_time{z}{l}= median(G0_Distance_time{1,z}{1,l},'omitnan');
        G0medianFishDistance_time(z)= median(cell2mat(G0medianDistance_time{1,z}),'omitnan');
        G0_Distance_SEM(z)=std(cell2mat(G0medianDistance_time{1,z}))/sqrt(length(Fish_G0));
        
        %BoutDuration
        G0_BoutDuration_time{z}{l}= BoutDuration_time{Fish_G0(l)}{z}; 
        G0medianBoutDuration_time{z}{l}= median(G0_BoutDuration_time{1,z}{1,l},'omitnan'); % median BD for fish l
        G0medianFishBoutDuration_time(z)= median(cell2mat(G0medianBoutDuration_time{1,z}),'omitnan');%median of all BD of all fish in period Z
        G0_SEM(z)=std(cell2mat(G0medianBoutDuration_time{1,z}))/sqrt(length(Fish_G0));
        
        %Speed
        G0_Speed_time{z}{l}= Speed_time{Fish_G0(l)}{z};
        G0medianSpeed_time{z}{l}= median(G0_Speed_time{1,z}{1,l},'omitnan');
        G0medianFishSpeed_time(z)= median(cell2mat(G0medianSpeed_time{1,z}),'omitnan');
        G0_Speed_SEM(z)=std(cell2mat(G0medianSpeed_time{1,z}))/sqrt(length(Fish_G0));
        
        %NumberOfOscillations 
        G0_NumberOfOscillations_time{z}{l}= NumberOfOscillations_time{Fish_G0(l)}{z};
        G0meanNumberOfOscillations_time{z}{l}= mean(G0_NumberOfOscillations_time{1,z}{1,l},'omitnan');
        G0meanFishNumberOfOscillations_time(z)= mean(cell2mat(G0meanNumberOfOscillations_time{1,z}),'omitnan');
        G0_NumberOfOscillations_SEM(z)=std(cell2mat(G0meanNumberOfOscillations_time{1,z}))/sqrt(length(Fish_G0));

        %TBF
        G0_TBF_time{z}{l}= TBF_time{Fish_G0(l)}{z}; 
        G0medianTBF_time{z}{l}= median(G0_TBF_time{1,z}{1,l},'omitnan');
        G0medianFishTBF_time(z)= median(cell2mat(G0medianTBF_time{1,z}),'omitnan');
        G0_TBF_SEM(z)=std(cell2mat(G0medianTBF_time{1,z}))/sqrt(length(Fish_G0));


        
    end
    
    
end % period z

%%  Amplitude G2/G0
G2_MaxAmp_BendAmplitude=[];
 
for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G2)
        %fprintf("-- Fish %d --\n",l);
        i
 
        G2_MaxAmp_BendAmplitude{z}{l}= MaxAmp_BendAmplitude{z}{Fish_G2(l)};
          G2_medianBout_MaxAmplitude{z}{l}= median(cell2mat(G2_MaxAmp_BendAmplitude{1,z}{1,l})); % median of all bout for fish l
          G2_medianFish_BendAmplitude{z}= median(cell2mat(G2_medianBout_MaxAmplitude{1,z}));%median of all fish in period Z
          G2_SEM(z)=std(cell2mat(G2_medianBout_MaxAmplitude{1,z}))/sqrt(length(Fish_G2));
     
    end
    
end;



G0_MaxAmp_BendAmplitude=[];
 
for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G0)
        %fprintf("-- Fish %d --\n",l);
        i
 
        G0_MaxAmp_BendAmplitude{z}{l}= MaxAmp_BendAmplitude{z}{Fish_G0(l)};
          G0_medianBout_MaxAmplitude{z}{l}= median(cell2mat(G0_MaxAmp_BendAmplitude{1,z}{1,l})); % median of all bout for fish l
          G0_medianFish_BendAmplitude{z}= median(cell2mat(G0_medianBout_MaxAmplitude{1,z}));%median of all fish in period Z
          G0_SEM(z)=std(cell2mat(G0_medianBout_MaxAmplitude{1,z}))/sqrt(length(Fish_G0));
     
    end
    
end;


save('workspace_allParameter.mat')



%% All parameters plot

h2=figure(2)

subplot(2,3,1)
title(['BoutDistance']);hold on;

plot(1:(Period), G2medianFishDistance_time,'bo-');hold on;
plot(1:(Period), G0medianFishDistance_time,'ro-');hold on;

errorbar(1:(Period), G2medianFishDistance_time, G2_SEM,'b'); hold on;
errorbar(1:(Period), G0medianFishDistance_time, G0_SEM,'r'); hold on;

xticks(1:5);
xlim ([0 5])
%legend('-/-','+/+');

xlabel("min");
ylabel('BoutDistance (mm)');hold on;

%grid();
hold off;


subplot(2,3,2)
title(['BoutDuration']);hold on;
 
plot(1:(Period), G2medianFishBoutDuration_time,'bo-');hold on;
plot(1:(Period), G0medianFishBoutDuration_time,'ro-');hold on;

errorbar(1:(Period), G2medianFishBoutDuration_time, G2_SEM,'b'); hold on;
errorbar(1:(Period), G0medianFishBoutDuration_time, G0_SEM,'r'); hold on;

xticks(1:5);
xlim ([0 5])

xlabel("min");hold on;
ylabel('BoutDuration (s)');hold on;
%legend('-/-','+/+');
%grid();
hold off;


subplot(2,3,3)
title(['BoutSpeed']);hold on;

plot(1:(Period), G2medianFishSpeed_time,'bo-');hold on;
plot(1:(Period), G0medianFishSpeed_time,'ro-');hold on;
errorbar(1:(Period), G2medianFishSpeed_time, G2_SEM,'b'); hold on;
errorbar(1:(Period), G0medianFishSpeed_time, G0_SEM,'r'); hold on;

xticks(1:5);
xlim ([0 5])
 
xlabel("min");
ylabel('BoutSpeed (mm/s)');hold on;

%legend('-/-','+/+');
%grid();
hold off;


subplot(2,3,4)
title(['Number Of Oscillations']);hold on;
 
plot(1:(Period), G2meanFishNumberOfOscillations_time,'bo-');hold on;
plot(1:(Period), G0meanFishNumberOfOscillations_time,'ro-');hold on;
errorbar(1:(Period), G2meanFishNumberOfOscillations_time, G2_SEM,'b'); hold on;
errorbar(1:(Period), G0meanFishNumberOfOscillations_time, G0_SEM,'r'); hold on;
xticks(1:5);
xlim ([0 5])
 
xlabel("min");hold on;
ylabel('Number Of Oscillations');hold on;
%legend('-/-','+/+');
%grid();
hold off;


subplot(2,3,5)
title(['TBF']);hold on;
 
plot(1:(Period), G2medianFishTBF_time,'bo-');hold on;
plot(1:(Period), G0medianFishTBF_time,'ro-');hold on;

errorbar(1:(Period), G2medianFishTBF_time, G2_SEM,'b'); hold on;
errorbar(1:(Period), G0medianFishTBF_time, G0_SEM,'r'); hold on;
xticks(1:5);
xlim ([0 5])
 
xlabel("min");hold on;
ylabel('TBF (Hz)');hold on;
%legend('-/-','+/+');
%grid();
hold off;


subplot(2,3,6)

title(['Maxi Bend Amplitude']);hold on;
 
plot(1:(Period), cell2mat(G2_medianFish_BendAmplitude),'bo-');hold on; 
plot(1:(Period), cell2mat(G0_medianFish_BendAmplitude),'ro-');hold on;

errorbar(1:(Period), cell2mat(G2_medianFish_BendAmplitude), G2_SEM,'b'); hold on;
errorbar(1:(Period), cell2mat(G0_medianFish_BendAmplitude), G0_SEM,'r'); hold on;
xticks(1:5);
xlim ([0 5])
xlabel("min");
ylabel('Maxi Bend Amplitude (Deg.)');hold on;
%legend('-/-','+/+');
hold off;


saveas(h2,['OverTime6Parameters.fig'])
saveas(h2,['OverTime6Parameters.epsc'])


