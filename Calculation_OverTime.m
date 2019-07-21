 %% Clean workspace and load data
close all
clear variables
clc
load('Analysis_results_SST_4manip_SlowSwim.mat')
%% BoutDuration over time calculation
       
            nFrames= 60634;% EscapeWindow(1); using EscapeWindow(1) for first 5 min calculations
            fps= unique([datasetPerBout(:).fps]);
            TimeWindow= 60; %in sec
            Period= nFrames/(fps*TimeWindow); %Period = (210000/(350*60);
% G2---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
 
G2_BoutDuration_time=[];
G2medianBoutDuration_time=[];
G2medianFishBoutDuration_time=[];
 
G2_Total_BoutDuration_time=[];
G2medianFishTotalBD_time=[];
for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G2)
        %fprintf("-- Fish %d --\n",l);
        l
       %while [datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1) & [datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z;
        if z==1;
            idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2_BoutDuration_time{z}{l} = [datasetPerBout(idx_TimeWindow).BoutDuration];
        else
            
         try
            idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2_BoutDuration_time{z}{l} = [datasetPerBout(idx_TimeWindow).BoutDuration];
       
            
         catch
            numel(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty cell array all IBIs for time' num2str(z) 'Fish' num2str(l)]);
            G2_BoutDuration_time{z}{l}=0;
         end;
          
        end;
 
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
G0medianBoutDuration_time=[];
G0medianFishBoutDuration_time=[];
 
G0_Total_BoutDuration_time=[];
G0medianFishTotalBD_time=[];

for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G0)
        %fprintf("-- Fish %d --\n",l);
        l
       %while [datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1) & [datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z;
        if z==1;
            idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G0_BoutDuration_time{z}{l} = [datasetPerBout(idx_TimeWindow).BoutDuration];
        else
            
         try
            idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G0_BoutDuration_time{z}{l} = [datasetPerBout(idx_TimeWindow).BoutDuration];
       
            
         catch
            numel(find( ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty cell array all IBIs for time' num2str(z) 'Fish' num2str(l)]);
            G0_BoutDuration_time{z}{l}=0;
         end;
          
        end;
 
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
% figure(1) BoutDuration over time
h1=figure(1);
title(['BoutDuration over time']);hold on;

plot(1:(Period), G2medianFishBoutDuration_time,'bo-');hold on;
errorbar(1:(Period), G2medianFishBoutDuration_time, G2_SEM,'b'); hold on;

plot(1:(Period), G0medianFishBoutDuration_time,'ro-');hold on;
errorbar(1:(Period), G0medianFishBoutDuration_time, G0_SEM,'r'); hold on;

xlabel("min");hold on;
ylabel('BoutDuration in sec');hold on;
%grid();
hold off;
 
% saveas(h1,['BoutDuration over time.fig'])
% saveas(h1,['BoutDuration over time.png'])
%% figure(2) Total Duration per min
h2=figure(2);
title(['Total Duration per min']);hold on;

plot(1:(Period), G2medianFishTotalBD_time,'bo-');hold on;
errorbar(1:(Period), G2medianFishTotalBD_time, G2_Total_SEM,'b'); hold on;

plot(1:(Period), G0medianFishTotalBD_time,'ro-');hold on;
errorbar(1:(Period), G0medianFishTotalBD_time, G2_Total_SEM,'r'); hold on;

xlabel("min");
ylabel('TotalBoutDuration');hold on;
%grid();
hold off;

% saveas(h2,['Total Duration per min.fig'])
% saveas(h2,['Total Duration per min.png'])


%% Distance over time calculations

            nFrames=60634; %EscapeWindow(1);
            fps= unique([datasetPerBout(:).fps]);
            TimeWindow= 60; %in sec
            Period= nFrames/(fps*TimeWindow); %Period = (210000/(350*60);
% G2---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
 
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
       %while [datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1) & [datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z;
        if z==1;
            idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2_Distance_time{z}{l} = [datasetPerBout(idx_TimeWindow).TotalDistance];
        else
            
         try
            idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2_Distance_time{z}{l} = [datasetPerBout(idx_TimeWindow).TotalDistance];
       
            
         catch
            numel(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
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
       %while [datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1) & [datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z;
        if z==1;
            idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G0_Distance_time{z}{l} = [datasetPerBout(idx_TimeWindow).TotalDistance];
        else
            
         try
            idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G0_Distance_time{z}{l} = [datasetPerBout(idx_TimeWindow).TotalDistance];
       
            
         catch
            numel(find( ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
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

%% figure(3) BoutDistance plot
h3=figure(3);
title(['BoutDistance over time']);hold on;
%G2
plot(1:(Period), G2medianFishDistance_time,'bo-');hold on;
errorbar(1:(Period), G2medianFishDistance_time, G2_SEM,'b'); hold on;
%G0
plot(1:(Period), G0medianFishDistance_time,'ro-');hold on;
errorbar(1:(Period), G0medianFishDistance_time, G0_SEM,'r'); hold on;

xlabel("min");
ylabel('BoutDistance in mm');hold on;
%grid();
hold off;

% saveas(h3,['Bout Duration per min.fig'])
% saveas(h3,['Bout Duration per min.png'])

%% figure(4) Total Distance plot
h4=figure(4);
title(['Total Distance over time']);hold on;
%G2
plot(1:(Period), G2medianFishTotalDistance_time,'bo-');hold on;
errorbar(1:(Period), G2medianFishTotalDistance_time, G2_Total_SEM,'b'); hold on;
%G0
plot(1:(Period), G0medianFishTotalDistance_time,'ro-');hold on;
errorbar(1:(Period), G0medianFishTotalDistance_time, G0_Total_SEM,'r'); hold on;

xlabel("min");
ylabel('Total Distance in mm');hold on;
hold off;

% saveas(h4,['Total Distance per min.fig'])
% saveas(h4,['Total Distance per min.png'])
%% nBout per min calculation       
            nFrames=60634;
            fps=100;
            TimeWindow= 60; %in sec
            Period= nFrames/(fps*TimeWindow); %Period = (210000/(350*60);
% G2---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
G2nBoutPerMin=[];
G2MedianPerMin=[];
for l=1:length(Fish_G2)
    %fprintf("-- Fish %d --\n",l);
    
    for z= 1:Period;
        %fprintf(" %d min %d\n",z);
        
        %find bout position in the time window
        G2nBoutPerMin{z}{l}= length(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] < (fps*TimeWindow).*z) ));
        
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
        G0nBoutPerMin{z}{l}= length(find( ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] < (fps*TimeWindow).*z) ));
        
        G0MedianPerMin(z)= median( cell2mat(G0nBoutPerMin{z}));
        G0_SEM(z)=std(cell2mat(G0nBoutPerMin{z}))/sqrt(length(Fish_G0));
    end;
end;
 

%% figure(5) nBout over section time
h5=figure(5);
title ('median nBout over section time');hold on;

plot(1:Period, G2MedianPerMin,'bo-');hold on;
errorbar(1:Period, G2MedianPerMin, G2_SEM,'b'); hold on;

plot(1:Period, G0MedianPerMin,'ro-');
errorbar(1:Period, G0MedianPerMin, G0_SEM,'r'); hold on;
xlabel("min");
ylabel('nBout');
%grid();
hold off;
legend('-/-','+/+');

% saveas(h5,['Median nBoutPerMin.fig'])
% saveas(h5,['Median nBoutPerMin.png'])  

%% IBI over time

            nFrames=60634;
            fps= unique([datasetPerBout(:).fps]);
            TimeWindow= 60; %in sec
            Period= nFrames/(fps*TimeWindow); %Period = (210000/(350*60);
% G2---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%

G2IBI_time=[];
G2medianIBI_time=[];
G2medianFishIBI_time=[];
for z= 1:Period-1;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G2)
        %fprintf("-- Fish %d --\n",l);
        l
       %while [datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1) & [datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z;
        if z==1;
            idx_TimeWindow= index{Fish_G2(l)}(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2_IBI_time{z}{l} = [datasetPerBout(idx_TimeWindow).InstantaneousIBI];
        else
            
         try
            idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2_IBI_time{z}{l} = [datasetPerBout(idx_TimeWindow).InstantaneousIBI];
       
            
         catch
            numel(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty cell array all IBIs for time' num2str(z) 'Fish' num2str(l)]);
            G2_IBI_time{z}{l}=10;
         end;
          
        end;
 
       %end;
        if isempty(G2_IBI_time{z}{l});
            G2_IBI_time{z}{l}=10;
        end;
       
         G2medianIBI_time{z}{l}= log(median(G2_IBI_time{1,z}{1,l},'omitnan'));
         G2medianFishIBI_time(z)= median(cell2mat(G2medianIBI_time{1,z}),'omitnan');
         G2_SEM(z)=std(cell2mat(G2medianIBI_time{1,z}))/sqrt(length(Fish_G2));
    end;
    
end;
 
% G0---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%

G0IBI_time=[];
G0medianIBI_time=[];
G0medianFishIBI_time=[];
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
        
       G0medianIBI_time{z}{l}= log(median(G0IBI_time{1,z}{1,l},'omitnan'));
       G0medianFishIBI_time(z)= median(cell2mat(G0medianIBI_time{1,z}),'omitnan');
       G0_SEM(z)=std(cell2mat(G0medianIBI_time{1,z}))/sqrt(length(Fish_G0));
    end;
end;

% h6=figure(6) IBI over time
h6=figure(6);
title ('IBI over section time');hold on;

plot(1:(Period-1), G2medianFishIBI_time,'bo-');hold on;
%errorbar(1:(Period-1), G2medianFishIBI_time, G2_SEM,'b'); hold on;
       
plot(1:(Period-1), G0medianFishIBI_time,'ro-');hold on;
%errorbar(1:(Period-1), G0medianFishIBI_time, G0_SEM,'r'); hold on;
xlabel("min");
ylabel('logIBI');hold on;
%grid();
hold off;
 
% saveas(h6,['IBIs per 10s.fig'])
% saveas(h6,['IBIs per 10s.png'])

%% Amplitude per min       
            nFrames= 60634;% EscapeWindow(1); %
            fps= unique([datasetPerBout(:).fps]);
            TimeWindow= 60; %in sec
            Period= nFrames/(fps*TimeWindow); %Period = (210000/(350*60);
% G2---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
 
G2_BendAmplitude=[];
G2_medianAmp_BendAmplitude=[];
G2_medianBout_BendAmplitude=[];
G2_medianFish_BendAmplitude=[];
 
 
for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G2)
        %fprintf("-- Fish %d --\n",l);
        l
 
         if numel(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty cell array all IBIs for time' num2str(z) 'Fish' num2str(l)]);
            G2_BendAmplitude{z}{l}{h}=0;
            G2_medianAmp_BendAmplitude{z}{l}{h}=0;
            G2_medianBout_BendAmplitude{z}{l}=0;
            G2_medianFish_BendAmplitude{z}=0;
          
         else
         idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            
            for h=1:length(idx_TimeWindow)
                h
            display(['currently processing period ' num2str(z)])    
            display(['currently processing fish ' num2str(l)])
            display(['currently processing bout number ' num2str(h)])
            
            G2_BendAmplitude{z}{l}{h}= 57.2958*[datasetPerBout(idx_TimeWindow(h)).Bend_Amplitude]; % all Amplitude in bout h, for fishl, in period z
            G2_medianAmp_BendAmplitude{z}{l}{h}=max(abs(G2_BendAmplitude{z}{l}{h}));%median Amplitude of bout h, for fishl, in period z
            end
            
         end
          G2_medianBout_BendAmplitude{z}{l}= median(cell2mat(G2_medianAmp_BendAmplitude{1,z}{1,l})); % median of all bout for fish l
          G2_medianFish_BendAmplitude{z}= median(cell2mat(G2_medianBout_BendAmplitude{1,z}));%median of all fish in period Z
          G2_SEM(z)=std(cell2mat(G2_medianBout_BendAmplitude{1,z}))/sqrt(length(Fish_G2));
         
        
    end
    
end;
 
% G0---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
G0_BendAmplitude=[];
G0_medianAmp_BendAmplitude=[];
G0_medianBout_BendAmplitude=[];
G0_medianFish_BendAmplitude=[];
 
 
for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G0)
        %fprintf("-- Fish %d --\n",l);
        l
 
         if numel(find( ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty Bout for time' num2str(z) 'Fish' num2str(l)]);
            G0_BendAmplitude{z}{l}{h}=0;
            G0_medianAmp_BendAmplitude{z}{l}{h}=0;
            G0_medianBout_BendAmplitude{z}{l}=0;
            G0_medianFish_BendAmplitude{z}=0;
          
         else
         idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            
            for h=1:length(idx_TimeWindow)
                h
            display(['currently processing period ' num2str(z)])    
            display(['currently processing fish ' num2str(l)])
            display(['currently processing bout number ' num2str(h)])
            
            G0_BendAmplitude{z}{l}{h}= 57.2958*[datasetPerBout(idx_TimeWindow(h)).Bend_Amplitude];% all Amplitude in bout h, for fishl, in period z
            G0_medianAmp_BendAmplitude{z}{l}{h}=max(abs(G0_BendAmplitude{z}{l}{h}));%median Amplitude of bout h, for fishl, in period z
            end
 
         end
          G0_medianBout_BendAmplitude{z}{l}= median(cell2mat(G0_medianAmp_BendAmplitude{1,z}{1,l})); % median of all bout for fish l
          G0_medianFish_BendAmplitude{z}= median(cell2mat(G0_medianBout_BendAmplitude{1,z}));%median of all fish in period Z
          G0_SEM(z)=std(cell2mat(G0_medianBout_BendAmplitude{1,z}))/sqrt(length(Fish_G0));
     
    end
    
end;
 
%  figure(7) Amplitude
h7=figure(7);
title(['Median Bend_Amplitude per min']);hold on;
 
plot(1:(Period), cell2mat(G2_medianFish_BendAmplitude),'bo-');hold on;
errorbar(1:(Period), cell2mat(G2_medianFish_BendAmplitude), G2_SEM,'b'); hold on;
 
plot(1:(Period), cell2mat(G0_medianFish_BendAmplitude),'ro-');hold on;
errorbar(1:(Period), cell2mat(G0_medianFish_BendAmplitude), G0_SEM,'r'); hold on;
xlabel("min");
ylabel('Degree');hold on;
 
hold off;
 
% saveas(h7,['Bend_Amplitude per min.fig'])
% saveas(h7,['Bend_Amplitude per min.png'])


%% InstantaneousTBF over time 
   
            nFrames= 60634;% EscapeWindow(1); %
            fps= unique([datasetPerBout(:).fps]);
            TimeWindow= 60; %in sec
            Period= nFrames/(fps*TimeWindow); %Period = (210000/(350*60);
% G2---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
 
G2_TBF=[];
G2_meanAmp_TBF=[];
G2_meanBout_TBF=[];
G2_meanFish_TBF=[]; 
 
for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G2)
        %fprintf("-- Fish %d --\n",l);
        l
 
         if numel(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty Bout for time' num2str(z) 'Fish' num2str(l)]);
            G2_TBF{z}{l}{h}=0;
            G2_meanAmp_TBF{z}{l}{h}=0;
            G2_meanBout_TBF{z}{l}=0;
            G2_meanFish_TBF{z}=0;
          
         else
         idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            
            for h=1:length(idx_TimeWindow)
                h
            display(['currently processing period ' num2str(z)])    
            display(['currently processing fish ' num2str(l)])
            display(['currently processing bout number ' num2str(h)])
            
            G2_TBF{z}{l}{h}= [datasetPerBout(idx_TimeWindow(h)).InstantaneousTBF]; % all Amplitude in bout h, for fishl, in period z
            G2_meanAmp_TBF{z}{l}{h}=mean(abs(G2_TBF{z}{l}{h}));%mean Amplitude of bout h, for fishl, in period z
            end
         
          
         end
 
       
         
          G2_meanBout_TBF{z}{l}= mean(cell2mat(G2_meanAmp_TBF{1,z}{1,l})); % mean of all bout for fish l
          G2_meanFish_TBF{z}= mean(cell2mat(G2_meanBout_TBF{1,z}));%mean of all fish in period Z
          G2_SEM(z)=std(cell2mat(G2_meanBout_TBF{1,z}))/sqrt(length(Fish_G2));
         
        
    end
    
end;
% G1---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
 
% G0---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
G0_TBF=[];
G0_meanAmp_TBF=[];
G0_meanBout_TBF=[];
G0_meanFish_TBF=[];
 
 
for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G0)
        %fprintf("-- Fish %d --\n",l);
        l
 
         if numel(find( ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty Bout for time' num2str(z) 'Fish' num2str(l)]);
            G0_TBF{z}{l}{h}=0;
            G0_meanAmp_TBF{z}{l}{h}=0;
            G0_meanBout_TBF{z}{l}=0;
            G0_meanFish_TBF{z}=0;
          
         else
         idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            
            for h=1:length(idx_TimeWindow)
                h
            display(['currently processing period ' num2str(z)])    
            display(['currently processing fish ' num2str(l)])
            display(['currently processing bout number ' num2str(h)])
            
            G0_TBF{z}{l}{h}= [datasetPerBout(idx_TimeWindow(h)).InstantaneousTBF]; % all Amplitude in bout h, for fishl, in period z
            G0_meanAmp_TBF{z}{l}{h}=mean(abs(G0_TBF{z}{l}{h}));%mean Amplitude of bout h, for fishl, in period z
            end

         end
 
          G0_meanBout_TBF{z}{l}= mean(cell2mat(G0_meanAmp_TBF{1,z}{1,l})); % mean of all bout for fish l
          G0_meanFish_TBF{z}= mean(cell2mat(G0_meanBout_TBF{1,z}));%mean of all fish in period Z
          G0_SEM(z)=std(cell2mat(G0_meanBout_TBF{1,z}))/sqrt(length(Fish_G0));
         
         
    end
    
end;

%% h8=figure(8) mean TBF
h8=figure(8);
title(['Mean TBF per min']);hold on;
%plot
plot(1:(Period), cell2mat(G2_meanFish_TBF),'bo-');hold on;
errorbar(1:(Period), cell2mat(G2_meanFish_TBF), G2_SEM,'b'); hold on;
 
plot(1:(Period), cell2mat(G0_meanFish_TBF),'ro-');hold on;
errorbar(1:(Period), cell2mat(G0_meanFish_TBF), G0_SEM,'r'); hold on;
xlabel("min");
ylabel('Hz');hold on;
hold off;
c% saveas(h8,['Mean TBF per min.fig'])
% saveas(h8,['Mean TBF per min.png'])

%%  Median NumberOfOscillations over time calculation
       
            nFrames= 60634;% EscapeWindow(1); using EscapeWindow(1) for first 5 min calculations
            fps= unique([datasetPerBout(:).fps]);
            TimeWindow= 60; %in sec
            Period= nFrames/(fps*TimeWindow); %Period = (210000/(350*60);
% G2---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
 
G2_NumberOfOscillations_time=[];
G2medianNumberOfOscillations_time=[];
G2medianFishNumberOfOscillations_time=[];
 
for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G2)
        %fprintf("-- Fish %d --\n",l);
        l
       %while [datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1) & [datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z;
        if z==1;
            idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2_NumberOfOscillations_time{z}{l} = [datasetPerBout(idx_TimeWindow).NumberOfOscillations];
        else
            
         try
            idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2_NumberOfOscillations_time{z}{l} = [datasetPerBout(idx_TimeWindow).NumberOfOscillations];
       
            
         catch
            numel(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
            disp(['empty cell array all IBIs for time' num2str(z) 'Fish' num2str(l)]);
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
       %while [datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1) & [datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z;
        if z==1;
            idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G0_NumberOfOscillations_time{z}{l} = [datasetPerBout(idx_TimeWindow).NumberOfOscillations];
        else
            
         try
            idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G0_NumberOfOscillations_time{z}{l} = [datasetPerBout(idx_TimeWindow).NumberOfOscillations];
       
            
         catch
            numel(find( ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
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
% figure(1) NumberOfOscillations over time
h1=figure(1);
title(['NumberOfOscillations over time']);hold on;
 
plot(1:(Period), G2medianFishNumberOfOscillations_time,'bo-');hold on;
errorbar(1:(Period), G2medianFishNumberOfOscillations_time, G2_SEM,'b'); hold on;
 
plot(1:(Period), G0medianFishNumberOfOscillations_time,'ro-');hold on;
errorbar(1:(Period), G0medianFishNumberOfOscillations_time, G0_SEM,'r'); hold on;
 
xlabel("min");hold on;
ylabel('NumberOfOscillations in sec');hold on;
%grid();
hold off;
 
% saveas(h1,['NumberOfOscillations over time.fig'])
% saveas(h1,['NumberOfOscillations over time.png'])

%% Mean NumberOfOscillations over time calculation
       
            nFrames= 60634;% EscapeWindow(1); using EscapeWindow(1) for first 5 min calculations
            fps= unique([datasetPerBout(:).fps]);
            TimeWindow= 60; %in sec
            Period= nFrames/(fps*TimeWindow); %Period = (210000/(350*60);
% G2---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
 
G2_NumberOfOscillations_time=[];
G2meanNumberOfOscillations_time=[];
G2meanFishNumberOfOscillations_time=[];

for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G2)
        %fprintf("-- Fish %d --\n",l);
        l
       %while [datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1) & [datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z;
        if z==1;
            idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2_NumberOfOscillations_time{z}{l} = [datasetPerBout(idx_TimeWindow).NumberOfOscillations];
        else
            
         try
            idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2_NumberOfOscillations_time{z}{l} = [datasetPerBout(idx_TimeWindow).NumberOfOscillations];
       
            
         catch
            numel(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
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
       %while [datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1) & [datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z;
        if z==1;
            idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G0_NumberOfOscillations_time{z}{l} = [datasetPerBout(idx_TimeWindow).NumberOfOscillations];
        else
            
         try
            idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G0_NumberOfOscillations_time{z}{l} = [datasetPerBout(idx_TimeWindow).NumberOfOscillations];
       
            
         catch
            numel(find( ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
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
% figure(1) NumberOfOscillations over time
h1=figure(1);
title(['NumberOfOscillations over time']);hold on;
 
plot(1:(Period), G2meanFishNumberOfOscillations_time,'bo-');hold on;
errorbar(1:(Period), G2meanFishNumberOfOscillations_time, G2_SEM,'b'); hold on;
 
plot(1:(Period), G0meanFishNumberOfOscillations_time,'ro-');hold on;
errorbar(1:(Period), G0meanFishNumberOfOscillations_time, G0_SEM,'r'); hold on;
 
xlabel("min");hold on;
ylabel('NumberOfOscillations in sec');hold on;
%grid();
hold off;
 
% saveas(h1,['NumberOfOscillations over time.fig'])
% saveas(h1,['NumberOfOscillations over time.png'])


%% TBF(NumberOfOscillations/Duration) over time calculation
       
            nFrames= 60634;% EscapeWindow(1); using EscapeWindow(1) for first 5 min calculations
            fps= unique([datasetPerBout(:).fps]);
            TimeWindow= 60; %in sec
            Period= nFrames/(fps*TimeWindow); %Period = (210000/(350*60);
% G2---------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
 
G2_TBF_time=[];
G2medianTBF_time=[];
G2medianFishTBF_time=[];
 
for z= 1:Period;
    z %fprintf(" %d min %d\n",z);
    for l=1:length(Fish_G2)
        %fprintf("-- Fish %d --\n",l);
        l
       %while [datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1) & [datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z;
        if z==1;
            idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2_TBF_time{z}{l} = [datasetPerBout(idx_TimeWindow).NumberOfOscillations]/[datasetPerBout(idx_TimeWindow).BoutDuration];
        else
          
         try
            idx_TimeWindow= allindex{Fish_G2(l)}(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G2_TBF_time{z}{l} = [datasetPerBout(idx_TimeWindow).NumberOfOscillations]/[datasetPerBout(idx_TimeWindow).BoutDuration];
       
            
         catch
            numel(find( ([datasetPerBout( allindex{Fish_G2(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G2(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
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
       %while [datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1) & [datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z;
        if z==1;
            idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G0_TBF_time{z}{l} = [datasetPerBout(idx_TimeWindow).NumberOfOscillations]/[datasetPerBout(idx_TimeWindow).BoutDuration];
        else
            
         try
            idx_TimeWindow= allindex{Fish_G0(l)}(find( ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ) );
            G0_TBF_time{z}{l} = [datasetPerBout(idx_TimeWindow).NumberOfOscillations]/[datasetPerBout(idx_TimeWindow).BoutDuration];
       
            
         catch
            numel(find( ([datasetPerBout( allindex{Fish_G0(l)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetPerBout( allindex{Fish_G0(l)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
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
%% figure(1) TBF over time
h1=figure(1);
title(['TBF over time']);hold on;
 
plot(1:(Period), G2medianFishTBF_time,'bo-');hold on;
errorbar(1:(Period), G2medianFishTBF_time, G2_SEM,'b'); hold on;
 
plot(1:(Period), G0medianFishTBF_time,'ro-');hold on;
errorbar(1:(Period), G0medianFishTBF_time, G0_SEM,'r'); hold on;
 
xlabel("min");hold on;
ylabel('TBF in Hz');hold on;
%grid();
hold off;


 













