%% Clean workspace and load data
close all
clear variables
clc
load('matlab_workspace_SST_4manip_SlowSwim.mat')

%% good swimers
    
% set new datasetSlowSwimBouts

% GoodSwimers=find(~([datasetPerFish(:).MeanIBI]> 10));
% datasetGoodFish=datasetPerFish(GoodSwimers);
datasetGoodFish=datasetPerFish;

Fish = unique([datasetGoodFish(:).Condition]);
NumberFish=length(Fish);


%% find EscapeWindow
EscapeWindow = [30000 30350];


EscapeBout=find((([datasetPerBout(:).BoutStart]> EscapeWindow(1)) & ([datasetPerBout(:).BoutStart]< EscapeWindow(2)) ))   ;
SlowSwimBouts=find(~(([datasetPerBout(:).BoutStart]> EscapeWindow(1)) & ([datasetPerBout(:).BoutStart]< EscapeWindow(2)) ) );
datasetSlowSwimBouts= datasetPerBout(SlowSwimBouts);

%% create outputTables

output_BoutFrequency= struct( 'Variable', [], 'Clutch', [], 'Fish', [], 'FishGeno',[], ...
    'Period1', [], 'Period2',[],'Period3',[],'Period4',[],'Period5',[]);


for i=1:NumberFish;

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

            nFrames= 60285;
            fps= unique([datasetSlowSwimBouts(:).fps]);
            TimeWindow= 120; %in sec
            Period= nFrames/(fps*TimeWindow); %Period = (60034/(100*60);

%% IBI over time
 
IBI_time=[];
medianIBI_time=[];
medianBoutFreq_time=[];
 

    for i=1:NumberFish
   
        for z= 1:Period
    
    display([' BoutFrequency currently processing fish ' num2str(i)])
    display([' BoutFrequency currently processing period ' num2str(z)])
     
        
    
    if  numel(index{Fish(i)})==0;
        IBI_time{Fish(i)}{z}=nan;

        
    else if z==1;
            idx_TimeWindow= index{Fish(i)}(find( ([datasetSlowSwimBouts( allindex{Fish(i)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish(i)}).BoutStart] < (fps*TimeWindow).*z) ) );
            IBI_time{Fish(i)}{z} = [datasetSlowSwimBouts(idx_TimeWindow).InstantaneousIBI];
            
        else
            
            try
                idx_TimeWindow= allindex{Fish(i)}(find( ([datasetSlowSwimBouts( allindex{Fish(i)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish(i)}).BoutStart] < (fps*TimeWindow).*z) ) );
                IBI_time{Fish(i)}{z} = [datasetSlowSwimBouts(idx_TimeWindow).InstantaneousIBI];
                
            catch
                numel(find( ([datasetSlowSwimBouts( allindex{Fish(i)} ).BoutStart] > (fps*TimeWindow).*(z-1)) & ([datasetSlowSwimBouts( allindex{Fish(i)}).BoutStart] < (fps*TimeWindow).*z) ))==0 ;
                disp(['empty cell array all IBIs for time' num2str(z) 'Fish' num2str(i)]);
                IBI_time{Fish(i)}{z}=nan;    % G2_IBI_time: period 1 has 68 fish, fish1 has 127 IBI
            end
            
        end
       
    end
    
 
         medianIBI_time{Fish(i)}{z}= median(IBI_time{1,Fish(i)}{1,z},'omitnan');
         medianBoutFreq_time{Fish(i)}(z)= 1/(medianIBI_time{1,Fish(i)}{1,z});% median of all boutfreq for each Fish

    end % NumberFish i
 
    end
%% output data

    for i=1:NumberFish

        output_BoutFrequency(i).Variable=['SlowSwim_BoutFrequency'];
        output_BoutFrequency(i).Clutch = datasetGoodFish(i).NTrial;
        output_BoutFrequency(i).Fish= Fish(i);
        output_BoutFrequency(i).FishGeno = datasetGoodFish(i).Genotype;
        output_BoutFrequency(i).Period1= medianBoutFreq_time{Fish(i)}(1);
        output_BoutFrequency(i).Period2= medianBoutFreq_time{Fish(i)}(2);
        output_BoutFrequency(i).Period3= medianBoutFreq_time{Fish(i)}(3);
        output_BoutFrequency(i).Period4= medianBoutFreq_time{Fish(i)}(4);
        output_BoutFrequency(i).Period5= medianBoutFreq_time{Fish(i)}(5);
        
    end

output_BoutFrequency=struct2table(output_BoutFrequency);
writetable(output_BoutFrequency);
%% BoutFrequency G2/G0

 for z= 1:Period
     
G2_IBI_time=[];
 
    for l=1:length(Fish_G2)
        %fprintf("-- Fish %d --\n",l);
        
        l
        G2_IBI_time{z}{l}= IBI_time{Fish_G2(l)}{z};
        
        G2medianIBI_time{z}{l}= median(G2_IBI_time{1,z}{1,l},'omitnan');
        G2medianBoutFreq_time{z}{l}=1/(G2medianIBI_time{z}{l});
        
        G2medianFishBoutFreq_time(z)= median(cell2mat(G2medianBoutFreq_time{1,z}),'omitnan');
        
        G2_BoutFreq_SEM(z)=std(cell2mat(G2medianBoutFreq_time{1,z}),'omitnan')/sqrt(length(Fish_G2));



       
    end
    
  
G0_IBI_time=[];
    
    for l=1:length(Fish_G0)
        %fprintf("-- Fish %d --\n",l);
        l
        G0_IBI_time{z}{l}= IBI_time{Fish_G0(l)}{z};
        
        G0medianIBI_time{z}{l}= median(G0_IBI_time{1,z}{1,l},'omitnan');
        G0medianBoutFreq_time{z}{l}=1/(G0medianIBI_time{z}{l});
        
        G0medianFishBoutFreq_time(z)= median(cell2mat(G0medianBoutFreq_time{1,z}),'omitnan');
        
        G0_BoutFreq_SEM(z)=std(cell2mat(G0medianBoutFreq_time{1,z}),'omitnan')/sqrt(length(Fish_G0));
        
    end
    
    
end % period z

save('OverTimeSlowSwim_4manip_BoutFrequency_Workspace.mat')
%% Bout Frequency plot

h1=figure(1)

title ('Bout Frequency over section time');hold on;

plot(1:(Period), G2medianFishBoutFreq_time,'bo-');hold on;     
plot(1:(Period), G0medianFishBoutFreq_time,'ro-');hold on;

legend('+/+','-/-');

errorbar(1:(Period), G2medianFishBoutFreq_time, G2_BoutFreq_SEM,'b'); hold on; 
errorbar(1:(Period), G0medianFishBoutFreq_time, G0_BoutFreq_SEM,'r'); hold on;

xlabel("min");
ylabel('Bout Frequency(Hz)');hold on;

%grid();
hold off;

saveas(h1,['OverTimeSlowSwim_4manip_BoutFrequency.fig'])
saveas(h1,['OverTimeSlowSwim_4manip_BoutFrequency.epsc'])