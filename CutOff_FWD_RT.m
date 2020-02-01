function [DatasetPreEscape_GoodSwimmers_FWD,DatasetPreEscape_GoodSwimmers_RT] = CutOff_FWD_RT(datasetPerBout)

%% add MaxBendAmp in the superstructure

for a=1:length(datasetPerBout)
    datasetPerBout(a).MaxBendAmp= max(abs([datasetPerBout(a).AmplitudeOfAllBends]));
    datasetPerBout(a).MedianBendAmp= median(abs([datasetPerBout(a).AmplitudeOfAllBends]));
end
 
DatasetPreEscape_GoodSwimmers_FWD= datasetPerBout(find([datasetPerBout(:).MaxBendAmp]< 25));   
 
DatasetPreEscape_GoodSwimmers_RT= datasetPerBout(find([datasetPerBout(:).MaxBendAmp]> 25)); 
  

end

