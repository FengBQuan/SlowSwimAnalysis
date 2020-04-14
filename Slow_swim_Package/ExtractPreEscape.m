function[DatasetPreEscape]= ExtractPreEscape(datasetPerBout)
%% Extract DataSet PreEscape(5min)

% find EscapeWindow: when you have sevral trials, the EscapeWindow could be differents, to trait all trials at same time, I use the mini of EscapeWindow1, maxi of EscapeWindow2 EscapeWindow2
EscapeWindow = [min(unique([datasetPerBout(:).EscapeWindow1])) max(unique([datasetPerBout(:).EscapeWindow2]))];

DatasetPreEscape=datasetPerBout(find([datasetPerBout(:).BoutEnd]< EscapeWindow(1))); %30037

end

