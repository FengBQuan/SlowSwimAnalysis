function [DatasetPreEscape_GoodSwimmers,GoodSwimmers] = SelecteGoodSwimmers(DatasetPreEscape,datasetPerFish)
%% selecte good swimmers

Fish = unique([datasetPerFish(:).Condition]);
NumberFish=length(Fish);

Fish_temp=Fish;

for i=1:NumberFish;
    
   % Calculate other parameters by using all index
    allindex{Fish(i)}= find(~([DatasetPreEscape(:).Condition]-Fish(i)));
    
    if length(allindex{Fish(i)})==0 | length(allindex{Fish(i)})<30;
        
        disp(['Bad swimmers for fish' num2str(Fish(i))]);
        
        allindex{Fish(i)}=[];
        Fish_temp(i)=nan;
        
    elseif Fish(i)== 204;
        allindex{Fish(i)}=[];
        Fish_temp(i)=nan;
        disp(['anormal fish' num2str(Fish(i))]);
           
    end      
                    
end

Fish_temp(isnan(Fish_temp))=[];
GoodSwimmers=Fish_temp;
%set dataset PreEscape_GoodSwimmers
DatasetPreEscape_GoodSwimmers= [DatasetPreEscape(cell2mat(allindex))];

end

