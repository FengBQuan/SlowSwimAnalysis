function [WT_FWD,WT_RT,Homo_FWD,Homo_RT] = Plot_ProportionPie(datasetPerBout_FWD1,datasetPerBout_RT1)
%% Proportion of FWD vs RT
nb_wtBouts_FWD1=length(find([datasetPerBout_FWD1(:).Genotype]==2));
nb_homoBouts_FWD1=length(find([datasetPerBout_FWD1(:).Genotype]==0));

nb_wtBouts_RT1=length(find([datasetPerBout_RT1(:).Genotype]==2));
nb_homoBouts_RT1=length(find([datasetPerBout_RT1(:).Genotype]==0));

WT_FWD=nb_wtBouts_FWD1/(nb_wtBouts_RT1+nb_wtBouts_FWD1)
WT_RT=nb_wtBouts_RT1/(nb_wtBouts_RT1+nb_wtBouts_FWD1)

Homo_FWD=nb_homoBouts_FWD1/(nb_homoBouts_RT1+nb_homoBouts_FWD1)
Homo_RT=nb_homoBouts_RT1/(nb_homoBouts_RT1+nb_homoBouts_FWD1)
%%
h1=figure(1);

subplot(1,2,1)
title('wtBout_FWDvsRT');hold on;
X = [nb_wtBouts_FWD1 nb_wtBouts_RT1];hold on
%labels = {'ForwardSwim','RoutineTurn'};
pie(X);hold on
%pie(X,labels)

subplot(1,2,2)
title('HomoBout_FWDvsRT');hold on;
Y = [nb_homoBouts_FWD1 nb_homoBouts_RT1];hold on
pie(Y)

%saveas(h1,'Proportion_FWDvsRT.fig')
end

