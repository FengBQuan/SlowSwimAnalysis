%% this function can plot scatter of 2 groups, boxplot and mean, errobar
function [] = boxplot_mean(Gp1,Gp2)
  
%scatter and box plot
X1= ones(1,length(Gp1)); 
X2 = ones(1,length(Gp2))*2; 
group = [X1'; X2'];hold on;


sz=18;
scatter(X1, Gp1,sz,'b','filled');
scatter(X2, Gp2,sz,'r','filled');
jitter2;
 
boxplot([Gp1'; Gp2'],group);
hold on;

% mean and errorbar
mean_gp1= mean(Gp1,'omitnan');
SEM_gp1= std(Gp1,'omitnan')/sqrt(numel(Gp1));

mean_gp2= mean(Gp2,'omitnan');
SEM_gp2= std(Gp2,'omitnan')/sqrt(numel(Gp2));

plot(1,mean_gp1,'d','LineWidth',2,'MarkerSize',8,...
'MarkerEdgeColor','k');hold on;
plot(2,mean_gp2,'d','LineWidth',2,'MarkerSize',8,...
'MarkerEdgeColor','k');hold on;


% errorbar(1, mean_gp1, SEM_gp1,'k');hold on;
% errorbar(2, mean_gp2, SEM_gp2,'k');hold on;

xticklabels({'+/+','-/-'});
xlim([0.5 2.5]);
hold off;

end