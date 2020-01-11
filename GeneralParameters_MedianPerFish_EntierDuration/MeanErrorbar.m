function [] = MeanErrorbar(gp1,gp2)

%Plot WT
X_gp1= ones(1,numel(gp1));
X_gp2= ones(1,numel(gp2))*2;

group = [X_gp1'; X_gp2'];

% [median_G2_pre,SEM_G2_pre]=grpstats(cell2mat(G2mean_pre),G2groupe_pre,{'median','sem'});
% [median_G2_post,SEM_G2_post]=grpstats(cell2mat(G2mean_post),G2groupe_post,{'median','sem'});
mean_gp1= mean(gp1,'omitnan');
SEM_gp1= std(gp1,'omitnan')/sqrt(numel(gp1));

mean_gp2= mean(gp2,'omitnan');
SEM_gp2= std(gp2,'omitnan')/sqrt(numel(gp2));

%%%%%%To delete if use with boxplot
sz=20;
scatter(X_gp1, gp1,sz,'b','filled');hold on; % ,'filled'
scatter(X_gp2, gp2,sz,'r','filled');hold on; % 
jitter1;hold on;
%%%%%%%

plot(1,mean_gp1,'d','LineWidth',2,'MarkerSize',8,...
'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5]);hold on;
plot(2,mean_gp2,'d','LineWidth',2,'MarkerSize',8,...
'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5]);hold on;

errorbar(1, mean_gp1, SEM_gp1,'k');hold on;
errorbar(2, mean_gp2, SEM_gp2,'k');hold on;


hold off;

end

