function [] = boxplot_4groupes(gp1,gp2,gp3,gp4);


%Plot WT
X_gp1= ones(1,numel(cell2mat(gp1)));
X_gp2= ones(1,numel(cell2mat(gp2)))*2;
X_gp3= ones(1,numel(cell2mat(gp3)))*3;
X_gp4= ones(1,numel(cell2mat(gp4)))*4;

group = [X_gp1'; X_gp2'; X_gp3'; X_gp4'];

% [median_G2_pre,SEM_G2_pre]=grpstats(cell2mat(G2mean_pre),G2groupe_pre,{'median','sem'});
% [median_G2_post,SEM_G2_post]=grpstats(cell2mat(G2mean_post),G2groupe_post,{'median','sem'});

sz=15;
scatter(X_gp1, cell2mat(gp1),sz,'b');hold on;
scatter(X_gp2, cell2mat(gp2),sz,'b','filled');hold on;
scatter(X_gp3, cell2mat(gp3),sz,'r');hold on;%,'filled'
scatter(X_gp4, cell2mat(gp4),sz,'r','filled');hold on;
jitter1;

boxplot([cell2mat(gp1)'; cell2mat(gp2)';cell2mat(gp3)'; cell2mat(gp4)'],group); hold on;%,cell2mat(gp3)'; cell2mat(gp4)'

% Annotations
%grid();
%ylabel('Bout duration(in sec)'); 
%xlabel('Fish');
xlim([0.5 4.5]);
%ylim([0.1 0.5]);
xticks([1 2 3 4]);    
%xticklabels({'+/+','+/-','-/-'});
xticklabels({'+/+','+/+','-/-','-/-'});hold off;


% xlim([0.5 8]);
% %ylim([0.05 0.25]);
% hold on;
% errorbar(1, median_G2_pre, SEM_G2_pre,'k');hold on;
% errorbar(5, median_G2_post, SEM_G2_post,'k');
end