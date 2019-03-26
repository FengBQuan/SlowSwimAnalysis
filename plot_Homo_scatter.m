function [G0groupe_pre,G0groupe_post] = plot_Homo_scatter(G0mean_pre,G0mean_post);
%Plot Homo
G0groupe_pre= ones(1,numel(cell2mat(G0mean_pre)))*3;
G0groupe_post= ones(1,numel(cell2mat(G0mean_post)))*7;
 
[median_G0_pre,SEM_G0_pre]=grpstats(cell2mat(G0mean_pre),G0groupe_pre,{'median','sem'});
[median_G0_post,SEM_G0_post]=grpstats(cell2mat(G0mean_post),G0groupe_post,{'median','sem'});
sz=20;
scatter(G0groupe_pre, cell2mat(G0mean_pre),sz,'r','filled');hold on;
scatter(G0groupe_post, cell2mat(G0mean_post),sz,'r','filled');
jitter1;
xlim([0.5 8]);
%ylim([0.05 0.25]);
hold on;
errorbar(3, median_G0_pre, SEM_G0_pre,'k');hold on;
errorbar(7, median_G0_post, SEM_G0_post,'k');hold off;
end

