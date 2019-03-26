function [G2groupe_pre,G2groupe_post] = plot_WT_scatter(G2mean_pre,G2mean_post);


%Plot WT
G2groupe_pre= ones(1,numel(cell2mat(G2mean_pre)));
G2groupe_post= ones(1,numel(cell2mat(G2mean_post)))*5; 
 
[median_G2_pre,SEM_G2_pre]=grpstats(cell2mat(G2mean_pre),G2groupe_pre,{'median','sem'});
[median_G2_post,SEM_G2_post]=grpstats(cell2mat(G2mean_post),G2groupe_post,{'median','sem'});
%[median_G2_Ratio,SEM_G2_Ratio]=grpstats(G2_Ratio,G2groupe_post,{'median','sem'});
sz=20;
scatter(G2groupe_pre, cell2mat(G2mean_pre),sz,'b','filled');hold on;
scatter(G2groupe_post, cell2mat(G2mean_post),sz,'b','filled');
jitter1;
xlim([0.5 8]);
%ylim([0.05 0.25]);
hold on;
errorbar(1, median_G2_pre, SEM_G2_pre,'k');hold on;
errorbar(5, median_G2_post, SEM_G2_post,'k');
end

