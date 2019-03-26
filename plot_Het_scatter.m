function [G1groupe_pre,G1groupe_post] = plot_Het_scatter(G1mean_pre,G1mean_post);

%Plot Het
G1groupe_pre= ones(1,numel(cell2mat(G1mean_pre)))*2;
G1groupe_post= ones(1,numel(cell2mat(G1mean_post)))*6;
 
[median_G1_pre,SEM_G1_pre]=grpstats(cell2mat(G1mean_pre),G1groupe_pre,{'median','sem'});
[median_G1_post,SEM_G1_post]=grpstats(cell2mat(G1mean_post),G1groupe_post,{'median','sem'});
sz=20;
scatter(G1groupe_pre, cell2mat(G1mean_pre),sz,'g','filled');hold on;
scatter(G1groupe_post, cell2mat(G1mean_post),sz,'g','filled');
jitter1;
xlim([0.5 8]);
%ylim([0.05 0.25]);
hold on;
errorbar(2, median_G1_pre, SEM_G1_pre,'k');hold on;
errorbar(6, median_G1_post, SEM_G1_post,'k');

end

