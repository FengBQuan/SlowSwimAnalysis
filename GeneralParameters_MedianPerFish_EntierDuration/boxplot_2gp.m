function [] = boxplot_2gp(Gp1,Gp2)
  
%scatter plot
X1= ones(1,length(Gp1));
 
X2 = ones(1,length(Gp2))*2;
 
group = [X1'; X2'];

 
hold on;

sz=10;
scatter(X1, Gp1,sz,'b','filled');

scatter(X2, Gp2,sz,'r','filled');

jitter1;
 
boxplot([Gp1'; Gp2'],group);

hold on;
 
end