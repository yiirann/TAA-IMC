clear all;
clc

% per=0.1;
% [X , T , ind , label , viewNum , clusters] = loaddataset_Yalebala(per);
[X,gt] = loadYale();
% gt=label;

for p =0.1:0.1:0.5
    [X,inds] = bala_incomplete(X,p,gt);
    filename = sprintf("Yale_%.1f_baladata.mat", p);
    save(filename, "X","gt", "inds");
    fprintf("Generated dataset with missing rate p = %.1f\n", p);
end

