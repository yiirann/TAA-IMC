clear all;
clc


[X,gt] = loadBDGPfea();
X = X';
for p = 0.1:0.1:0.5
    [X,inds] = U_balanced_incomplete(X,p);
    filename = sprintf("BDGP_%.1f_unbaladata.mat", p);
    save(filename, "gt", "inds");
    %save("calall_0.1_unbaladata.mat","gt","inds");
    fprintf("Generated dataset with missing rate p = %.1f\n", p);
end
