function [X gt] = loadBDGPfea()
    load("BDGP_fea.mat")
    for i = 1:length(X);
        X{i} = X{i}';
    end
    gt = Y;
end

