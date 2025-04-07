function [X gt] = loadYale()
    load('Yale.mat');
    X{1}=data{1};
    X{2}=data{2};
    X{3}=data{3};
    gt = double(truelabel{1}');
end

