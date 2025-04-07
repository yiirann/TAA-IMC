function [X_miss, splitInd] = splitunbala(X,P)
M = length(X);
X1 = X{1};
N = size(X1,2);
splitInd = ones(N,M);
indCell = cell(M,1);
for i=1:M
p = P{i};
delNum{i} = round(p * N);
indCell{i} = randperm(N);
splitInd(indCell{i}(1:delNum{i}),i) = 0;%在第i个视图，使缺失的样本为0
counter(i) = delNum{i} + 1;
end
while 1
    zerosInd = find(sum(splitInd,2)==0);
    if size(zerosInd,1) == 0 
        break;
    else
        i = randi(M);
        splitInd(zerosInd(1),i) = 1;     
        splitInd(indCell{i}(counter(i)),i) = 0; 
        counter(i) = counter(i) + 1;
    end        
end
for j = 1 : M
    XX = X{j}';
    XX = NormalizeFea(XX,1);
    ind_0 = find(splitInd(:,j) == 0);  
    XX(ind_0,:) = 0; %缺失样本处为0
    Y{j} = XX';
end
X_miss = Y;
end