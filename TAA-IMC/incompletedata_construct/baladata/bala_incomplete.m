function  [X ind] = bala_incomplete(X,p,gt)
truth = gt;
numsample = size(truth,1);
viewNum = length(X); %view number  
%clusters = length(unique(truth));     %cluster number
T = cell(viewNum,1);
%label=truth;
ind = cell(1,10);
for jj = 1:10
   ind{jj} = ones(numsample, viewNum);
   %remove some instances to obtain incomplete data
   ind{jj} = splitDigitData(ind{jj}, p, 1 );
end
   f = randi(10);
   ind = ind{f}; 
for j = 1 : viewNum
    XX = X{j}';
    XX = NormalizeFea(XX,1);
    ind_0 = find(ind(:,j) == 0);  
    XX(ind_0,:) = 0; %缺失样本处为0
    Y{j} = XX';
end
X = Y;

end