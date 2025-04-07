function [X, ind_folds] = U_balanced_incomplete(X1,r)
V = size(X1,2);

if V==2
    p{1} = 0.2*r;
    p{2} = 1*r;
end
if V==3
    p{1} = 0.2*r;
    p{2} = 1*r;
    p{3} = 1.8*r;
end
if V==4
    p{1} = 0.25*r;
    p{2} = 0.75*r;
    p{3} = 1.25*r;
    p{4} = 1.75*r;
end
if V==5
    p{1} = 0.25*r;
    p{2} = 0.75*r;
    p{3} = 1*r;
    p{4} = 1.25*r;
    p{5} = 1.75*r;
end
if V==6
    p{1} = 0.25*r;
    p{2} = 0.5*r;
    p{3} = 0.75*r;
    p{4} = 1*r;
    p{5} = 1.5*r;
    if r ~=0.5   
        p{6} = 2*r;
    end

end
[X,ind_folds] = splitunbala(X1,p);

end

