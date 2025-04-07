function [UU,obj] = update_TAAIMC(X,truthF,numanchor,ind,lambda1,lambda2)
%%% initialize
epson = 1e-6;
m = numanchor;
numclass = length(unique(truthF));
numview = length(X);
numsample = size(truthF,1);
n = numsample;
alpha = ones(1,numview)/numview;
P = cell(1,numview);   % m  * m
for i = 1:numview
    dv = size (X{i},1);
    A{i} = zeros(dv,m);         % d  * m
    Z{i} = zeros(m,n);
    Z{i}(:,1:m) = eye(m);
    P{1,i} = eye(m,m);
    J{i} = zeros(m,numsample); % m  * n
    Y{i} = zeros(m,numsample); % m  * n
end
mu = 1e-5;
maxmu = 1e10;
rho = 1.5;
missingindex = constructA(ind);% miss:1x6cell 1x2386 存在样本位置为1
tao = zeros(1,numview);
%%
sX = [m, numsample, numview];
converge_Z_J=[];
Isconverg = 0;
iter = 0;
%%
while(Isconverg == 0)
    %% optimize A
    for iv = 1:numview
        part1 = X{iv} * Z{iv}';
        [Unew,~,Vnew] = svd(part1,'econ');
        A{iv} = Unew*Vnew';
    end
    
    %% optimize Zv
    for iv=1:numview
        C1 = alpha(iv) * ind(:,iv)'+lambda2+0.5*mu;%1x1474
        C2 = alpha(iv) * A{iv}' * X{iv} + lambda2 * P{iv}' * Z{1} + 0.5*mu*P{iv}' * (J{iv} - Y{iv}/mu);%7x1474
        for ii=1:numsample
            idx = 1:numanchor;
            ut = C2(idx,ii)./(C1(ii));  %Z是mxn
            Z{iv}(idx,ii) = EProjSimplex_new(ut);
        end
    end
    clear C1 C2
    %% optimize 排序矩阵P
    for iv=1:numview
        temp1 = 2*lambda2*Z{iv}*Z{1}'+mu*Z{iv}*(J{iv} - Y{iv}/mu)';
        P{1,iv} = shiyishi(temp1);
    end
    clear temp1
    %% update alpha
    for iv =1:numview
        tao(iv) = norm(X{iv} - A{iv} * (Z{iv}.*repmat(missingindex{iv},m,1)),'fro');
        alpha(iv)=sqrt(1/tao(iv));
    end
    %% optimize J
    for i=1:numview
        ZP{i} = P{i}*Z{i};
    end
    ZP_tensor = cat(3, ZP{:,:});
    Y_tensor = cat(3, Y{:,:});
    z = ZP_tensor(:);
    y = Y_tensor(:);
    [j, ~] = wshrinkObj(z+1/mu*y, lambda1/mu,sX,0,3);
    J_tensor = reshape(j, sX);
    for i=1:numview
        J{i} = J_tensor(:,:,i);
    end
    %% optimize Y
    for i=1:numview
        Y{i} = Y{i} + mu * (P{i}*Z{i} - J{i});
    end
    
    max_Z_J=0;
    Isconverg = 1;
    for i = 1:numview
        if (norm(P{i}*Z{i} - J{i}, inf) > epson)
            history.norm_Z_J = norm(P{i}*Z{i} - J{i}, inf);
            Isconverg = 0;
            max_Z_J = max(max_Z_J, history.norm_Z_J);
        end
        
    end
    converge_Z_J=[converge_Z_J max_Z_J];
    %---------- Update mu--------------%
    mu = min(mu*rho, maxmu);
    iter = iter + 1;
    if (iter==100)
        Isconverg = 1;
    end
    obj(iter) = max_Z_J;
end

affinity1 = zeros(m,numsample);
for i = 1:numview
    affinity1 = affinity1 + ZP{i};
end
Z1 = affinity1/numview;

[UU,~,V]=svd(Z1','econ');
UU = UU(:,1:numclass);