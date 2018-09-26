function indices = tetris_processing_affinities(A,columnInds,K,opts)

N = size(A,1);

indices = zeros(N,1);
numCols = size(A,2);
tenOrder = size(columnInds,1) + 1;

W = zeros(N);
for i = 1:numCols
    W(:,columnInds(:,i)) = W(:,columnInds(:,i)) + repmat(A(:,i),1,tenOrder-1);
end
degrees = sum(W,2);

index_sort = 1:N;

if opts.alpha>0

    [~, index_sort] = sort(degrees, 'ascend');
    if opts.alpha < 1 % percentage
        opts.alpha = ceil(N*opts.alpha);
    end

    W = W(index_sort(opts.alpha+1:N),index_sort(opts.alpha+1:N));
    degrees = degrees(index_sort(opts.alpha+1:N),1); 
    N = N - opts.alpha;

end

if opts.normalizeW
  
    I0 = (degrees == 0);
    degrees(I0) = 1; %ones(sum(I0),1);
    W = W./repmat(degrees,1,N); 

end

[Ui,~,~,flag] = svds(W,K);
if ((flag) || size(Ui,2)<K)
    [Ui,Si,~] = svd(W);
    [~,ISi] = sort(diag(Si),'descend');
    Ui = Ui(:,ISi(1:K));
end 

U = Ui(:,1:K);

if opts.normalizeU
    %U = unitize(U,2);
    rowNorms = sqrt(sum(U.^2,2));
    rowNorms(rowNorms==0) = 1;
    U = U./repmat(rowNorms,1,K);
end

%find initial centers
tol = 1e-8;
seeds = zeros(K,K);

switch opts.seedType
    
    case 'hard'

        u0 = mean(U,1);
        [um,ind_m] = max(sum((U-repmat(u0,N,1)).^2,2));
        seeds(1,:) = U(ind_m(1),:);

        k = 1;
        U1 = U(sum((U - repmat(seeds(1,:),N,1)).^2,2)>tol,:);
        while k < K && size(U1,1)>0
            [um,ind_m] = max(sum((repmat(U1,1,k)-repmat(reshape(seeds(1:k,:)',[],1)',size(U1,1),1)).^2,2));
            %[um,ind_m] = min(max(U1*seeds(1:k,:)',[],2));
            k = k+1;
            seeds(k,:) = U1(ind_m(1),:);
            U1 = U1(sum((U1 - repmat(seeds(k,:),size(U1,1),1)).^2,2)>tol,:);
        end
        
    case 'soft'
        
        u0 = mean(U,1);
        w = sum((U-repmat(u0,N,1)).^2,2);
        w = w/sum(w);
        seeds(1,:) = U(randsample(N,1,true,w),:);

        k = 1;
        U1 = U(sum((U - repmat(seeds(1,:),N,1)).^2,2)>tol,:);
        while k < K && size(U1,1)>0
            w = sum((repmat(U1,1,k)-repmat(reshape(seeds(1:k,:)',[],1)',size(U1,1),1)).^2,2);
            w = w/sum(w);
            %[um,ind_m] = min(max(U1*seeds(1:k,:)',[],2));
            k = k+1;
            seeds(k,:) = U1(randsample(size(U1,1),1,true,w),:);
            U1 = U1(sum((U1 - repmat(seeds(k,:),size(U1,1),1)).^2,2)>tol,:);
        end
end

if k<K
    indicesKmeans = ceil(K*rand(N,1));
else
    indicesKmeans = kmeans(U,K,'start',seeds,'EmptyAction','drop');
end

indices(index_sort(opts.alpha+1:end),1) = indicesKmeans;
