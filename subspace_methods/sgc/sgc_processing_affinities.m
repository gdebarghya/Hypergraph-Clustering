function [indices,trackSubspace] = sgc_processing_affinities(A,K,trackSubspace,opts)

N = size(A,1);
degrees = A*sum(A,1).';

if opts.normalizeW 
    I0 = (degrees == 0);
    degrees(I0) = 1; %ones(sum(I0),1);
    A = repmat(1./sqrt(degrees),1,size(A,2)).*A;
end

if ~trackSubspace.flag
    [Ui,~,~] = svds(A,K);
    trackSubspace.U = Ui(:,1:K);
    trackSubspace.flag = true;
else
    b = sum(A.^2,1);
    [~,ind] = sort(b,'descend');
    A = A(:,ind);
    for i = 1:size(A,2)
        %% GROUSE update
      
        p = A(:,i);
        w = pinv(trackSubspace.U)*p;
        r = p - trackSubspace.U*w; 
        q = trackSubspace.U*w;
        sig = 2*norm(r)*norm(w); %svds(-2*r*w',1);
        eta = 0.001;
        if norm(r)>0 
            r = r/norm(r); end
        if norm(q)>0 
            q = q/norm(q); end
        if norm(w)>0 
            w = w/norm(w); end
        trackSubspace.U = trackSubspace.U + (sin(sig*eta)*r + (cos(sig*eta)-1)*q)*w';
    end
end        
% trackSubspace.U = orth(trackSubspace.U);

%find initial centers
tol = 1e-8;
seeds = zeros(K,K);
U = trackSubspace.U;

% if opts.normalizeU
%     %U = unitize(U,2);
%     rowNorms = sqrt(sum(U.^2,2));
%     rowNorms(rowNorms==0) = 1;
%     U = U./repmat(rowNorms,1,K);
% end

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

%figure; do_plot_data(U, indicesKmeans);

indices = indicesKmeans;
