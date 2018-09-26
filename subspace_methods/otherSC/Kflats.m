function A = Kflats(Y,d,K)
    % Run K-flat for the columns of Y
    % 
    % Y : (p*N matrix) Data matrix whose columns are the data points
    % d : The dimension of the flats
    % K : The number of flats
    
    [p,N] = size(Y);

    Iold = zeros(1,N);
    I    = floor(rand(1,N)*K)+1;
    
    M    = zeros(p,K);
    D    = cell(1,K);
    while(nnz(Iold ~= I))
        Iold = I;

        for i=1:K
            Y2 = Y(:,I==i);
            M(:,i) = mean(Y2,2);
            [U,~,~] = svd(Y2 - repmat(M(:,i),1,size(Y2,2)));
            D{i} = U(:,1:d);
        end
        

        dist = zeros(K,N);
        for i=1:K
            Y2 = Y - repmat(M(:,i),1,N);
            dist(i,:) = sum((Y2 - (D{i}*(D{i}'*Y2))).^2,1);
        end
        
        [~,I] = min(dist,[],1);
    end
    
    A = I;
end