function [NSE, WtErr] = computeNSE(Z,A,K)
    % Compute NSE, the proportion of the points that have an incorrect neighbor
    % (those with different labels). If there are weights on the
    % neighborhood matrix, check only the neighbors with the K maximum weights
    % in magnitude for each point.
    %
    % Z : Neighborhood matrix
    % A : True labels of data points
    % K : Number of effective neighbors (largest weights)
    
    N = size(Z,1);
    is_zero_one_matrix = isequal((Z == 0) | (Z == 1), ones(N,N));
    
    C = 0; Z = Z - diag(diag(Z));
    for i=1:N
        if is_zero_one_matrix
            I = find(Z(:,i));
        else
            [~,I] = sort(abs(Z(:,i)),1,'descend');
            I = I(1:min(K,nnz(Z(:,i))));
        end
        C = C + (sum(A(i) == A(I)) == numel(I));
    end    
    
    NSE = (N-C)/N;
    
    M = zeros(N,N);
    for i=1:N
        M(:,i) = (A ~= A(i));
    end
    WtErr = sum(sum(abs(M.*Z))) / sum(sum(abs(Z)));
end