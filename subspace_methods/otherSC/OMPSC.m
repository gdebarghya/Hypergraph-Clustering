function Z = OMPSC(Y,K)
    % Run OMP-SSC for the columns of Y.
    % 
    % Y : (p*N matrix) Data matrix whose columns are the data points
    % K : The number of neighbors to be collected for each point
    % (If K = 0, then each point find neighbors until the residual < 1e-4)
    
    [p,N] = size(Y);

    if K == 0
        K = N;
    end
    
    Z = zeros(N,N);
    for i=1:N
        y = Y(:,i);
        r = y;
        I = [];
 
        P = zeros(p,0);
        while(norm(r) > 1e-4 && numel(I) <= K)
            c = abs(r'*Y); c(i) = 0;
            [~,J] = max(c);
            I = [I J];
            P = orth(Y(:,I));
            r = y - P*(P'*y);
        end

        Z(I,i) = abs(Y(:,I)\Y(:,i));
    end

end