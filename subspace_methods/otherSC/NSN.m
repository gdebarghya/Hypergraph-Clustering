function [Z,U] = NSN(Y,K,kmax,epsilon)
    % Run NSN for the columns of Y
    % 
    % Y : (p*N matrix) Data matrix whose columns are the data points
    % K : The number of neighbors to be collected for each point
    % kmax : The maximum dimension of the subspace
    
    if (nargin < 4)
        epsilon = 1e-4;
    end
    if (nargin < 3)
        kmax = K+1;
    end
    
    [p,N] = size(Y);
    
    % Normalize the data points
    Y = Y./repmat(sqrt(sum(Y.^2,1)),p,1);
    U = cell(1,N);
    
    Z = zeros(N,N);
    for i=1:N
        U{i} = zeros(p,kmax);
        
        % Find K neighbors of Y(:,i)  
        I = [i]; U{i}(:,1) = Y(:,i); P = (U{i}(:,1)'*Y).^2; P(i) = 0;
        for k=1:(kmax-1)
            % Select the closest point to current subspace
            [~,J] = max(P);
            I = [I J];
            
            % Update subspace and projections 
            Unew = Y(:,J)-U{i}(:,1:k)*(U{i}(:,1:k)'*Y(:,J));
            if norm(Unew) > 1e-4
                Unew = Unew/norm(Unew);
                U{i}(:,k+1) = Unew;
                P = P + (Unew'*Y).^2;
            end
            
            P(I) = 0;
        end

        [~,J] = sort(P,2,'descend');
        I = [I J(1:(K-kmax+1))];
        
        Z(I,i) = 1;

        % Collect additional neighbor points lying on the current space  
        if epsilon > 0
            Z(1-P<epsilon,i) = .9;
        end
    end
end