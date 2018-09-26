% Code by Debarghya Ghoshdastidar for the paper arXiv:1606.06516

clear all; close all; clc;

n_settings = 16;

K_exp = 3;               % # partitions
n_exp = repmat(10:10:40,1,4);             % # vertices in each partition
p_exp = reshape(repmat(0.25:0.05:0.4,4,1),1,n_settings);                  
q_exp = 0.2*ones(1,n_settings);                      

n_algo = 6;          
n_trial = 50;               

E_exp = zeros(n_algo,n_settings,n_trial);     % Clustering error
tStart = tic;
for i_trial = 1:n_trial          % # trials
for i_exp = 1:n_settings
    
    disp(' ')
    disp(['Trial ' int2str(i_trial) '/' int2str(n_trial) ' : Setting ' int2str(i_exp) '/' int2str(n_settings)])
    
    K = K_exp;
    p = p_exp(i_exp);
    q = q_exp(i_exp);
    n = n_exp(i_exp)*ones(1,K);  
    N = sum(n);                  % Total # of samples
    
    %% Adjacency tensor and incidence matrix
    trueLabel = zeros(1,N);           % True labels for sample points
    count = 0;
    for i = 1:K
        trueLabel((count+1):(count+n(i))) = i;
        count = count+n(i);
    end
    
    A = zeros(N,N,N);
    edg = zeros(nchoosek(N,3),3);
    edgwt = ones(nchoosek(N,3),1);
    edgcount = 0;
    
    for i1 = 1:N
        if mod(i1,50)==0
            disp(int2str(i1)) 
        end
        for i2 = (i1+1):N
            for i3 = (i2+1):N
                if (min(trueLabel([i1 i2 i3]))==max(trueLabel([i1 i2 i3])))
                    prob = p;
                else
                    prob = q;
                end
                if (rand(1)>prob)
                    continue
                end
                
                A(i1,i2,i3) = 1; A(i1,i3,i2) = 1;
                A(i2,i1,i3) = 1; A(i2,i3,i1) = 1;
                A(i3,i2,i1) = 1; A(i3,i1,i2) = 1;
                
                edgcount = edgcount + 1;
                edg(edgcount,:) = [i1 i2 i3];
            end
        end
    end
    
    edg = edg(1:edgcount,:);
    edgwt = edgwt(1:edgcount,1);
    disp('Hypergraph constructed...')
    
    %% TTM
    display('TTM..'); i_algo = 1; 
    W = sum(A,3);
    W(isnan(W)) = 0;
    D = sqrt(sum(W,2)); D(D==0) = 1; dd = 1./(D*D');
    W = W.*dd;
    [vec,val] = eig(W,'nobalance');
    temp = sortrows([diag(val) vec'],-1);
    evecs = temp(1:K,2:end)';
    for i = 1:N
        if (norm(evecs(i,:))>0)
            evecs(i,:) = evecs(i,:)./norm(evecs(i,:));
        end
    end
    idx = kmeans(evecs,K,'emptyaction','singleton','replicates',10);
    E_exp(i_algo,i_exp,i_trial) = computeCE(idx,trueLabel);
    clear W D dd vec val temp evecs idx
 
    
    %% NH-Cut
    display('NH-Cut..'); i_algo = 2; 
    H = zeros(N,edgcount);
    for i = 1:edgcount
        H(edg(i,:),i) = 1;
    end
    W = (H.*repmat(edgwt',N,1))*(H'./3);
    W(isnan(W)) = 0;
    D = sqrt(sum(W,2)); D(D==0) = 1; dd = 1./(D*D');
    W = eye(N) - W.*dd;
    [vec,val] = eig(W,'nobalance');
    temp = sortrows([diag(val) vec'],1);
    evecs = temp(1:K,2:end)';
    for i = 1:N
        if (norm(evecs(i,:))>0)
            evecs(i,:) = evecs(i,:)./norm(evecs(i,:));
        end
    end
    idx = kmeans(evecs,K,'emptyaction','singleton','replicates',10);
    E_exp(i_algo,i_exp,i_trial) = computeCE(idx,trueLabel);
    clear H W D dd vec val temp evecs idx

    
    %% HOSVD
    display('HOSVD..'); i_algo = 3; 
    W = reshape(A,N,N^2);
    W(isnan(W)) = 0;
    W = W*W';
    D = sqrt(sum(W,2)); D(D==0) = 1; dd = 1./(D*D');
    W = W.*dd;
    [vec,val] = eig(W,'nobalance');
    temp = sortrows([diag(val) vec'],-1);
    evecs = temp(1:K,2:end)';
    for i = 1:N
        if (norm(evecs(i,:))>0)
            evecs(i,:) = evecs(i,:)./norm(evecs(i,:));
        end
    end
    idx = kmeans(evecs,K,'emptyaction','singleton','replicates',10);
    E_exp(i_algo,i_exp,i_trial) = computeCE(idx,trueLabel);
    clear W D vec val temp evecs idx
    
    %% SNTF
    display('SNTF..'); i_algo = 4; 
    W = reshape(A,N,N^2);
    W(isnan(W)) = 0;
    % normalization
    for t = 1:50
        a = diag(1./(sum(W,2).^(1/3)));
        a(isnan(a)) = 0;
        W = a*W;
        W = reshape(W',N,N^2);
        W = a*W;
        W = reshape(W',N,N^2);
        W = a*W;
        W = reshape(W',N,N^2);
    end
    
    G = rand(N,K);
    for j = 1:N
        G(j,:) = G(j,:)/sum(G(j,:));
    end
    
    for t = 1:50
        B1 = zeros(N,K);
        B = zeros(N,K);
        GG = G'*G;
        for r = 1:K
            B1(:,r) = W*kron(G(:,r),G(:,r));
            for s = 1:N
                B(s,r) = ((GG(r,:) - G(s,r)*G(s,:)).^2)*(G(s,:))';
            end
        end
        G = G.*B1./B;
        for j = 1:N
            G(j,:) = G(j,:)/sum(G(j,:));
        end
    end
    
    [~,idx] = max(G,[],2);
    E_exp(i_algo,i_exp,i_trial) = computeCE(idx,trueLabel);
    clear  W a G A B GG idx

    
    %% HGT
    display('HGT..'); i_algo = 5; 
    HGT = struct('nVertices',N,'nEdges',edgcount,'edges',edg,'w',edgwt);
    [idx,~] = hgtClustering(HGT,K,50);
%     E_exp(i_algo,i_exp,i_trial) = computeCE(idx,trueLabel);
    E_exp(i_algo,i_exp,i_trial) = (computeCE(idx(idx~=0),trueLabel(idx~=0))*sum(idx~=0) + sum(idx==0))/N;
    clear HGT idx
    

    %% h-METIS
    display('h-METIS..'); i_algo = 6; 
    cd hmetis-1.5-osx-i686
    str = 'hypergraph_planted.hgr';
    fid = fopen(str,'w');
    fprintf(fid,'%d %d\n',size(edg,1),N);
    fprintf(fid,'%d %d %d\n',edg');
    fclose(fid);
%     dlmwrite(str,[size(edg,1) N],'delimiter',' ')
%     dlmwrite(str,edg,'-append','delimiter',' ')
    system(['./shmetis ',str,' 3 1']);
    idx = dlmread([str,'.part.3']); idx = idx+1;
    cd ..
    E_exp(i_algo,i_exp,i_trial) = computeCE(idx,trueLabel);
    clear edg1 idx

    %%
    disp('     TTM      NH-Cut     HOSVD     SNTF      HGT     hMETIS')
    disp(mean(E_exp(:,:,1:i_trial),3).')
    toc(tStart)
end
end

disp('     TTM      NH-Cut     HOSVD     SNTF      HGT    hMETIS')
mE_exp = mean(E_exp,3).'

algo = {'TTM', 'NH-Cut', 'HOSVD', 'SNTF', 'HGT', 'hMETIS'};
figure;
set(gcf, 'Position', [100 100 1000 200]);
for i_algo=1:n_algo
    subplot(1,n_algo,i_algo);
    h = imagesc(flipud(reshape(mE_exp(:,i_algo),4,4)'),[0 1]); colormap(flipud(gray));
    pos = get(gca, 'Position');
    set(gca, 'Position', [pos(1) pos(2)+.1 pos(3) pos(4)-.2] );
    set(gca, 'xTick', 1:4);
    set(gca, 'xTickLabel', 10:10:40);
    set(gca, 'yTick', 1:4);
    set(gca, 'yTickLabel', 0.2:-0.05:0.05);
    title(algo{i_algo},'FontSize',14);
    if i_algo == 1
        ylabel('Probability gap (p-q)','FontSize',16);
    end
    if i_algo == 3
        xlabel('Number of points for each cluster (n/k)','FontSize',16);
    end
end
colorbar;