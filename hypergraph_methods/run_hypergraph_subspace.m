% Code by Debarghya Ghoshdastidar for the paper arXiv:1606.06516

clear all; close all; clc;
n_settings = 16;

p = 3;                       % Ambient dimension
d = 1;                       % subspace dimension
K_exp = 3;               % # subspaces
n_exp = repmat(10:10:40,1,4);
o_exp = reshape(repmat(0:0.05:0.15,4,1),1,n_settings);
             % # sample points for each subspace

n_algo = 6;          
n_trial = 50;               

E_exp = zeros(n_algo,n_settings,n_trial);     % Clustering error
tStart = tic;
for i_trial = 1:n_trial          % # trials
for i_exp = 1:n_settings
    
    disp(' ')
    disp(['Trial ' int2str(i_trial) '/' int2str(n_trial) ' : Setting ' int2str(i_exp) '/' int2str(n_settings)])
    
    o = o_exp(i_exp);            
    K = K_exp;
    n = n_exp(i_exp)*ones(1,K);  
    N = sum(n);                  % Total # of samples

    %% True subspace generation   
    Vss = normc(rand(p,K) - 0.5);       % True lines
    Vn = o*normc(rand(p,N) - 0.5);      % Noise
    
    %% Data point generation
    trueLabel = zeros(1,N);           % True labels for sample points
    count = 0;
    for i = 1:K
        trueLabel((count+1):(count+n(i))) = i;
        count = count+n(i);
    end
        
    Y  = zeros(p,N);           % Sample points
    for i=1:N
        Y(:,i) = 2*(rand(1)-0.5)*Vss(:,trueLabel(i)) + Vn(:,i);
    end
    disp('Data generated.')
    
    %% Adjacency tensor and incidence matrix
    A = inf*ones(N,N,N);
    edg = zeros(nchoosek(N,3),3);
    edgwt = inf*ones(nchoosek(N,3),1);
    edgcount = 0;
    
    Af = computing_polar_curvatures(Y',combnk(1:N,2).',1);
    combcount = 0;
    for i1 = 1:N
        for i2 = (i1+1):N
            combcount = combcount+1;
            for i3 = (i2+1):N
                edgcount = edgcount + 1;
                edg(edgcount,:) = [i1 i2 i3];
                edgwt(edgcount,1) = Af(i3,combcount);
                
                A(i1,i2,i3) = edgwt(edgcount,1); A(i1,i3,i2) = edgwt(edgcount,1);
                A(i2,i1,i3) = edgwt(edgcount,1); A(i2,i3,i1) = edgwt(edgcount,1);
                A(i3,i2,i1) = edgwt(edgcount,1); A(i3,i1,i2) = edgwt(edgcount,1);                
            end
        end
    end
    [~,ind] = sort(edgwt,'ascend');
    optSigma = edgwt(ind(round(length(edgwt)/K^2)));  
    % this sigma is suggested in (Chen & Lerman, IJCV-2009)
    A = exp(-(A.^2)./(2*optSigma^2));%double(A<=optSigma);%
    edgwt = exp(-(edgwt.^2)./(2*optSigma^2));%double(Hwt<=optSigma);%
    
    
    %% TTM
    display('TTM..'); i_algo = 1; 
    W = sum(A,3);
    W(isnan(W)) = 0;
    D = sqrt(sum(W,2)); D(D==0) = 1; dd = 1./(D*D');
    W = W.*dd;
%     [evecs,~] = eigs(W,K,'lm');
    [vec,val] = eig(W);
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
%     [evecs,~] = eigs(W,K,'sm');
    [vec,val] = eig(W);
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
%     [evecs,~,~] = svds(W,K);
    [vec,val,~] = eig(W);
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
    for t = 1:10
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
    
    for t = 1:20
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
    [idx,~] = hgtClustering(HGT,K,10);
    E_exp(i_algo,i_exp,i_trial) = (computeCE(idx(idx~=0),trueLabel(idx~=0))*sum(idx~=0) + sum(idx==0))/N;
    clear HGT idx
    
    
    %% h-METIS
    display('h-METIS..'); i_algo = 6; 
    str = 'hypergraph_subspace.hgr';
    edg1 = [round(100*(edgwt./max(edgwt))) edg];
    edg1 = edg1(edg1(:,1)>0,:);  
    edg1 = sortrows(edg1,-1);
    cd hmetis-1.5-osx-i686
    fid = fopen(str,'w');
    fprintf(fid,'%d %d %d\n',size(edg1,1),N,1);
    fprintf(fid,'%d %d %d %d\n',edg1');
    fclose(fid);
%     dlmwrite(str,[size(edg1,1) N 1],'delimiter',' ')
%     dlmwrite(str,edg1,'-append','delimiter',' ')
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

disp('     TTM      NH-Cut     HOSVD     SNTF      HGT     hMETIS')
mE_exp = mean(E_exp,3).'

algo = {'TTM', 'NH-Cut', 'HOSVD', 'SNTF', 'HGT', 'hMETIS'};
figure;
set(gcf, 'Position', [100 100 1000 200]);
for i_algo=1:n_algo
    subplot(1,n_algo,i_algo);
    pos = get(gca, 'Position');
    set(gca, 'Position', [pos(1) pos(2)+.1 pos(3) pos(4)-.2] );
    h = imagesc(flipud(reshape(mE_exp(:,i_algo),4,4)'),[0 1]); colormap(flipud(gray));
    set(gca, 'xTick', 1:4);
    set(gca, 'xTickLabel', 10:10:40);
    set(gca, 'yTick', 1:4);
    set(gca, 'yTickLabel', 0.15:-0.05:0);
    title(algo{i_algo},'FontSize',14);
    if i_algo == 1
        ylabel('Noise level (\sigma_a)','FontSize',16);
    end
    if i_algo == 3
        xlabel('Number of points for each subspace (n/k)','FontSize',16);
    end
end
colorbar;
