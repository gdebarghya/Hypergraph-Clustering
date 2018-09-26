function [sampleLabels,averageL2Error] = ttm(X,d,K,OPTIONS)

%   Tensor Trace maximization
%
%   [sampleLabels, averageL2Error] = ttm(X,d,K) 
%   partitions the points 
%   in the N-by-D data matrix X into K clusters, each representing a d-flat.
%   Rows of X correspond to points, columns correspond to variables.  
%   TTM returns an N-by-1 vector sampleLabels containing the cluster 
%   label of each point and the averaged L2 error of the K detected 
%   clusters. Those with zero labels are detected outliers.  
%
%   [ ... ] = ttm(..., OPTIONS) allows you to specify optional parameters 
%   to control the algorithm.  OPTIONS is a structure array consisting of
%   the following fields:
%
%   'n' - number of points used for computing a curvature
%         default = d+2, but can be larger than d+2
%
%   'c' - number of columns sampled from the matrix A, default = K*100 
%
%   'normalizeW' - 0/1, if we normalize the matrix W. default = 1(YES)
%
%   'normalizeU' - 0/1, if we normalize the matrix U. default = 1(YES)
%
%   'sigma' - the tuning parameter used in the affinity tensor
%             can be specified by user when findOptimalSigma = 0; 
%             otherwise the algorithm will infer its optimal value from data
%
%   'alpha' - number (if >=1) or percentage (if <1) of outliers in the data; 
%             default = 0
%
%   'initialLabels' - initial labelling of data, not a required parameter for scc.
%                     This option allows scc to improve clusters obtained by 
%                     other algorithms (e.g., k-flats, GPCA)

%   'seedType' - method of selecting initial seeds for kmeans clustering
%       'hard' - the deterministic procedure presented in the IJCV paper
%       'soft' - probabilistic procedure for selecting initial seeds
%
%
%   (c)2016 Debarghya Ghoshdastidar
%   Last updated on 23/02/2016.
%   Most of the supporting codes are directly borrowed for the implementation
%   of the SCC algorithm by G. Chen and G. Lerman (IJCV, 2009).
%   If you have any questions please email debarghya.g@csa.iisc.ernet.in
%
%   Most Relevant Publications:
%   1. Uniform Hypergraph Partitioning: Provable Tensor Methods and Sampling Techniques, 
%      D. Ghoshdastidar and A. Dukkipati, arXiv:1602.06516.

%% Set default values for OPTIONS
ABSOLUTE_MINIMUM = 1e-15;

if nargin < 4
    OPTIONS = struct();
end

if ~isfield(OPTIONS,'n') || OPTIONS.n < d+2 ... 
        || (OPTIONS.n > d+2 && d == 0)
    OPTIONS.n = d+2;
end

if ~isfield(OPTIONS,'c')
    OPTIONS.c = K*100;
else
    if mod(OPTIONS.c,K)>0
        OPTIONS.c = K*ceil(OPTIONS.c/K);
        warning('The value of the parameter c in OPTIONS has been modified to be an integer multiple of K!'); %#ok<WNTAG>
    end
end

if ~isfield(OPTIONS,'normalizeW')
    OPTIONS.normalizeW = 1;
end

if ~isfield(OPTIONS,'normalizeU')
    OPTIONS.normalizeU = 1;
end

if ~isfield(OPTIONS,'seedType')
    OPTIONS.seedType = 'hard';
end

if ~isfield(OPTIONS,'alpha')
    OPTIONS.alpha = 0;
end

N = size(X,1);

%% initialization

% auxiliary 
averageL2Error = Inf;
sampleLabels = zeros(N,1);

if ~isfield(OPTIONS, 'initialLabels')
    % consider given data as one cluster
    sampleLabels1 = ones(N,1);
else
    % only for improving clusters obtained by other algorithms
    sampleLabels1 = OPTIONS.initialLabels;
end

averageL2Error1 = computing_average_L2_error(X, d, sampleLabels1);

%% main body of the code
q = OPTIONS.n; % q = d+2 by default

sampleLabels = sampleLabels1;
averageL2Error = averageL2Error1;

%if isfield(OPTIONS,'sampledColumns')
%    sampledColumns = OPTIONS.sampledColumns;
%    OPTIONS = rmfield(OPTIONS,'sampledColumns');
%else
%   sampledColumns = sampling_columns(sampleLabels1,K,OPTIONS);
%end

sampledColumns = sampling_columns(sampleLabels1,K,OPTIONS);

polarCurv = computing_polar_curvatures(X,sampledColumns,d);
polarCurv_sorted = sort(reshape(polarCurv,1,N*OPTIONS.c));

q = max(q-1,1);


if isfield(OPTIONS,'sigma')
    sigma = OPTIONS.sigma;
else
    sigma = max(polarCurv_sorted(1,ceil((N-OPTIONS.n+1)*OPTIONS.c/K)),ABSOLUTE_MINIMUM);
    %sigma = max(polarCurv_sorted(1,ceil((N-OPTIONS.n+1)*OPTIONS.c/K^(OPTIONS.n-1))),ABSOLUTE_MINIMUM);
end

isigma = 1/(2*sigma);
A = exp(-polarCurv*isigma);
sampleLabels1 = ttm_processing_affinities(A,sampledColumns,K,OPTIONS);
averageL2Error1 = computing_average_L2_error(X, d*ones(K,1), sampleLabels1);


if averageL2Error1 < averageL2Error 
    averageL2Error = averageL2Error1;
    sampleLabels = sampleLabels1;
end
%figure; do_plot_data(X,sampleLabels); title('clusters obtained by SCC');