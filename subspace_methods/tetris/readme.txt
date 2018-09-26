1. GENERAL INFORMATION

This toolbox contains the Matlab implementation of the 
TeTrIS algorithm (Tensor Trace maximization with Iterative Sampling) proposed in:
 
Debarghya Ghoshdastidar and Ambedkar Dukkipati (2016). Uniform Hypergraph Partitioning: Provable Tensor Methods and Sampling Techniques. arXiv:1602.06516.


The implementation is based on that of the Spectral Curvature Clustering algorithm 
by Guangliang Chen and Gilad Lerman (IJCV 2009). 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

2. USAGE AND HELP

The main file is tetris.m, which performs clustering of affine subspaces. The simplest way to use them is as follows:

[sampleLabels, averageL2Error] = tetris(X,d,K);

where,

X: N-by-D data matrix 
d: intrinsic dimension
K: number of clusters

sampleLabels: cluster labels of the data points
averageL2Error: averaged L2 error of the detected clusters.

For more detailed description, as well as choices of optional parameters, 
please type in the Matlab command window

help tetris

If you have any questions please email Debarghya Ghoshdastidar at:
 debarghya.g@csa.iisc.ernet.in
(NOTE: email id may change based on affiliation. Please check present email id at Debarghya’s webpage)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COPYRIGHT NOTICE:

The implementation is entirely based on the available codes for SCC algorithm, and have been used for research purpose only. The authors of SCC do not specify any restrictions regarding usage of these codes.

On our part, we do not impose any restriction the use or modification of this material for research purpose. However, if this implementation is used, we request you to cite the relevant paper mention above. 

For use of this material for commercial purpose, one must check the relevant copyright or policy imposed by the Indian Institute of Science, where this research was conducted.
