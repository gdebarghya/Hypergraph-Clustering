1. GENERAL INFORMATION

This toolbox contains the Matlab implementation of the 
Sparse Grassmann Clustering (SGC) approach proposed in:
 
Suraj Jain and Venu Madhav Govindu (2013). Efficient Higher-Order Clustering on the Grassmann Manifold. In 2013 IEEE International Conference on Computer Vision.

This is not the implementation used by the authors of the above paper. The implementation is based on that of the Spectral Curvature Clustering (SCC) algorithm by Guangliang Chen and Gilad Lerman (IJCV 2009). The code was written solely for research purpose, in particular for a comparative study of SGC, SCC and TeTrIS (arxiv:1602.06516).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

2. USAGE AND HELP

The main file is sgc.m, which performs clustering of affine subspaces. The simplest way to use them is as follows:

[sampleLabels, averageL2Error] = sgc(X,d,K);

where,

X: N-by-D data matrix 
d: intrinsic dimension
K: number of clusters

sampleLabels: cluster labels of the data points
averageL2Error: averaged L2 error of the detected clusters.

For more detailed description, as well as choices of optional parameters, 
please type in the Matlab command window

help sgc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COPYRIGHT NOTICE:

We do not hold any copyright for SGC algorithm, and we do not guarantee that this implementation resembles the original codes for SGC algorithm, which is not available.

The implementation is entirely based on the available codes for SCC algorithm, and have been used for research purpose only. The authors of SCC do not specify any restrictions regarding usage of these codes.

On our part, we do not impose any restriction the use or modification of this material for research purpose. For use of this material for commercial purpose, one must check the relevant copyright or policy imposed by the Indian Institute of Science, where this research was conducted.
