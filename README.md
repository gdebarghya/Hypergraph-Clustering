# Hypergraph-Clustering
MATLAB codes for tensor based methods for hypergraph partitioning and subspace clustering 

The repostory contains all implementation associated with the paper [1]. 
This also includes implementations of methods proposed in [2,3,4].

1. D. Ghoshdastidar and A. Dukkipati. Uniform hypergraph partitioning: Provable tensor methods and sampling techniques. 
*Journal of Machine Learning Research* 18(50), pp. 1−41, 2017. [paper](http://jmlr.org/papers/volume18/16-100/16-100.pdf)
1. D. Ghoshdastidar and A. Dukkipati. Consistency of spectral hypergraph partitioning under planted partition model. *The Annals of Statistics*, 45(1), pp. 289-315, 2017. [preprint](https://arxiv.org/pdf/1505.01582.pdf)
1. D. Ghoshdastidar and A. Dukkipati. A Provable Generalized Tensor Spectral Method for Uniform Hypergraph Partitioning. 
In *Proceedings of the 32nd International Conference on Machine Learning (ICML)*, PMLR 37:400-409, 2015.[paper](http://proceedings.mlr.press/v37/ghoshdastidar15.pdf)
1. D. Ghoshdastidar and A. Dukkipati. Consistency of spectral partitioning of uniform hypergraphs under planted partition model. In *Advances in Neural Processing Systems (NIPS)*, 2014. [paper](http://papers.nips.cc/paper/5314-consistency-of-spectral-partitioning-of-uniform-hypergraphs-under-planted-partition-model.pdf)

Please cite [1] if you use these codes/results in your work. If the methods in [2,4] are used, please cite accordingly.

## Copyright
Copyright(c) 2017 [Debarghya Ghoshdastidar](https://gdebarghya.github.io)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Contents of the repository
The repository contains Matlab implementation of all methods and experiments, and data files used by the codes. 
The codes grouped based on the problem considered:

#### Subspace clustering
The main codes for the generating the results in the paper are:
* `run_sc_synthetic.m`: Comparison of methods on synthetic data
* `run_sc_motionseg.m`: Comparison of methods on Hopkins 155 motion segmentation benchmark
* `tetris/tetris.m`: main code for the proposed Tetris subspace clustering algorithm

The codes are based on implementations by Dohyung Park for the paper [Greedy Subspace Clustering](arxiv.org/pdf/1410.8864.pdf). 
We have made few modifications in the experimental setup on synthetic data. 
We have also included the codes for the other methods.

#### Hypergraph partitioning
The main codes for the generating the results in the paper are:
* `run_consistency.m`: Comparison of TTM, HOSVD and NH-Cut under planted model
* `run_hypergraph_planted.m`: Comparison of different methods for 3-uniform planted model
* `run_hypergraph_subspace.m`: Comparison of different methods for line clustering problem
* `planted_hypergraph.m`: Contains codes for the proposed algorithms TTM [1,3], HOSVD [4] and NH-Cut [2] 

We have adapted portions of the experimental setup from the implementations by Dohyung Park for the paper [Greedy Subspace Clustering](arxiv.org/pdf/1410.8864.pdf).
We have also included the codes for the other hypergraph partitioning methods, such as hMETIS (Karypis & Kumar, VLSI 2000), HGT (Rota Bulo & Pelillo, TPAMI 2013) and our implementations of SNTF (Shashua, Zass & Hasan, ECCV 2006).

