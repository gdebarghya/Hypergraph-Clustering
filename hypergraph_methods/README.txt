This folder contains codes for comparing several hypergraph partitioning algorithms.
These codes have been used for the results given in:

Debarghya Ghoshdastidar and Ambedkar Dukkipati (2016).
Uniform Hypergraph Partitioning: Provable Tensor Methods and Sampling Techniques. 
arXiv:1602.06516.

The main codes for the generating the results in the paper are
run_consistency.m: Comparison of TTM, HOSVD and NH-Cut under planted model
run_hypergraph_planted.m : Comparison of different methods for 3-uniform planted model
run_hypergraph_subspace.m : Comparison of different methods for line clustering problem

We have adapted portions of the experimental setup from the implementations by Dohyung Park for the NIPS-2014 paper:
Greedy Subspace Clustering (arxiv:1410.8864)

We have also included the codes for the other hypergraph partitioning methods, such as hMETIS (Karypis & Kumar, VLSI 2000), HGT (Rota Bulo & Pelillo, TPAMI 2013) and our implementations of SNTF (Shashua, Zass & Hasan, ECCV 2006).


