# Adaptive Randomized Pivoting

Code and data for experiments around adaptive randomized pivoting (ARP) and related column subset selection methods.

Paper:
- Adaptive randomized pivoting and volume sampling
- arXiv: https://arxiv.org/abs/2510.02513

## Repository Structure

- code/: Core MATLAB implementations and MEX-accelerated routines.
- tests/: Scripts used to reproduce timing and accuracy experiments.
- data/: Saved experiment outputs and processed datasets.
- figs/: Saved MATLAB figure files.
- utils/: Helper utilities (plotting, setup, sparse matrix helpers).

## Adaptive randomized pivoting

Adaptive randomized pivoting (ARP) is an algorithm for construction an interpolative low-rank approximation to a matrix $A$.
That is, it finds an interpolation matrix $W$ and an index set $S$ for which $A \approx W \cdot A(S,:)$.
The adaptive randomized pivoting was originally discovered by [Alice Cortinovis](https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://sites.google.com/view/alicecortinovis/home&ved=2ahUKEwjev5vAntCTAxXpOTQIHUUpMN8QFnoECCwQAQ&usg=AOvVaw0AdxojA7J5UwHURmdHpJBP) and [Daniel Kressner](https://people.epfl.ch/daniel.kressner?lang=en); see their original paper [here](https://arxiv.org/abs/2412.13992).
In my paper, I develop connections between ARP and volume sampling, leading to faster implementations of the ARP method.

## Quick Start (MATLAB)

To run ARP, execute the following MATLAB command:

```matlab
addpath("utils")
addpath("code")
arp_startup

A = randn(2000, 500);   % Example data matrix
k = 50;                 % Target subset size

[idx, W] = arp(A, k);
Ahat = W * A(idx, :);
relerr = norm(A - Ahat, "fro") / norm(A, "fro")
```

The function ARP outputs an interpolation matrix `W` and a set of rows indices `idx`, which define a low-rank approximation `Ahat = W * A(idx,:)`.

We give several versions of the ARP method; see the paper for discussion of different ARP method variants.

```matlab
[idx, W] = arp(A, k, "rejection");              % Default ARP variant
[idx, W] = arp(A, k, "ck");                     % CK baseline variant
[idx, W] = arp(A, k, "rejection", "sketchy"); % Sketchy interpolation
[idx, W] = arp(A, k, "rejection", "optimal"); % Optimal interpolation
```

## Reproducing Experiments

From the tests/ directory, run scripts such as:

- test_arp_dense.m
- test_arp_sparse.m
- test_arp_accuracy.m

These scripts generate .mat outputs in data/ and figures in figs/.

## Rebuilding MEX Files

To compile from source in MATLAB, run the following commands:

```matlab
cd code
mex rejection_helper.cpp
mex sparsestack.c
```

## Citation

If you use this repository, please cite:

```bibtex
@misc{epperly2025adaptiverandomizedpivotingvolume,
  title={Adaptive randomized pivoting and volume sampling},
  author={Ethan N. Epperly},
  year={2025},
  eprint={2510.02513},
  archivePrefix={arXiv},
  primaryClass={stat.ML},
  url={https://arxiv.org/abs/2510.02513}
}
```
