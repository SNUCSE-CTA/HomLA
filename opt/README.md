# Homomorphic Computation of Local Alignment

A tool that computes an optimal local alignment between two encrypted genomic sequences with affine gap penalty.

## Requirement

OpenMP<br/>
[TFHE](https://tfhe.github.io/)<br/>

## Compile and run
```
make
```
```
./opt_tfhe inputFileX inputFileY [M S Go Ge]
```

## Sample run

```
./opt_tfhe ../x.fasta ../y.fasta
```
or
```
./opt_tfhe ../x.fasta ../y.fasta 10 -4 -5 -2
```
## Reference
[M. Bataa, S. Song, K. Park, M. Kim, J.H. Cheon, and S. Kim. Homomorphic Computation of Local Alignment, Proceedings of the IEEE BIBM, pp. 2167-2174, 2020.](https://ieeexplore.ieee.org/document/9313199)
