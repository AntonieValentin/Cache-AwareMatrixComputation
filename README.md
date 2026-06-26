# Cache-Aware Matrix Computation

A C project focused on implementing and optimizing a fixed matrix expression. The goal is to compute the same mathematical result in three different ways: a plain manual version, a BLAS-based version, and an optimized manual version.

The project is mainly about understanding how implementation choices affect performance on modern processors, especially for matrix multiplication and memory-heavy code.

## Computation

The implemented expression is:

```text
C = A^T * B
D = C * C^T
y = 0
for i = 0 to N - 1:
    y += D * row(C, i)
y = y + x
```

`A` and `B` are square `N x N` matrices of `double` values, while `x` and `y` are vectors of size `N`. Since `D = C * C^T`, the implementation also takes into account that `D` is symmetric.

## Implementations

The project contains three versions of the same computation.

### Naive Version

The `neopt` version is a straightforward manual implementation. It keeps the algorithm close to the mathematical formula and is useful as a clear baseline for correctness and performance comparisons.

### BLAS Version

The `blas` version uses BLAS routines for the main matrix and vector operations. This version shows how much performance can be gained by relying on highly optimized numerical libraries.

### Optimized Manual Version

The `opt_m` version improves the manual implementation without changing the mathematical complexity. The optimizations are focused on better memory access patterns, avoiding unnecessary work, reusing computed values, and reducing cache misses where possible.

## Project Structure

* `solver_neopt.c` - baseline manual implementation;
* `solver_blas.c` - implementation using BLAS routines;
* `solver_opt.c` - optimized manual implementation;
* `README.md` - project description and implementation notes.

## Performance Focus

The interesting part of the project is not only getting the correct result, but also comparing how each version behaves for larger matrix sizes.

The optimized version is designed with cache behavior in mind. Matrix code can become slow when it repeatedly accesses memory in an unfriendly order, so loop structure and data reuse matter a lot.

## Notes

This project is a practical exercise in performance-oriented programming. The same formula can be written in several correct ways, but the runtime can be very different depending on how well the implementation uses the CPU cache and available numerical routines.

The final result is a compact comparison between direct C code, library-based computation, and manual optimization.
