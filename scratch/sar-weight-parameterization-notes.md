# SAR Weight Parameterization Notes

## SAR Precision

A simultaneous autoregressive (SAR) spatial random effect is often written as

```text
u = rho W u + e
```

or equivalently with precision

```text
Q = (I - rho W)' (I - rho W)
```

The interpretation of `rho` depends directly on how the spatial weights matrix
`W` is constructed.

## Raw Adjacency Weights

With raw adjacency weights, `W_ij = 1` when areas `i` and `j` are neighbors and
`W_ij = 0` otherwise. The SAR relationship is then

```text
u_i = rho * sum_j W_ij u_j + e_i
```

Here, `rho` multiplies the raw neighbor sum. This means `rho` is effectively an
edge-level feedback coefficient, and the total neighbor influence depends on how
many neighbors an area has.

For example, if `rho = 0.08`:

```text
area with 3 neighbors: total neighbor weight ~= 3 * 0.08 = 0.24
area with 8 neighbors: total neighbor weight ~= 8 * 0.08 = 0.64
```

Thus, the same `rho` implies stronger spatial feedback for highly connected
areas than for sparsely connected areas.

The valid range of `rho` also depends on the graph. For raw adjacency, the upper
bound is governed by the eigenvalues of `W`; roughly, positive `rho` must stay
below `1 / lambda_max(W)` to avoid singularity or instability. This bound varies
across spatial graphs.

## Row-Normalized Weights

With row-normalized weights, each non-island row of `W` sums to 1. The SAR
relationship becomes

```text
u_i = rho * average_neighbor_value_i + e_i
```

or, more generally, dependence on a weighted neighborhood average.

In this parameterization, `rho` has a more stable interpretation: it is the
strength of dependence on the neighborhood average rather than on the raw
neighbor sum. A value such as `rho = 0.6` has broadly comparable meaning for
areas with 3 neighbors and areas with 8 neighbors.

Row-normalization also gives a simpler and more portable parameter scale. For
typical row-normalized adjacency matrices, constraining

```text
-1 < rho < 1
```

is a natural and interpretable default.

## Practical Recommendation

For most SAR areal models, row-normalized weights are preferable because:

- `rho` has a clearer interpretation as dependence on the neighborhood average.
- The meaning of `rho` is more comparable across areas with different numbers of
  neighbors.
- The same model specification is more portable across different graphs.
- The usual `(-1, 1)` parameter bound is more defensible.

Raw adjacency weights can still be useful for exact software comparisons or for
models where degree-dependent neighbor accumulation is intentionally desired.
However, with irregular areal graphs, raw weights can make spatial dependence
partly reflect graph degree rather than purely spatial association.
