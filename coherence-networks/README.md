# coherence-networks

Graph-theoretic summary of scalp EEG coherence networks, with the thresholding sensitivity checks used in the paper.

## What the method does

Magnitude-squared coherence between all electrode pairs gives a symmetric weighted matrix per visit and frequency band. The matrix is thresholded either proportionally (the strongest 15 to 25% of edges are kept) or reduced to its maximum spanning tree, which is threshold-free and insensitive to changes in mean coherence across visits. Global efficiency, characteristic path length, clustering coefficient and Louvain modularity are computed on the binarised graph with the Brain Connectivity Toolbox; small-world sigma is computed against degree-preserving null networks. Regional coherence means (left, right, cross-hemispheric, midline) are taken from the unthresholded matrix using the 10-20 label convention.

In the paper, global efficiency fell by about 10% at the post-implantation visit and recovered by week 24, in every subject, and the effect was broadband and distributed. The same pipeline is run on debiased weighted phase lag index matrices as a volume-conduction control by passing `'Negative', 'zero'`.

## Files

| File | Purpose |
|---|---|
| `coherence_network_metrics.m` | One coherence matrix in, global metrics and regional means out. |
| `example_usage.m` | Synthetic 19-channel network. |

## Requirements

MATLAB R2020b or later and the Brain Connectivity Toolbox (Rubinov and Sporns, 2010) on the path.

## Usage

```matlab
out = coherence_network_metrics(W, labels);                       % 25% proportional threshold
out = coherence_network_metrics(W, labels, 'Threshold', 'mst');   % maximum spanning tree
out = coherence_network_metrics(W, labels, 'SmallWorld', true);   % adds sigma
```

