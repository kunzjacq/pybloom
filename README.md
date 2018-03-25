## Script to compute the false positive probability of Bloom filter variants

Based on ideas from the publication

'[Cache-, Hash- and Space-Efficient Bloom Filters](http://algo2.iti.kit.edu/documents/cacheefficientbloomfilters-jea.pdf)', F. Putze, P. Sanders and J. Singler

The explanations below assume a working knowledge of Bloom filters.

Structures considered are sets of small Bloom filters selected by a hash function.
an element is tested to belong to the structure by testing it against the filter
associated to it. The small size of the filter makes memory accesses faster. For
instance, on may want to use a 512-bit filter to make it fit in a cache line.

A filter may be implemented as several independent, smaller filters.
A 512-bit filter may for instance be replaced by 8 64-bit filters.
With smaller filters, the bit positions corresponding to elements to insert or
test can be pre-computed: the resulting bit masks can be stored in a mask
set of reasonable size.

When using precomputed masks, the use of several small filters, or
'cascading', enables to simulate a larger mask set than what would be possible
without cascading. e.g. if one is limited to a mask set of size 2^8, a
cascading of 2 is almost as good as a cascading of 1 with words of double
size, double hamming weight, and mask set of size 2^16.
However cascading does *not* enable to beat larger word size when there is no
limitation on mask set size (i.e. when masks are generated pseudo-randomly from
the inserted values).

### computing the false positive rate of a fixed structure

A structure with small bloom filters of 64-bits, each containing on
average 4 values, with masks of hamming weight 6, and a set of masks of
size 2^16 (thus requiring 16 hashed bits from the input value for mask
selection), one gets a collision probability of 0.37%, as is computed by
```local_bf.false_positive_proba(bloom_size=64, mask_weight=6, avg_loading=4, log2_mask_set_size=16, cascading_factor=1)```

The mask table has 2^16 values of 64 bits: it uses 2^19 bytes.

false_positive_proba returns 4 probabilities, the most interesting one being the
last one (compounded probability with cascading, and finite mask set). See
function documentation for more details.

there is an implicit parameter F that limits the maximum number of values in
the filter to its average value avg_loading plus F times its standard deviation
(which is sqrt(avg_loading)). Values larger than this threshold are not taken
into account in the collision probability computation. F defaults to 10 which
should always be equivalent to infinity.

A cascade 4 of small structures of size 16, each with 2^8 masks of hamming weight 3
has false positive rate is 1.05%, as computed by
```local_bf.false_positive_proba(bloom_size=16, mask_weight=3, avg_loading=4, log2_mask_set_size=8, cascading_factor=4)```

This structure requires 4\*8 = 32 hashed bits per value for mask selection;
the mask table is composed of 2^8 16-bit words, i.e. 2^10 bytes. The same table
can be reused for the different cascaded filters, as long as the randomness
used is different.

An optimiser is provided to find the structure with lowest false positive probability
under constraints.

### Using the optimizer

Constraint on total storage size per value, total access size to check one
element, mask set size and cascading:
```local_bf.optimiser(log2_storage_size=4,max_log2_access_size=6,max_log2_cascading=3, max_log2_mask_set_size=8)```

one can additionally impose a constraint on the word size:
```local_bf.optimiser(log2_storage_size=6, max_log2_access_size=8, max_log2_filter_size=5, max_log2_mask_set_size=8,max_log2_cascading=4)```

or one can impose a constraint on the mask set storage size (log, in bits):
```local_bf.optimiser(log2_storage_size=6, max_log2_access_size=8, max_log2_filter_size=6, max_log2_mask_set_size=8, max_log2_cascading=2,max_log2_mask_storage_bit_size=14)```

to obtain very low false positive rates, it is often useful to enable a filter
set size larger than the storage size allocated per element. Compare:
```local_bf.optimiser(log2_storage_size=7,max_log2_access_size=6,max_log2_cascading=3, max_log2_mask_set_size=8,max_log2_mask_storage_bit_size=13)```

which yields a 512-bit structure with f.p. rate of 1.275e-9 with finite masks,
and average loading 4,

and
```local_bf.optimiser(log2_storage_size=7,max_log2_access_size=6,max_log2_cascading=3, max_log2_mask_set_size=8,max_log2_mask_storage_bit_size=13, max_log2_filterset_size=7)```

where a constraint on the total filter set size was added: it cannot be larger
than 128 bits. As a result, the f.p. rate climbs to 8.7e-8.

see optimiser parameter definition for general usage.

The code uses of the following result:
the limit of the log-probability that there are u values in a filter after
n elements are inserted at random into m filters, n -> infinity, n/m = a is constant,
is
-log(u!) - a + u * log(a)

this is obtained by a Taylor development of binomial(n,u) (1/m)^u (1-1/m)^n-u
when m tends to infinity with m = a*n.