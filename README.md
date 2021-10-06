# Computation of False Positive Probability for Bloom filter variants

This project is based on ideas from the publication '[Cache-, Hash- and Space-Efficient Bloom Filters](http://algo2.iti.kit.edu/documents/cacheefficientbloomfilters-jea.pdf)', by F. Putze, P. Sanders and J. Singler.

Bloom filters enable to check whether elements belong to a set, with no false negatives and a computable false positive (f.p.) probability; see for instance https://en.wikipedia.org/wiki/Bloom_filter. For a given f.p. probability, Bloom filters storage requirements as a function of the set size are modest: the storage per elements in bits, *a*, and the f.p. probability, *p*, are related by *-ln(p) ≈ a × ln(2)<sup>2</sup> ≈ 0.5 × a*. However, testing an element in a Bloom filter requires several memory reads at independent, pseudorandom positions in the structure. This makes Bloom filter implementations rather cache-unfriendly and slow. Alternate structures where memory accesses are grouped together, at the expense of the memory requirements or f.p. probability, are therefore desirable for the uses where speed matters most.

## Better memory locality with a hash table of Bloom filters

The paper cited above considers structures that are collections of small Bloom filters where elements to be stored are associed to a particular filter by a hash function. An element is then tested to belong to a structure representing a set by testing it against the filter associated to it. For instance, one may want to use a collection of 512-bit Bloom filters to make each filter fit in a cache line.

## Cascading smaller filters

As a second optimization step, each Bloom filter may be approximated by several independent, smaller filters.
The insertion (resp. test) of an element is done by inserting (resp. testing) the element in each smaller filter, using (pseudo-)independently computed bit positions.
A 512-bit filter may for instance be approximated by 8 64-bit filters. Formally, the set represented by the collection of small filters is the intersection of the sets corresponding to each individual filter.

The f.p. probability of such set of smaller filters is worse than with a monolithic Bloom filter of the same size. Sufficently small filters are nevertheless amenable to a simplification of their implementation:
pre-computed *bit masks* with some fixed hamming weight can be chosen once and for all and stored in a table; then, when an element is inserted or tested, a mask index in the table is chosen pseudorandomly. As a result, all bits to set or test corresponding to the element processed are chosen in one operation. When a filter fits into a machine word, operations on a filter reduce to simple logical operations. When several small filters are used, the same mask table can be used for all filters, with mask indexes computed (pseudo-) independently for different filters.

To sum up, a set of small filters trades some f.p. probability in exchange for a simpler and faster implementation, with more predictible memory accesses.

Using several small filters as an approximation of a larger one is called 'cascading' filters in the sequel of this document.

## Computing the false positive rate of a fixed structure: examples

All code examples assume the module 'bloom' is imported at toplevel with:

```from bloom import *```

A structure using 64-bit Bloom filters, each containing on average 4 values, with masks of hamming weight 6, and a set of masks of size 2<sup>16</sup> (thus requiring 16 hashed bits from the input value for mask selection), has a f.p. probability of 0.37%, as computed by

```false_positive_proba(bloom_size=64, mask_weight=6, avg_loading=4, log2_mask_set_size=16, cascading_factor=1)```

Since the mask table contains 2<sup>16</sup> 64-bit values, it uses 2<sup>19</sup> bytes of storage.

`false_positive_proba` returns 4 probabilities, the most interesting one being the last one (compounded probability with cascading, and finite mask set). See function documentation for more details.

A cascade of 4 16-bit filters, each with 2<sup>8</sup> masks of hamming weight 3 has false positive rate equal to 1.05%, as computed by

```false_positive_proba(bloom_size=16, mask_weight=3, avg_loading=4, log2_mask_set_size=8, cascading_factor=4)```

This structure requires 4 × 8 = 32 hashed bits per value for mask selection; the mask table is composed of 2<sup>8</sup> 16-bit words, i.e. 2<sup>10</sup> bytes. The same table can be reused for the different cascaded filters, as long as the randomness used for mask selection in each filter is different.

`false_positive_proba` has an optional parameter *F* that defines the maximum number of values in the filter that is considered during computations. This maximum number *U* is defined as *a + F * s*, where *s = a<sup>1/2</sup>* is the loading standard deviation. Values larger than *U* are not taken into account in the f.p. probability computation. *F* defaults to 10 which should always be equivalent to infinity for practical purposes.

An optimiser is provided to find the structure with lowest false positive probability under constraints.

## Using the optimizer

Given constraints on a filter, or on a filter set composed of several smaller cascaded Bloom filters, the optimizer outputs a set of settings satisfying the constraints yielding the lowest false positive rate.

For instance, it can be given constraints on

 - total storage size per value
 - total access size to test one element
 - maximum mask set size
 - maximum cascading

In the following example, total storage size per value is limited to 2<sup>4</sup>=16 bits, access size to test one element is at most 2<sup>6</sup>=64 bits, maximum cascading is 2<sup>3</sup>=8 levels, maximum number of masks is 2<sup>8</sup>=256:

```optimiser(log2_storage_size=4,max_log2_access_size=6,max_log2_cascading=3, max_log2_mask_set_size=8)```

The result is

```
** search parameters **
storage size per element: 2^4 = 16
Maximum access size to check one element: 2^6 = 64 bits = 8 bytes
Maximum mask set size: : 2^8
Minimum Hamming weight of masks: 1
Maximum filter size: same as access size
Maximum number of iterations: 2^3 = 8

Maximum number of elements per filter considered: 88
480 settings analyzed; 446 eligible settings found

**** best result: ****

collision probability
 random masks: 1.250e-03
 finite mask set: 2.522e-03
 [in one filter, with finite mask set: 1.958e-01]

settings:
 word size: 64
 cascading factor: 4
 loading factor (average number of values per filter): 2^4 = 16
 hamming weight of masks: 2
 mask set size: 2^8.000 <= binomial(64, 2) = 2^10.977
 storage size = word size * cascading factor / loading = 64 * 4 / 16 = 16

Misc info:
 Mask set log storage size, in bits: 14.00; in bytes: 11.00
 random bits required per element insertion: 32.00, i.e. 4.00 bytes

Bloom Filter comparison:
 optimal number of hash functions for a unique BF with the same storage size per value: 11
 collision probability of this optimal BF: 4.587e-04
 practical vs optimal collision probability ratio: 5.5
```

Therefore one sees that a false positive probability of 2.522e-03 can be obtained with these constraints, and the settings to obtain this probability are given.

A constraint on the smaller filter size can be added:

```optimiser(log2_storage_size=6, max_log2_access_size=8, max_log2_filter_size=5, max_log2_mask_set_size=8,max_log2_cascading=4)```

or one can impose a constraint on the mask set storage size (log, in bits):

```optimiser(log2_storage_size=6, max_log2_access_size=8, max_log2_filter_size=6, max_log2_mask_set_size=8, max_log2_cascading=2,max_log2_mask_storage_bit_size=14)```

See optimiser parameter definition for general usage.

### Considerations about the way to obtain efficient constructs

To obtain low false positive rates, it is often useful to enable a filter set size larger than the storage size allocated per element. (remember that "filter set" refers to a collection of smaller filters assembled through cascading to form a larger structure analogous to a plain Bloom filter).

Compare:

```optimiser(log2_storage_size=7,max_log2_access_size=6,max_log2_cascading=3, max_log2_mask_set_size=8,max_log2_mask_storage_bit_size=13)```

which yields a 512-bit filter set with f.p. rate of 1.275e-9 with finite masks, and average loading (=average number of elements stored in a filter set) of 4 (which is equal to the filter set size divided by the storage size per element), and

```optimiser(log2_storage_size=7,max_log2_access_size=6,max_log2_cascading=3, max_log2_mask_set_size=8,max_log2_mask_storage_bit_size=13,max_log2_filterset_size=7)```

where a constraint on the filter set size was added: it cannot be larger than 128 bits. As a result, the false positive rate climbs to 8.7e-8: **this is a 68 times higher f.p. probability, for the same storage size per element.** The average loading factor is now 1.

The degraded behavior of filter sets is due to the fact that they do not enable a good averaging of the number of values inserted into them. As a result, there is a significant probability that a filter set is overcrowded (because of the random behavior of the hash function that distributes elements between filter sets). Since the same number of elements is inserted in every individal filter of the filter set, they are all simultaneously overcrowded, which is bad for the overall false positive probability.

## False positive probability computation

To compute the false positive probability, `false_positive_proba` only needs to model the behavior of one Bloom filter or one set of smaller cascaded Bloom filters under a variable load. Indeed, if *n* elements are inserted at random into m filters, the average number of elements in a filter set is *a = n/m* and the probability to have *u* elements in the filter is

*p = binomial(n,u) (1/m)<sup>u</sup> (1-1/m)<sup>(n-u)</sup>.*

For large *n*, *log(p)* tends to *-log(u!) - a + u * log(a)*.

More precisely, this is the limit of the log-probability that there are *u* values in a filter after *n* elements are inserted at random into *m* filters when *n* tends to infinity and *n/m = a* is constant. This is obtained by a Taylor development of *log(p)* starting with the exact formula above, when *n* tends to infinity and *m = a × n*.
