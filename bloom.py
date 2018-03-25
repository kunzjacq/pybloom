#!/usr/bin/python3

import math

def optimiser(log2_storage_size, max_log2_access_size, max_log2_cascading,
    max_log2_filterset_size = 0,
    max_log2_mask_set_size = 0,
    max_log2_mask_storage_bit_size=0, min_hamming_weight=1, max_hamming_weight=0, \
    max_log2_filter_size=0, F=10):
  '''
  variables with name log2 are logs in base 2.
  a constraint set to 0 is ignored.
   only log2_storage_size, max_log2_access_size, max_log2_cascading must be set.
  log2_storage_size: log of number of bits allocated per element, i.e.
   word_size times number of words allocated per value times cascading factor.
  max_log2_mask_set_size: log of the maximum max set size.
  max_log2_cascading: the log of the maximum number of times the whole structure
   is repeated.
  max_log2_mask_storage_bit_size: log of the maximum storage for mask set, in bits.
  max hamming weight: maximum hamming weight considered for masks.
  max_log2_access_size: log of the maximum number of bits accessed to check a
   particular value. Equal to filter size * cascading factor
  max_log2_filter_size: the log of the maximum filter size.
  F : number of standard deviations of the average number of elements present in a
  filter that are considered in computions. We ensure that the maximum number
  of masks is above F * the average number of elements. Furthermore, the probability
  that the number of elements in a filter is above that threshold is neglected.
  '''
  print_search_parameters(log2_storage_size, max_log2_access_size, max_log2_mask_set_size, \
    max_log2_mask_storage_bit_size, min_hamming_weight, max_hamming_weight, \
    max_log2_filter_size, max_log2_cascading)
  def max_log2_set_size(log2maxq, log2_word_size):
    res = log2maxq
    if max_log2_mask_set_size > 0:
      res = min(res, max_log2_mask_set_size)
    if max_log2_mask_storage_bit_size > 0:
      res = min(res, max_log2_mask_storage_bit_size - log2_word_size)
    return res

  results = []
  settings_analyzed = 0
  min_log2_filter_size = 3 #minimum filter size: 1 byte
  # one filter must fit in one memory access
  if max_log2_filter_size == 0:
    max_log2_filter_size = max_log2_access_size
  else:
  	max_log2_filter_size = min(max_log2_filter_size, max_log2_access_size)
  if max_hamming_weight == 0:
    max_hamming_weight = 2**max_log2_access_size
  max_log2_loading_global = max_log2_filter_size + max_log2_cascading - log2_storage_size
  max_u_global = int(2**max_log2_loading_global + F * math.sqrt(2**max_log2_loading_global))
  print('Maximum number of elements per filter considered: {}'.format(max_u_global))
  max_log_table = max(2**max_log2_access_size, max_u_global)
  log_table = [math.log(i) if i > 0 else 0 for i in range(0, max_log_table + 1)]

  for log2_c in range(0, max_log2_cascading + 1):
    for log2_b in range(min_log2_filter_size, max_log2_filter_size + 1):
      log2_loading = log2_b + log2_c - log2_storage_size
      if max_log2_filterset_size > 0 and log2_c + log2_b > max_log2_filterset_size:
      	continue
      c = 2**log2_c
      b = 2**log2_b
      if min_hamming_weight >= min(b, max_hamming_weight):
        continue
      avg_loading = float(2**log2_loading)
      umax_f = avg_loading + F * math.sqrt(avg_loading) #mean + F times standard dev
      log2_min_num_masks = math.log(umax_f, 2)
      umax = int(umax_f)
      log_binom_num = 0.
      log_binom_den = 0.
      precomp = 1.
      precomp_fact = 1 - 1./b
      for k in range(1, min(b, max_hamming_weight) + 1):
        p = 0.
        log_binom_num += log_table[b + 1 - k]
        log_binom_den += log_table[k]
        logmaxq = log_binom_num - log_binom_den
        #logmaxq is equal to log(binomial coefficient(b, k))
        log2maxq = logmaxq / log_table[2]
        m = max_log2_set_size(log2maxq, log2_b)
        settings_analyzed += 1
        precomp *= precomp_fact
        # precomp is equal to (1 - 1 / b)**k
        if k < min_hamming_weight or log2maxq < log2_min_num_masks or m <= 0:
          continue
        logfactorial = 0
        logavgloading = log2_loading * log_table[2]
        fk = float(k)
        precomp_pow_u = precomp
        for u in range(1, umax + 1):
          logfactorial += log_table[u]
          # logfactorial  is equal to math.log(u!)
          logproba_u = -logfactorial - avg_loading + u * logavgloading
          # logproba_u is the log-probability that u elements were inserted in before
          # current element in the same filter as current element
          # random mask case:
          # logproba_bit_coll = log-probability that the bits chosen for the
          # mask of the new value are already set because of the masks of the
          # u values already inserted
          # note: this formula is approximate, especially if k is not small compared to b.
          # indeed, it assumes that bits of previous mask values are chosen independently
          # (they are not since all bits are distinct) and that bits of the current mask
          # are chosen independently (same problem).
          logproba_bit_coll = fk * math.log(1 - precomp_pow_u)
          #logproba_bit_coll = k * math.log(1 - (1 - 1 / b)**(k * u))
          precomp_pow_u *= precomp
          p += math.exp(logproba_u  + c*logproba_bit_coll)

        results.append((log2_b, c, k, log2_loading, log2maxq, m, p))
        #print('log word size: {}; cascading: {}; hamming: {}; log2-loading: {}; log2-mask set size: {:.2f}; coll. proba for random masks: {:.2g}'.format(b, c, k, log2_loading, logmaxq, p))
  if len(results) == 0:
    print('no setting was found for the requested log storage size / for the mask set size and storage requirements')
  else:
    print('{} settings analyzed; {} eligible settings found'.format(settings_analyzed, len(results)))
    print()
    full_results = {}

    s = sorted(results, key = lambda setting:setting[6])
    full_results[s[0]] = false_positive_probability_w_maskset(s[0], F)
    target_proba = full_results[s[0]][3]

    # limit the number of calls to false_positive_probability_w_maskset by
    # computing it only for entries with collision probability for random maskset
    # below the target probability
    # (which is a collision probability with finite maskset for which we know a setting exists)

    for setting in s[1:]:
      proba_random_masks_approx = setting[6]
      # factor 2 below: we assume that the approximate probabilities computed above are
      # not off by a factor more than 2
      if proba_random_masks_approx < 2*target_proba:
        full_results[setting] = false_positive_probability_w_maskset(setting, F)
        target_proba = min(target_proba, full_results[setting][3])

    x_best = min(full_results, key=lambda x:full_results.get(x)[3])

    print('**** best result: ****')
    print_result(x_best, full_results[x_best])

def false_positive_probability_w_maskset(setting, F):
  log2_b, c, mask_weight, logloading, _, log2_mask_set_size, _ = setting
  return false_positive_proba(2**log2_b, mask_weight, float(2**logloading), log2_mask_set_size, c, False, F)

def false_positive_proba(bloom_size, mask_weight, avg_loading, log2_mask_set_size, cascading_factor, verbose = True, F=10):
  '''
  computes the false positive probabilities (f.p.p.) associated with the
  parameters provided.
  returns 4 values (a,b,c,d)
  a = f.p.p. of one bloom filter of size 'bloom_size', and an average number of
   elements equal to 'avg_loading' (the actual distribution of the loading is
   simulated), with pseudo-random masks.
  b = f.p.p. of 'cascading_factor' bloom filters with pseudo-random masks.
  c = f.p.p. of one bloom filter as before but with a finite (random) set
   of masks, of size 2^'log2_mask_set_size'.
  d = f.p.p. of  'cascading_factor' bloom filters of size 'bloom_size',
   and an average number of elements equal to 'avg_loading'.

  d is the interesting value for practical purposes. it is not equal to
   c^cascading_factor because the number of elements in all bloom filters
   are the same, hence occurrences of collisions in filters are correlated.
  '''
  max_log2_mask_set_size = logbinomial(bloom_size, mask_weight)/math.log(2)
  if  max_log2_mask_set_size < log2_mask_set_size:
    if verbose:
      print("** WARNING: log2_mask_set_size too large, reducing it to its maximum value {} **".format(max_log2_mask_set_size))
    log2_mask_set_size = max_log2_mask_set_size

  b = bloom_size
  c = cascading_factor
  proba = 0.
  proba_cascaded = 0.
  proba_finite_maskset = 0.
  proba_finite_maskset_cascaded = 0.
  umax = int(avg_loading + F * math.sqrt(avg_loading))
  mask_set_size = int(2**log2_mask_set_size)
  distinct_elts_proba = num_masks_inserted(umax, mask_set_size)
  logfact = 0.

  proba_bit_coll_exact = inclusion_probability(mask_weight, b, umax)
  proba_bit_coll_exact_distinct_masks = inclusion_probability_distinct_masks(mask_weight, b, umax)
  #s_proba_u = math.exp(-avg_loading)

  for u in range(1, umax + 1):
    logfact += math.log(u)
    #logfact is equal to log(u!)
    logproba_u = -logfact - avg_loading + u * math.log(avg_loading)
    #s_proba_u += math.exp(logproba_u)
    # log-probability that there are u elements inserted in current
    # (group of) filters, fixed-mask-set case
    proba_finite_maskset_coll_in_one_filter_with_u_elements = 0
    logproba_bit_coll = math.log(proba_bit_coll_exact[u])

    for v in range(1, min(u, mask_set_size) + 1):
      # proba_v is is the probability that v distinct masks were used in filter
      # when u values were inserted.
      proba_v = distinct_elts_proba[u][v]
      if proba_v > 0:
        if v == mask_set_size:
          proba_finite_maskset_coll_in_one_filter_with_u_elements += proba_v
        else:
          mask_coll_proba = v / mask_set_size
          proba_finite_maskset_coll_in_one_filter_with_u_elements += proba_v * \
            ((1 - mask_coll_proba) * proba_bit_coll_exact_distinct_masks[v] + mask_coll_proba)

    proba                         += math.exp(logproba_u + logproba_bit_coll)
    proba_cascaded                += math.exp(logproba_u + c * logproba_bit_coll)
    proba_u = math.exp(logproba_u)
    proba_finite_maskset          += proba_u * proba_finite_maskset_coll_in_one_filter_with_u_elements
    proba_finite_maskset_cascaded += proba_u * proba_finite_maskset_coll_in_one_filter_with_u_elements**c
  #print(s_proba_u)
  return proba, proba_cascaded, proba_finite_maskset, proba_finite_maskset_cascaded

def num_masks_inserted(umax, q):
  '''
  outputs a list l where l[u][i] (u <= umax, i <= min(umax,q)) is the probability
  to obtain i distinct elements when drawing u times an element at random in
  a fixed set of size q.
  '''
  bound = min(umax, q) + 1
  res = [[] for v in range(0, umax + 1)]
  # probabilities after 0 drawings
  output_distrib = [1. if v == 0 else 0. for v in range(0, bound)]
  res[0] = output_distrib
  for u in range(1, umax + 1):
    #probabilities after u drawings
    output_distrib_next = \
    [(output_distrib[v - 1] * (1 - (v - 1) / q) if v >= 1 else 0) +\
     output_distrib[v] * v / q for v in range(0, bound)]
    res[u] = output_distrib_next
    output_distrib = output_distrib_next
  return res

log_table = [math.log(i) if i > 0 else 0 for i in range(0, 512 + 1)]

def logbinomial(n,p):
  r = 0.
  for i in range(0, p):
    r += log_table[n - i] - log_table[i + 1]
  return r

def num_bits(k, b):
  '''outputs res[i][j] = proba to increase the hamming weight by j bits in a b-bit
  vector where i bits are already set to 1 when setting k randomly chosen bits to 1
  '''

  res = [[] for j in range(0, b + 1)]
  res[0] = [0. if j != k else 1. for j in range(0, k + 1)]

  den = logbinomial(b, k)
  for i in range(1, b + 1):
    jmin = max(0, k - i)
    jmax = min(b - i, k)
    res[i] = [0. for j in range(0, k + 1)]
    for j in range(jmin, jmax + 1):
      res[i][j] = math.exp(logbinomial(i, k - j) + logbinomial(b - i, j) - den)
  return res

def inclusion_probability(mask_weight, b, umax):
  ''' returns p s.t. for 0 <= u <= umax,
   p[u] = probability to have 'mask_weight' random bits included in the union
    of u sets of 'mask_weight' random bits, where bits are distinct in each of
   the u sets, and all bits are chosen in a b-bit set.
   The bits sets are not assumed to be distinct.
   See inclusion_probability_distinct_masks for the case where bit sets are
   distinct.
  '''
  num_added_bits_to_1_law = num_bits(mask_weight, b)
  #res[0] = 0: the first group of bits set to 1 necessarily increases the number of bits set to 1...
  res = [0. for i in range(0, umax + 1)]
  # num_bits_to_1_law is the law of the number of bits set to 1 in a b-bit word
  # when u sets of k distincts random bits were set to 1.
  # initialized below for u=0.
  # in that case, the probability computed here overestimates the actual
  # inclusion probability
  # since it does not take into account the exclusion of equal masks.
  num_bits_to_1_law = [1. if i == 0 else 0. for i in range(0, b + 1)]
  for u in range(1, umax + 1):
    new_num_bits_to_1_law = [0. for i in range(0, b + 1)]
    for i in range(0, b + 1):
      for j in range(0, min(mask_weight, b - i) + 1):
        new_num_bits_to_1_law[i + j] += num_bits_to_1_law[i] * num_added_bits_to_1_law[i][j]
    num_bits_to_1_law = new_num_bits_to_1_law
    # probability to have mask_weight bits included in u groups =
    # sum for i = 0 ... b of
    # (probability to have i distinct bits set when setting u groups) * \
    # (probability to add 0 more bits when setting mask_weight random bits
    #  assuming i bits are already set)
    for i in range(0, b + 1):
      res[u] += num_bits_to_1_law[i]*num_added_bits_to_1_law[i][0]
  return res

def inclusion_probability_distinct_masks(mask_weight, b, umax):
  ''' returns p s.t. for 0 <= u <= umax,
   p[u] = probability to have 'mask_weight' random bits included in the union
    of u sets of 'mask_weight' random bits, where bits are distinct in each of
   the u sets, and all bits are chosen in a b-bit set.
   The bit sets of 'mask_weight' bits are assumed to be pairwise distinct
   (but not disjoint).
  '''
  num_added_bits_to_1_law = num_bits(mask_weight, b)
  #res[0] = 0
  res = [0. for i in range(0, umax + 1)]
  num_bits_to_1_law = [1. if i == 0 else 0. for i in range(0, b + 1)]
  den = 1 / (math.exp(logbinomial(b, mask_weight)))
  for u in range(1, umax + 1):
    #probability to have a mask identical to the existing masks
    #there are u-1 existing masks
    p = (u-1) * den
    new_num_bits_to_1_law = [-p*num_bits_to_1_law[i] for i in range(0, b + 1)]
    for i in range(0, b + 1):
      for j in range(0, min(mask_weight, b - i) + 1):
        new_num_bits_to_1_law[i + j] += num_bits_to_1_law[i] * num_added_bits_to_1_law[i][j]
    for i in range(0, b + 1):
        new_num_bits_to_1_law[i] *= 1/(1-p)
    num_bits_to_1_law = new_num_bits_to_1_law
    for i in range(0, b + 1):
      res[u] += num_bits_to_1_law[i]*num_added_bits_to_1_law[i][0]
    pprime = u * den
    res[u] -= pprime
    res[u] *= 1/(1-pprime)
  return res


def print_search_parameters(log2_storage_size, max_log2_access_size, \
    max_log2_mask_set_size, max_log2_mask_storage_bit_size, min_hamming_weight, max_hamming_weight, max_log2_filter_size, max_log2_cascading):
  print('** search parameters **')
  print('storage size per element: 2^{} = {}'.format(log2_storage_size, 2**log2_storage_size))
  print('Maximum access size to check one element: 2^{} = {} bits = {} bytes'.format(\
      max_log2_access_size, 2**max_log2_access_size, 2**(max_log2_access_size-3)))
  if max_log2_mask_set_size > 0:
   print('Maximum mask set size: : 2^{}'.format(max_log2_mask_set_size))
  if max_log2_mask_storage_bit_size > 0:
   print('Maximum mask storage size: : 2^{} bits'.format(max_log2_mask_storage_bit_size))
  if min_hamming_weight > 0:
    print('Minimum Hamming weight of masks: {}'.format(min_hamming_weight))
  if max_hamming_weight > 0:
    print('Maximum Hamming weight of masks: {}'.format(max_hamming_weight))
  print('Maximum filter size: ',end='')
  if max_log2_filter_size > 0:
    print('2^{} = {} bits = {} bytes'.format(max_log2_filter_size, 2**max_log2_filter_size, 2**(max_log2_filter_size-3)))
  else:
    print('same as access size')
  if max_log2_cascading > 0:
    print('Maximum number of iterations: 2^{} = {}'.format(max_log2_cascading, 2**max_log2_cascading))
  print()

def print_result(setting, probas):
  _, p_cascaded, proba_finite_q, proba_finite_q_cascaded = probas
  log2_b, c, k, log2_loading, log2_full_max_set_size, log2_mask_size, approx_proba_with_random_masks = setting
  b = 2**log2_b
  print('\ncollision probability')
  print(' random masks: {:.3e}'.format(p_cascaded))
  print(' finite mask set: {:.3e}'.format(proba_finite_q_cascaded))

  if c > 1:
    print(' [in one filter, with finite mask set: {:.3e}]'.format(proba_finite_q))

  print('\nsettings:')
  print(' word size: {}'.format(b))
  print(' cascading factor: {}'.format(c))
  print(' loading factor (average number of values per filter): 2^{} = {}'.format(log2_loading, 2**log2_loading))
  print(' hamming weight of masks: {}'.format(k))
  print(' mask set size: 2^{:.3f} <= binomial({}, {}) = 2^{:.3f}'.format( \
    log2_mask_size, b, k, log2_full_max_set_size))
  storage_size = b * c // 2**log2_loading
  print(' storage size = word size * cascading factor / loading = {} * {} / {} = {}'.format(\
    b, c, 2**log2_loading, storage_size))

  print('\nMisc info:')
  print(' Mask set log storage size, in bits: {:.2f}; in bytes: {:.2f}'.format(log2_mask_size + log2_b, log2_mask_size + log2_b - 3))
  print(' random bits required per element insertion: {:.2f}, i.e. {:.2f} bytes'.format(log2_mask_size * c, (log2_mask_size * c)/8))

  opt_hashes = int(storage_size * math.log(2) + 0.5)
  print('\nBloom Filter comparison:')
  print(' optimal number of hash functions for a unique BF with the same storage size per value: {}'.format(opt_hashes))
  p_bloom = (1-math.exp(-opt_hashes/storage_size))**opt_hashes
  print(' collision probability of this optimal BF: {:.3e}'.format(p_bloom))
  print(' practical vs optimal collision probability ratio: {:.3g}'.format(proba_finite_q_cascaded / p_bloom))
  print()

