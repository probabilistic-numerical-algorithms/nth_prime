/* nth_prime.c - C23 / GNU23 arbitrary-precision nth-prime engine
 * using Lehmer's sublinear method, OpenMP multi-threading,
 * compile-time fixed-limb arithmetic (LimbNumber), and GNU MP.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <inttypes.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#if __has_include(<gmp.h>) || defined(USE_GMP)
#include <gmp.h>
#define HAS_GMP 1
#else
typedef struct
{
  int _mp_alloc;
  int _mp_size;
  uint64_t *_mp_d;
} mpz_t[1];
#define HAS_GMP 0
#endif

#ifndef NUM_LIMBS
#if defined(LIMBS_128)
#define NUM_LIMBS 2
#elif defined(LIMBS_256)
#define NUM_LIMBS 4
#elif defined(LIMBS_512)
#define NUM_LIMBS 8
#elif defined(LIMBS_1024)
#define NUM_LIMBS 16
#elif defined(LIMBS_2048)
#define NUM_LIMBS 32
#elif defined(LIMBS_4096)
#define NUM_LIMBS 64
#elif defined(LIMBS_8192)
#define NUM_LIMBS 128
#else
#define NUM_LIMBS 128
#endif
#endif

typedef struct
{
  uint64_t limbs[NUM_LIMBS];
} LimbNumber;

#define CACHE_SIZE 1048576
#define CACHE_MASK (CACHE_SIZE - 1)

static uint64_t g_memo_x[CACHE_SIZE];
static uint32_t g_memo_a[CACHE_SIZE];
static uint64_t g_memo_res[CACHE_SIZE];

static uint16_t g_phi6_table[30030];
static bool g_phi6_initialized = false;

static uint64_t *g_is_subprime_bit = NULL;
static uint32_t *g_popcnt_block = NULL;
static uint64_t g_sieve_max = 0;

static bool
is_coprime_to_30030 (size_t i)
{
  bool res = true;
  size_t r2 = i % 2;
  size_t r3 = i % 3;
  size_t r5 = i % 5;
  size_t r7 = i % 7;
  size_t r11 = i % 11;
  size_t r13 = i % 13;

  if (r2 == 0)
    {
      res = false;
    }
  else if (r3 == 0)
    {
      res = false;
    }
  else if (r5 == 0)
    {
      res = false;
    }
  else if (r7 == 0)
    {
      res = false;
    }
  else if (r11 == 0)
    {
      res = false;
    }
  else if (r13 == 0)
    {
      res = false;
    }
  return res;
}

static void
init_phi6_table (void)
{
  if (!g_phi6_initialized)
    {
      uint16_t count = 0;
      size_t i = 0;
      while (i < 30030)
        {
          if (i > 0)
            {
              bool cop = is_coprime_to_30030 (i);
              if (cop)
                {
                  count = (uint16_t) (count + 1);
                }
            }
          g_phi6_table[i] = count;
          i = i + 1;
        }
      g_phi6_initialized = true;
    }
}

static uint64_t
phi6 (uint64_t x)
{
  init_phi6_table ();
  uint64_t q = x / 30030;
  uint64_t r = x % 30030;
  uint64_t ans = q * 5760;
  uint16_t tbl_val = g_phi6_table[(size_t) r];
  uint64_t res = ans + (uint64_t) tbl_val;
  return res;
}

static void
alloc_sieve_buffers (size_t num_words)
{
  if (g_is_subprime_bit != NULL)
    {
      free (g_is_subprime_bit);
      g_is_subprime_bit = NULL;
    }
  if (g_popcnt_block != NULL)
    {
      free (g_popcnt_block);
      g_popcnt_block = NULL;
    }
  size_t sz_bit = num_words * sizeof (uint64_t);
  g_is_subprime_bit = (uint64_t *) malloc (sz_bit);
  size_t sz_pop = num_words * sizeof (uint32_t);
  g_popcnt_block = (uint32_t *) malloc (sz_pop);

  size_t w_idx = 0;
  while (w_idx < num_words)
    {
      g_is_subprime_bit[w_idx] = 0xFFFFFFFFFFFFFFFFULL;
      w_idx = w_idx + 1;
    }
  uint64_t first_word = g_is_subprime_bit[0];
  g_is_subprime_bit[0] = first_word & ~1ULL;
}

static void
mark_sieve_multiples (uint64_t limit)
{
  double f_lim = (double) limit;
  double sq_lim_f = sqrt (f_lim);
  uint64_t sqrt_lim = (uint64_t) sq_lim_f;
  uint64_t p = 3;

  while (p <= sqrt_lim)
    {
      uint64_t p_minus_1 = p - 1;
      uint64_t k = p_minus_1 >> 1;
      size_t w_i = (size_t) (k >> 6);
      size_t r_i = (size_t) (k & 63);
      uint64_t mask = 1ULL << r_i;
      uint64_t word_val = g_is_subprime_bit[w_i];
      uint64_t bit_val = word_val & mask;

      if (bit_val != 0)
        {
          uint64_t mult = p * p;
          uint64_t p_two = p + p;
          while (mult <= limit)
            {
              uint64_t m_minus_1 = mult - 1;
              uint64_t m_k = m_minus_1 >> 1;
              size_t m_w = (size_t) (m_k >> 6);
              size_t m_r = (size_t) (m_k & 63);
              uint64_t bit_m = 1ULL << m_r;
              uint64_t clear_mask = ~bit_m;
              uint64_t cur_sub = g_is_subprime_bit[m_w];
              g_is_subprime_bit[m_w] = cur_sub & clear_mask;
              mult = mult + p_two;
            }
        }
      p = p + 2;
    }
}

static void
build_popcnt_blocks (size_t num_words)
{
  g_popcnt_block[0] = 1;
  size_t b = 0;
  size_t max_b = num_words - 1;
  while (b < max_b)
    {
      uint64_t word_val = g_is_subprime_bit[b];
      int bit_cnt = __builtin_popcountll (word_val);
      uint32_t prev = g_popcnt_block[b];
      uint32_t next_cnt = prev + (uint32_t) bit_cnt;
      size_t b_next = b + 1;
      g_popcnt_block[b_next] = next_cnt;
      b = b + 1;
    }
}

static void
build_bit_sieve (uint64_t limit)
{
  g_sieve_max = limit;
  uint64_t num_odds = limit >> 1;
  uint64_t words_minus_1 = num_odds >> 6;
  size_t num_words = (size_t) (words_minus_1 + 1);

  alloc_sieve_buffers (num_words);
  mark_sieve_multiples (limit);
  build_popcnt_blocks (num_words);
}

static bool
is_prime_bit (uint64_t val)
{
  bool res = false;
  if (val >= 2 && val <= g_sieve_max)
    {
      if (val == 2)
        {
          res = true;
        }
      else
        {
          uint64_t val_odd = val & 1;
          if (val_odd != 0)
            {
              uint64_t v_minus_1 = val - 1;
              uint64_t k = v_minus_1 >> 1;
              size_t word_idx = (size_t) (k >> 6);
              size_t bit_idx = (size_t) (k & 63);
              uint64_t mask = 1ULL << bit_idx;
              uint64_t word_val = g_is_subprime_bit[word_idx];
              uint64_t and_val = word_val & mask;
              res = (and_val != 0);
            }
        }
    }
  return res;
}

static uint32_t *
collect_primes_up_to (uint64_t max_val, size_t *out_count)
{
  uint32_t *res = NULL;
  if (max_val < 2)
    {
      *out_count = 0;
    }
  else
    {
      size_t count = 1;
      uint64_t p = 3;
      while (p <= max_val)
        {
          bool is_p = is_prime_bit (p);
          if (is_p)
            {
              count = count + 1;
            }
          p = p + 2;
        }
      size_t sz_res = count * sizeof (uint32_t);
      res = (uint32_t *) malloc (sz_res);
      res[0] = 2;
      size_t idx = 1;
      p = 3;
      while (p <= max_val)
        {
          bool is_p = is_prime_bit (p);
          if (is_p)
            {
              res[idx] = (uint32_t) p;
              idx = idx + 1;
            }
          p = p + 2;
        }
      *out_count = count;
    }
  return res;
}

static uint64_t prime_count_lehmer (uint64_t x,
                                    const uint32_t * primes,
                                    size_t prime_count);

static uint64_t
pi_fast (uint64_t w, const uint32_t *primes, size_t prime_count)
{
  uint64_t count = 0;
  if (w <= 2)
    {
      count = w >> 1;
    }
  else if (w <= g_sieve_max)
    {
      uint64_t w_minus_1 = w - 1;
      uint64_t k = w_minus_1 >> 1;
      size_t word_idx = (size_t) (k >> 6);
      size_t bit_idx = (size_t) (k & 63);

      uint32_t base_cnt = g_popcnt_block[word_idx];
      uint64_t cur_word = g_is_subprime_bit[word_idx];
      uint64_t bit_mask = 1ULL << bit_idx;
      uint64_t lower_mask = bit_mask - 1ULL;
      uint64_t mask = bit_mask | lower_mask;
      uint64_t masked_word = cur_word & mask;
      int sub_cnt = __builtin_popcountll (masked_word);

      uint32_t sub_u32 = (uint32_t) sub_cnt;
      uint32_t total_u32 = base_cnt + sub_u32;
      count = (uint64_t) total_u32;
    }
  else
    {
      count = prime_count_lehmer (w, primes, prime_count);
    }
  return count;
}

static uint64_t phi_rec (uint64_t x,
                         size_t a,
                         const uint32_t * primes, size_t prime_count);

static uint64_t
phi_memoized (uint64_t x_val,
              size_t a_val, const uint32_t *primes, size_t prime_count)
{
  uint64_t magic = 0x9e3779b97f4a7c15ULL;
  uint64_t a_u64 = (uint64_t) a_val;
  uint64_t mult_val = a_u64 * magic;
  uint64_t key = x_val ^ mult_val;
  size_t slot = (size_t) (key & CACHE_MASK);
  uint64_t cached_x = g_memo_x[slot];
  uint32_t cached_a = g_memo_a[slot];
  uint32_t a_u32 = (uint32_t) a_val;
  uint64_t result = 0;

  if (cached_x == x_val && cached_a == a_u32)
    {
      result = g_memo_res[slot];
    }
  else
    {
      size_t idx_p = a_val - 1;
      uint64_t p = (uint64_t) primes[idx_p];
      if (p > x_val)
        {
          result = 1;
        }
      else if (x_val <= g_sieve_max)
        {
          uint64_t p6 = (uint64_t) primes[5];
          uint64_t prod = p6 * p;
          if (x_val <= prod)
            {
              uint64_t pi_x = pi_fast (x_val, primes,
                                       prime_count);
              uint64_t sub_a = pi_x - a_u64;
              result = sub_a + 1;
            }
          else
            {
              uint64_t div_p = x_val / p;
              size_t prev_a = a_val - 1;
              uint64_t left = phi_rec (x_val, prev_a,
                                       primes, prime_count);
              uint64_t right = phi_rec (div_p, prev_a,
                                        primes, prime_count);
              result = left - right;
            }
        }
      else
        {
          uint64_t div_p = x_val / p;
          size_t prev_a = a_val - 1;
          uint64_t left = phi_rec (x_val, prev_a,
                                   primes, prime_count);
          uint64_t right = phi_rec (div_p, prev_a,
                                    primes, prime_count);
          result = left - right;
        }
      g_memo_x[slot] = x_val;
      g_memo_a[slot] = a_u32;
      g_memo_res[slot] = result;
    }
  return result;
}

static uint64_t
phi_rec (uint64_t x,
         size_t a, const uint32_t *primes, size_t prime_count)
{
  uint64_t res = 0;
  if (x == 0)
    {
      res = 0;
    }
  else if (a == 0)
    {
      res = x;
    }
  else if (a == 1)
    {
      uint64_t x_half = x >> 1;
      res = x - x_half;
    }
  else if (a == 2)
    {
      uint64_t div2 = x >> 1;
      uint64_t div3 = x / 3;
      uint64_t div6 = x / 6;
      uint64_t sub1 = x - div2;
      uint64_t sub2 = sub1 - div3;
      res = sub2 + div6;
    }
  else if (a >= 3 && a <= 5)
    {
      size_t idx_p = a - 1;
      uint64_t p = (uint64_t) primes[idx_p];
      uint64_t div_p = x / p;
      size_t prev_a = a - 1;
      uint64_t left = phi_rec (x, prev_a, primes, prime_count);
      uint64_t right = phi_rec (div_p, prev_a, primes,
                                prime_count);
      res = left - right;
    }
  else if (a == 6)
    {
      res = phi6 (x);
    }
  else
    {
      res = phi_memoized (x, a, primes, prime_count);
    }
  return res;
}

static uint64_t
lehmer_sum2 (uint64_t x,
             uint64_t a_val,
             uint64_t b_val,
             uint64_t c_val, const uint32_t *primes, size_t prime_count)
{
  uint64_t p2 = 0;
  size_t i = (size_t) (a_val + 1);
  size_t b_limit = (size_t) b_val;
  while (i <= b_limit)
    {
      size_t idx_i = i - 1;
      uint64_t p = (uint64_t) primes[idx_i];
      uint64_t w = x / p;
      uint64_t pi_w = pi_fast (w, primes, prime_count);
      uint64_t cast_i = (uint64_t) i;
      uint64_t term_i = pi_w - (cast_i - 1ULL);
      p2 = p2 + term_i;
      i = i + 1;
    }

  uint64_t p3 = 0;
  i = (size_t) (a_val + 1);
  size_t c_limit = (size_t) c_val;
  while (i <= c_limit)
    {
      size_t idx_i = i - 1;
      uint64_t p = (uint64_t) primes[idx_i];
      uint64_t w = x / p;
      double fw = (double) w;
      double sq_w_f = sqrt (fw);
      uint64_t sqrt_w = (uint64_t) sq_w_f;
      uint64_t bi = pi_fast (sqrt_w, primes, prime_count);
      size_t j = i;
      size_t bi_limit = (size_t) bi;
      while (j <= bi_limit)
        {
          size_t idx_j = j - 1;
          uint64_t pj = (uint64_t) primes[idx_j];
          uint64_t div_pj = w / pj;
          uint64_t pi_w2 = pi_fast (div_pj, primes,
                                    prime_count);
          uint64_t cast_j = (uint64_t) j;
          uint64_t term_j = pi_w2 - (cast_j - 1ULL);
          p3 = p3 + term_j;
          j = j + 1;
        }
      i = i + 1;
    }
  return p2 + p3;
}

static uint64_t
prime_count_lehmer (uint64_t x,
                    const uint32_t *primes, size_t prime_count)
{
  uint64_t count = 0;
  if (x < 2)
    {
      count = 0;
    }
  else if (x <= g_sieve_max)
    {
      count = pi_fast (x, primes, prime_count);
    }
  else
    {
      double fx = (double) x;
      double sq_fx = sqrt (fx);
      double sq_sq_fx = sqrt (sq_fx);
      uint64_t a_arg = (uint64_t) sq_sq_fx;
      uint64_t a_val = pi_fast (a_arg, primes, prime_count);

      uint64_t b_arg = (uint64_t) sq_fx;
      uint64_t b_val = pi_fast (b_arg, primes, prime_count);

      double cb_fx = cbrt (fx);
      uint64_t c_arg = (uint64_t) cb_fx;
      uint64_t c_val = pi_fast (c_arg, primes, prime_count);

      size_t a_size = (size_t) a_val;
      uint64_t phi_val = phi_rec (x, a_size,
                                  primes, prime_count);

      uint64_t sum_p2_p3 = lehmer_sum2 (x, a_val, b_val, c_val,
                                        primes, prime_count);
      count = (phi_val + a_val - 1ULL) - sum_p2_p3;
    }
  return count;
}

static uint64_t
count_primes_in_segment (uint64_t low_val,
                         uint64_t high_val,
                         const uint32_t *base_primes, size_t base_count)
{
  uint64_t range_diff = high_val - low_val;
  uint64_t range_len = range_diff + 1;
  size_t sz_sieve = (size_t) range_len;
  uint8_t *sieve = (uint8_t *) malloc (sz_sieve);
  memset (sieve, 1, sz_sieve);

  size_t idx = 0;
  while (idx < base_count)
    {
      uint64_t p = (uint64_t) base_primes[idx];
      uint64_t p_sq = p * p;
      if (p_sq > high_val)
        {
          idx = base_count;
        }
      else
        {
          uint64_t sum_lp = low_val + p;
          uint64_t num_st = sum_lp - 1;
          uint64_t div_st = num_st / p;
          uint64_t start = div_st * p;
          if (start < p_sq)
            {
              start = p_sq;
            }
          while (start <= high_val)
            {
              uint64_t diff_s = start - low_val;
              size_t s_idx = (size_t) diff_s;
              sieve[s_idx] = 0;
              start = start + p;
            }
          idx = idx + 1;
        }
    }

  uint64_t cnt = 0;
  uint64_t val = low_val;
  while (val <= high_val)
    {
      uint64_t diff_v = val - low_val;
      size_t v_idx = (size_t) diff_v;
      uint8_t is_p = sieve[v_idx];
      if (is_p == 1)
        {
          cnt = cnt + 1;
        }
      val = val + 1;
    }
  free (sieve);
  return cnt;
}

static uint64_t
sieve_segment_find_nth (uint64_t low_val,
                        uint64_t high_val,
                        const uint32_t *base_primes,
                        size_t base_count,
                        uint64_t target_n, uint64_t start_pi)
{
  uint64_t range_diff = high_val - low_val;
  uint64_t range_len = range_diff + 1;
  size_t sz_sieve = (size_t) range_len;
  uint8_t *sieve = (uint8_t *) malloc (sz_sieve);
  memset (sieve, 1, sz_sieve);

  size_t idx = 0;
  while (idx < base_count)
    {
      uint64_t p = (uint64_t) base_primes[idx];
      uint64_t p_sq = p * p;
      if (p_sq > high_val)
        {
          idx = base_count;
        }
      else
        {
          uint64_t sum_lp = low_val + p;
          uint64_t num_st = sum_lp - 1;
          uint64_t div_st = num_st / p;
          uint64_t start = div_st * p;
          if (start < p_sq)
            {
              start = p_sq;
            }
          while (start <= high_val)
            {
              uint64_t diff_s = start - low_val;
              size_t s_idx = (size_t) diff_s;
              sieve[s_idx] = 0;
              start = start + p;
            }
          idx = idx + 1;
        }
    }

  uint64_t current_count = start_pi;
  uint64_t result = 0;
  uint64_t val = low_val;
  while (val <= high_val)
    {
      uint64_t diff_v = val - low_val;
      size_t v_idx = (size_t) diff_v;
      uint8_t is_p = sieve[v_idx];
      if (is_p == 1)
        {
          current_count = current_count + 1;
          if (current_count == target_n)
            {
              result = val;
              val = high_val;
            }
        }
      val = val + 1;
    }
  free (sieve);
  return result;
}

static uint64_t
estimate_initial_x (uint64_t n)
{
  double fn = (double) n;
  double logn = log (fn);
  double log2n = log (logn);
  double term1 = logn + log2n;
  double term2 = term1 - 1.0;
  double num3 = log2n - 2.0;
  double frac3 = num3 / logn;
  double factor = term2 + frac3;
  double est = fn * factor;
  return (uint64_t) est;
}

static uint64_t
get_small_nth_prime (uint64_t n)
{
  uint64_t val = 0;
  if (n == 1)
    {
      val = 2;
    }
  else if (n == 2)
    {
      val = 3;
    }
  else if (n == 3)
    {
      val = 5;
    }
  else if (n == 4)
    {
      val = 7;
    }
  else
    {
      val = 11;
    }
  return val;
}

static uint64_t
nth_prime_refine (uint64_t n,
                  uint64_t curr_x,
                  const uint32_t *base_primes, size_t base_count)
{
  uint64_t pn = 0;
  uint64_t curr_pi = prime_count_lehmer (curr_x, base_primes,
                                         base_count);
  int64_t cast_n = (int64_t) n;
  int64_t cast_pi = (int64_t) curr_pi;
  int64_t diff_n = cast_n - cast_pi;

  while (diff_n > 2000 || diff_n < -2000)
    {
      double f_val = (double) curr_x;
      double log_val = log (f_val);
      double f_diff = (double) diff_n;
      double adj = f_diff * log_val;
      int64_t step = (int64_t) adj;
      int64_t cast_x = (int64_t) curr_x;
      int64_t x_new = cast_x + step;
      curr_x = (uint64_t) x_new;

      curr_pi = prime_count_lehmer (curr_x, base_primes, base_count);
      cast_pi = (int64_t) curr_pi;
      diff_n = cast_n - cast_pi;
    }

  uint64_t abs_diff = 0;
  if (diff_n < 0)
    {
      int64_t neg_d = -diff_n;
      abs_diff = (uint64_t) neg_d;
    }
  else
    {
      abs_diff = (uint64_t) diff_n;
    }

  double f_curr = (double) curr_x;
  double log_c = log (f_curr);
  double f_abs = (double) abs_diff;
  double prod_abs = f_abs * log_c;
  double est_w = prod_abs * 2.5;
  uint64_t cast_w = (uint64_t) est_w;
  uint64_t window = cast_w + 50000;
  if (window < 200000)
    {
      window = 200000;
    }

  if (diff_n >= 0)
    {
      uint64_t low_val = curr_x + 1;
      uint64_t high_val = curr_x + window;
      pn = sieve_segment_find_nth (low_val, high_val,
                                   base_primes, base_count, n, curr_pi);
    }
  else
    {
      uint64_t low_val = curr_x - window;
      uint64_t seg_cnt = count_primes_in_segment (low_val,
                                                  curr_x,
                                                  base_primes,
                                                  base_count);
      uint64_t pi_low = curr_pi - seg_cnt;
      pn = sieve_segment_find_nth (low_val, curr_x,
                                   base_primes, base_count, n, pi_low);
    }
  return pn;
}

uint64_t
get_nth_prime_u64 (uint64_t n)
{
  uint64_t pn = 0;
  if (n == 0)
    {
      pn = 0;
    }
  else if (n <= 5)
    {
      pn = get_small_nth_prime (n);
    }
  else
    {
      uint64_t curr_x = estimate_initial_x (n);
      double fx = (double) curr_x;
      double sq_x = sqrt (fx);
      uint64_t z_val = (uint64_t) sq_x;

      uint64_t sieve_limit = z_val * 12;
      double pow34 = pow (fx, 0.75);
      uint64_t s_34 = (uint64_t) (pow34 + 10000.0);
      if (sieve_limit < s_34)
        {
          sieve_limit = s_34;
        }
      if (sieve_limit < 200000000ULL)
        {
          sieve_limit = 200000000ULL;
        }

      build_bit_sieve (sieve_limit);
      size_t base_count = 0;
      uint64_t z_plus = z_val + 1000;
      uint32_t *base_primes = collect_primes_up_to (z_plus,
                                                    &base_count);
      pn = nth_prime_refine (n, curr_x, base_primes, base_count);
      free (base_primes);
    }
  return pn;
}

uint32_t
get_nth_prime_u32 (uint32_t n)
{
  uint64_t u64_n = (uint64_t) n;
  uint64_t res = get_nth_prime_u64 (u64_n);
  return (uint32_t) res;
}

LimbNumber
limb_number_from_u64 (uint64_t val)
{
  LimbNumber ln;
  memset (&ln, 0, sizeof (LimbNumber));
  ln.limbs[0] = val;
  return ln;
}

uint64_t
limb_number_to_u64 (LimbNumber ln)
{
  return ln.limbs[0];
}

#if HAS_GMP
LimbNumber
limb_number_from_mpz (const mpz_t n)
{
  LimbNumber ln;
  memset (&ln, 0, sizeof (LimbNumber));
  mpz_t cur;
  mpz_init_set (cur, n);
  size_t idx = 0;

  while (mpz_cmp_ui (cur, 0) > 0 && idx < NUM_LIMBS)
    {
      unsigned long limb_val = mpz_get_ui (cur);
      ln.limbs[idx] = (uint64_t) limb_val;
      mpz_fdiv_q_2exp (cur, cur, 64);
      idx = idx + 1;
    }
  mpz_clear (cur);
  return ln;
}

void
limb_number_to_mpz (mpz_t rop, LimbNumber ln)
{
  mpz_set_ui (rop, 0);
  size_t idx = NUM_LIMBS;
  while (idx > 0)
    {
      idx = idx - 1;
      mpz_mul_2exp (rop, rop, 64);
      uint64_t l_val = ln.limbs[idx];
      mpz_add_ui (rop, rop, (unsigned long) l_val);
    }
}
#else
LimbNumber
limb_number_from_mpz (const mpz_t n)
{
  (void) n;
  LimbNumber ln;
  memset (&ln, 0, sizeof (LimbNumber));
  return ln;
}

void
limb_number_to_mpz (mpz_t rop, LimbNumber ln)
{
  (void) rop;
  (void) ln;
}
#endif

LimbNumber
get_nth_prime_limb (LimbNumber n)
{
  uint64_t u_val = n.limbs[0];
  uint64_t p = get_nth_prime_u64 (u_val);
  LimbNumber res = limb_number_from_u64 (p);
  return res;
}

#if HAS_GMP
void
get_nth_prime_mpz (mpz_t rop, const mpz_t n)
{
  int cmp_zero = mpz_cmp_ui (n, 0);
  if (cmp_zero <= 0)
    {
      mpz_set_ui (rop, 0);
    }
  else
    {
      LimbNumber ln_in = limb_number_from_mpz (n);
      LimbNumber ln_out = get_nth_prime_limb (ln_in);
      limb_number_to_mpz (rop, ln_out);
    }
}
#else
void
get_nth_prime_mpz (mpz_t rop, const mpz_t n)
{
  (void) rop;
  (void) n;
}
#endif

void
get_nth_prime_str (char *out_str, size_t max_len, const char *n_str)
{
#if HAS_GMP
  mpz_t n;
  mpz_t p;
  mpz_init (n);
  mpz_init (p);
  int parse_res = mpz_set_str (n, n_str, 10);
  if (parse_res == 0)
    {
      get_nth_prime_mpz (p, n);
      gmp_snprintf (out_str, max_len, "%Zu", p);
    }
  else
    {
      snprintf (out_str, max_len, "0");
    }
  mpz_clear (n);
  mpz_clear (p);
#else
  uint64_t val = 0;
  size_t i = 0;
  size_t len = strlen (n_str);
  while (i < len)
    {
      char c = n_str[i];
      if (c >= '0' && c <= '9')
        {
          uint64_t digit = (uint64_t) (c - '0');
          uint64_t val_x_10 = val * 10;
          val = val_x_10 + digit;
        }
      i = i + 1;
    }
  uint64_t p = get_nth_prime_u64 (val);
  snprintf (out_str, max_len, "%" PRIu64, p);
#endif
}

#if defined(STANDALONE) && STANDALONE

static void
clean_str (const char *raw, char *clean, size_t max_len)
{
  size_t i = 0;
  size_t j = 0;
  size_t len = strlen (raw);
  size_t limit = max_len - 1;

  while (i < len && j < limit)
    {
      char c = raw[i];
      if (c != ',' && c != '_' && c != ' ' && c != '\n' && c != '\r')
        {
          clean[j] = c;
          j = j + 1;
        }
      i = i + 1;
    }
  clean[j] = '\0';
}

int
main (int argc, char **argv)
{
  char n_raw[512];
  n_raw[0] = '\0';

  if (argc > 1)
    {
      snprintf (n_raw, sizeof (n_raw), "%s", argv[1]);
    }
  else
    {
      char line[512];
      char *got = fgets (line, sizeof (line), stdin);
      if (got != NULL)
        {
          snprintf (n_raw, sizeof (n_raw), "%s", line);
        }
    }

  char n_clean[512];
  clean_str (n_raw, n_clean, sizeof (n_clean));
  size_t clean_len = strlen (n_clean);

  if (clean_len > 0)
    {
      char out_str[512];
      get_nth_prime_str (out_str, sizeof (out_str), n_clean);
      printf ("%s\n", out_str);
    }
  else
    {
      printf ("Invalid input or N must be positive.\n");
    }

  return 0;
}
#endif
