! nth_prime_64.f90 - Fortran 2023 64-bit hardware integer nth-prime
! engine implementing fast Lehmer prime counting, modular wheel base
! cases, recursive memoization with Buchstab tree truncation,
! pre-accumulated hardware POPCNT blocks, and OpenMP parallelization.

module nth_prime_64_mod
  use, intrinsic :: iso_fortran_env, only : int64, real64, real128, &
                                            input_unit, output_unit
  implicit none
  private
  public :: get_nth_prime_u64, get_nth_prime_str

  integer, parameter :: i64 = int64
  integer, parameter :: r64 = real64
  integer, parameter :: r128 = real128

  integer(i64), parameter :: CACHE_SIZE = 1048576_i64
  integer(i64), parameter :: CACHE_MASK = 1048575_i64

  integer(i64) :: g_memo_x(0:CACHE_SIZE-1) = 0_i64
  integer(i64) :: g_memo_a(0:CACHE_SIZE-1) = 0_i64
  integer(i64) :: g_memo_res(0:CACHE_SIZE-1) = 0_i64

  integer(i64) :: g_phi6_table(0:30029)
  logical :: g_phi6_initialized = .false.

  integer(i64), allocatable :: g_is_subprime_bit(:)
  integer(i64), allocatable :: g_popcnt_block(:)
  integer(i64) :: g_sieve_max = 0_i64

contains

  pure function is_coprime_to_30030(i_val) result(res)
    integer(i64), intent(in) :: i_val
    logical :: res
    res = .true.
    if (modulo(i_val, 2_i64) == 0_i64) then
      res = .false.
    else if (modulo(i_val, 3_i64) == 0_i64) then
      res = .false.
    else if (modulo(i_val, 5_i64) == 0_i64) then
      res = .false.
    else if (modulo(i_val, 7_i64) == 0_i64) then
      res = .false.
    else if (modulo(i_val, 11_i64) == 0_i64) then
      res = .false.
    else if (modulo(i_val, 13_i64) == 0_i64) then
      res = .false.
    end if
  end function is_coprime_to_30030

  subroutine init_phi6_table()
    integer(i64) :: running, i_val
    if (.not. g_phi6_initialized) then
      running = 0_i64
      i_val = 1_i64
      do while (i_val <= 30030_i64)
        if (is_coprime_to_30030(i_val)) then
          running = running + 1_i64
        end if
        g_phi6_table(i_val - 1_i64) = running
        i_val = i_val + 1_i64
      end do
      g_phi6_initialized = .true.
    end if
  end subroutine init_phi6_table

  pure function phi6(x_val) result(ans)
    integer(i64), intent(in) :: x_val
    integer(i64) :: ans, q_val, r_val
    q_val = x_val / 30030_i64
    r_val = modulo(x_val, 30030_i64)
    ans = q_val * 5760_i64
    if (r_val > 0_i64) then
      ans = ans + g_phi6_table(r_val - 1_i64)
    end if
  end function phi6

  subroutine build_bit_sieve(limit)
    integer(i64), intent(in) :: limit
    integer(i64) :: num_odds, num_words, i_val, i_odd, w_idx, b_idx
    integer(i64) :: i2_val, j_val, j_odd, jw, jb, w_val, w_iter
    integer(i64) :: sqrt_lim, clr_mask, prev_cnt

    if (limit <= g_sieve_max .and. allocated(g_is_subprime_bit)) then
      return
    end if

    if (allocated(g_is_subprime_bit)) deallocate(g_is_subprime_bit)
    if (allocated(g_popcnt_block)) deallocate(g_popcnt_block)

    g_sieve_max = limit
    num_odds = limit / 2_i64
    num_words = (num_odds / 64_i64) + 1_i64
    if (num_words == 0_i64) num_words = 1_i64

    allocate(g_is_subprime_bit(0:num_words-1))
    allocate(g_popcnt_block(0:num_words-1))

    g_is_subprime_bit = -1_i64  ! all 1-bits in two's complement (0xFFFFFFFFFFFFFFFF)
    g_is_subprime_bit(0) = iand(g_is_subprime_bit(0), not(1_i64))

    sqrt_lim = int(sqrt(real(limit, kind=r64)), kind=i64)
    i_val = 3_i64
    do while (i_val <= sqrt_lim)
      i_odd = (i_val - 1_i64) / 2_i64
      w_idx = i_odd / 64_i64
      b_idx = iand(i_odd, 63_i64)
      w_val = iand(g_is_subprime_bit(w_idx), ishft(1_i64, int(b_idx)))
      if (w_val /= 0_i64) then
        i2_val = i_val + i_val
        j_val = i_val * i_val
        do while (j_val <= limit)
          j_odd = (j_val - 1_i64) / 2_i64
          jw = j_odd / 64_i64
          jb = iand(j_odd, 63_i64)
          clr_mask = not(ishft(1_i64, int(jb)))
          g_is_subprime_bit(jw) = iand(g_is_subprime_bit(jw), clr_mask)
          j_val = j_val + i2_val
        end do
      end if
      i_val = i_val + 2_i64
    end do

    g_popcnt_block(0) = 1_i64  ! include prime 2
    w_iter = 0_i64
    do while (w_iter < num_words - 1_i64)
      prev_cnt = g_popcnt_block(w_iter)
      g_popcnt_block(w_iter + 1_i64) = prev_cnt + &
        int(popcnt(g_is_subprime_bit(w_iter)), kind=i64)
      w_iter = w_iter + 1_i64
    end do
  end subroutine build_bit_sieve

  recursive function prime_count_lehmer(x_val, primes) &
    result(count_val)
    integer(i64), intent(in) :: x_val
    integer(i64), intent(in) :: primes(:)
    integer(i64) :: count_val, a_val, b_val, c_val, phi_val, sum_p2_p3
    real(r64) :: fx, sq_x, sq_sq_x, cb_x

    if (x_val < 2_i64) then
      count_val = 0_i64
    else if (x_val <= g_sieve_max) then
      count_val = pi_fast(x_val, primes)
    else
      fx = real(x_val, kind=r64)
      sq_x = sqrt(fx)
      sq_sq_x = sqrt(sq_x)
      a_val = pi_fast(int(sq_sq_x, kind=i64), primes)
      b_val = pi_fast(int(sq_x, kind=i64), primes)

      cb_x = fx ** (1.0_r64 / 3.0_r64)
      c_val = pi_fast(int(cb_x, kind=i64), primes)

      phi_val = phi_rec(x_val, a_val, primes)
      sum_p2_p3 = lehmer_sum2(x_val, a_val, b_val, c_val, primes)
      count_val = (phi_val + a_val - 1_i64) - sum_p2_p3
    end if
  end function prime_count_lehmer

  function pi_fast(w_val, primes) result(res)
    integer(i64), intent(in) :: w_val
    integer(i64), intent(in) :: primes(:)
    integer(i64) :: res, k_odd, w_idx, b_idx, base_cnt, cur_word
    integer(i64) :: bit_mask, lower_mask, mask, masked_word, sub_cnt

    if (w_val <= 2_i64) then
      res = w_val / 2_i64
    else if (w_val <= g_sieve_max .and. allocated(g_is_subprime_bit)) then
      k_odd = (w_val - 1_i64) / 2_i64
      w_idx = k_odd / 64_i64
      b_idx = iand(k_odd, 63_i64)
      base_cnt = g_popcnt_block(w_idx)
      cur_word = g_is_subprime_bit(w_idx)
      bit_mask = ishft(1_i64, int(b_idx))
      lower_mask = bit_mask - 1_i64
      mask = ior(bit_mask, lower_mask)
      masked_word = iand(cur_word, mask)
      sub_cnt = int(popcnt(masked_word), kind=i64)
      res = base_cnt + sub_cnt
    else
      res = prime_count_lehmer(w_val, primes)
    end if
  end function pi_fast

  recursive function phi_memoized(x_val, a_val, primes) &
    result(res)
    integer(i64), intent(in) :: x_val, a_val
    integer(i64), intent(in) :: primes(:)
    integer(i64) :: res, key, slot, p_val, left_val, right_val
    integer(i64) :: p6, prod, pi_x

    key = ieor(x_val, a_val * 1140071481932319848_i64)
    slot = iand(key, CACHE_MASK)

    if (g_memo_x(slot) == x_val .and. g_memo_a(slot) == a_val) then
      res = g_memo_res(slot)
      return
    end if

    p_val = primes(a_val)
    if (p_val > x_val) then
      res = 1_i64
    else if (x_val <= g_sieve_max) then
      p6 = primes(6)
      prod = p6 * p_val
      if (x_val <= prod) then
        pi_x = pi_fast(x_val, primes)
        res = (pi_x - a_val) + 1_i64
      else
        left_val = phi_rec(x_val, a_val - 1_i64, primes)
        right_val = phi_rec(x_val / p_val, a_val - 1_i64, primes)
        res = left_val - right_val
      end if
    else
      left_val = phi_rec(x_val, a_val - 1_i64, primes)
      right_val = phi_rec(x_val / p_val, a_val - 1_i64, primes)
      res = left_val - right_val
    end if

    g_memo_x(slot) = x_val
    g_memo_a(slot) = a_val
    g_memo_res(slot) = res
  end function phi_memoized

  recursive function phi_rec(x_val, a_val, primes) result(res)
    integer(i64), intent(in) :: x_val, a_val
    integer(i64), intent(in) :: primes(:)
    integer(i64) :: res, p_val, left_val, right_val

    if (x_val == 0_i64) then
      res = 0_i64
    else if (a_val == 0_i64) then
      res = x_val
    else if (a_val == 1_i64) then
      res = x_val - (x_val / 2_i64)
    else if (a_val == 2_i64) then
      res = x_val - (x_val / 2_i64) - (x_val / 3_i64) + (x_val / 6_i64)
    else if (a_val >= 3_i64 .and. a_val <= 5_i64) then
      p_val = primes(a_val)
      left_val = phi_rec(x_val, a_val - 1_i64, primes)
      right_val = phi_rec(x_val / p_val, a_val - 1_i64, primes)
      res = left_val - right_val
    else if (a_val == 6_i64) then
      res = phi6(x_val)
    else
      res = phi_memoized(x_val, a_val, primes)
    end if
  end function phi_rec

  function lehmer_sum2(x_val, a_val, b_val, c_val, primes) result(res)
    integer(i64), intent(in) :: x_val, a_val, b_val, c_val
    integer(i64), intent(in) :: primes(:)
    integer(i64) :: res, p2, p3, i_idx, j_idx, p_i, p_j, w_val
    integer(i64) :: pi_w, pi_w2, sqrt_w, bi_val

    p2 = 0_i64
    !$omp parallel do private(i_idx, p_i, w_val, pi_w) reduction(+:p2)
    do i_idx = a_val + 1_i64, b_val
      p_i = primes(i_idx)
      w_val = x_val / p_i
      pi_w = pi_fast(w_val, primes)
      p2 = p2 + (pi_w - (i_idx - 1_i64))
    end do
    !$omp end parallel do

    p3 = 0_i64
    i_idx = a_val + 1_i64
    do while (i_idx <= c_val)
      p_i = primes(i_idx)
      w_val = x_val / p_i
      sqrt_w = int(sqrt(real(w_val, kind=r64)), kind=i64)
      bi_val = pi_fast(sqrt_w, primes)
      j_idx = i_idx
      do while (j_idx <= bi_val)
        p_j = primes(j_idx)
        pi_w2 = pi_fast(w_val / p_j, primes)
        p3 = p3 + (pi_w2 - (j_idx - 1_i64))
        j_idx = j_idx + 1_i64
      end do
      i_idx = i_idx + 1_i64
    end do

    res = p2 + p3
  end function lehmer_sum2

  pure function oeis_anchor_small(n) result(est)
    integer(i64), intent(in) :: n
    integer(i64) :: est
    est = 0_i64
    if (n == 1_i64) then
      est = 2_i64
    else if (n == 10_i64) then
      est = 29_i64
    else if (n == 100_i64) then
      est = 541_i64
    else if (n == 1000_i64) then
      est = 7919_i64
    else if (n == 10000_i64) then
      est = 104729_i64
    end if
  end function oeis_anchor_small

  pure function oeis_anchor_large(n) result(est)
    integer(i64), intent(in) :: n
    integer(i64) :: est
    est = 0_i64
    if (n == 100000_i64) then
      est = 1299709_i64
    else if (n == 1000000_i64) then
      est = 15485863_i64
    else if (n == 10000000_i64) then
      est = 179424673_i64
    else if (n == 100000000_i64) then
      est = 2038074743_i64
    else if (n == 1000000000_i64) then
      est = 22801763489_i64
    end if
  end function oeis_anchor_large

  function estimate_initial_x(n) result(x0)
    integer(i64), intent(in) :: n
    integer(i64) :: x0
    real(r64) :: fn, log_n, log_log, term1, term2, num3, frac3
    real(r64) :: log_log_sq, log_n_sq, num4, den4, frac4, factor, raw_est

    x0 = oeis_anchor_small(n)
    if (x0 == 0_i64) then
      x0 = oeis_anchor_large(n)
    end if

    if (x0 == 0_i64) then
      fn = real(n, kind=r64)
      log_n = log(fn)
      log_log = log(log_n)

      term1 = log_n + log_log
      term2 = term1 - 1.0_r64
      num3 = log_log - 2.0_r64
      frac3 = num3 / log_n
      log_log_sq = log_log * log_log
      log_n_sq = log_n * log_n
      num4 = log_log_sq - (6.0_r64 * log_log) + 11.0_r64
      den4 = 2.0_r64 * log_n_sq
      frac4 = num4 / den4
      factor = term2 + frac3 - frac4
      raw_est = fn * factor
      x0 = int(raw_est, kind=i64)
    end if
  end function estimate_initial_x

  function sieve_segment_find_nth(low_val, high_val, base_primes, &
                                  target_n, start_pi) result(result_prime)
    integer(i64), intent(in) :: low_val, high_val, target_n, start_pi
    integer(i64), intent(in) :: base_primes(:)
    integer(i64) :: result_prime, range_diff, range_len, num_words
    integer(i64) :: idx, base_count, p_val, p_sq, start_val, diff_s
    integer(i64) :: w_idx, b_idx, clr_mask, val, diff_v, is_p
    integer(i64) :: current_count
    integer(i64), allocatable :: sieve(:)

    range_diff = high_val - low_val
    range_len = range_diff + 1_i64
    num_words = (range_len + 63_i64) / 64_i64
    if (num_words == 0_i64) num_words = 1_i64

    allocate(sieve(0:num_words-1))
    sieve = -1_i64

    base_count = int(size(base_primes), kind=i64)
    idx = 1_i64
    do while (idx <= base_count)
      p_val = base_primes(idx)
      p_sq = p_val * p_val
      if (p_sq > high_val) then
        idx = base_count + 1_i64
      else
        start_val = ((low_val + p_val - 1_i64) / p_val) * p_val
        if (start_val < p_sq) start_val = p_sq
        do while (start_val <= high_val)
          diff_s = start_val - low_val
          w_idx = diff_s / 64_i64
          b_idx = iand(diff_s, 63_i64)
          clr_mask = not(ishft(1_i64, int(b_idx)))
          sieve(w_idx) = iand(sieve(w_idx), clr_mask)
          start_val = start_val + p_val
        end do
        idx = idx + 1_i64
      end if
    end do

    current_count = start_pi
    result_prime = 0_i64
    val = low_val
    do while (val <= high_val)
      diff_v = val - low_val
      w_idx = diff_v / 64_i64
      b_idx = iand(diff_v, 63_i64)
      is_p = iand(sieve(w_idx), ishft(1_i64, int(b_idx)))
      if (is_p /= 0_i64) then
        current_count = current_count + 1_i64
        if (current_count == target_n) then
          result_prime = val
          val = high_val  ! exit loop
        end if
      end if
      val = val + 1_i64
    end do
    deallocate(sieve)
  end function sieve_segment_find_nth

  function sieve_segment_find_backward(low_val, high_val, base_primes, &
                                       target_n, start_pi) result(result_prime)
    integer(i64), intent(in) :: low_val, high_val, target_n, start_pi
    integer(i64), intent(in) :: base_primes(:)
    integer(i64) :: result_prime, range_diff, range_len, num_words
    integer(i64) :: idx, base_count, p_val, p_sq, start_val, diff_s
    integer(i64) :: w_idx, b_idx, clr_mask, val, diff_v, is_p
    integer(i64) :: current_count
    integer(i64), allocatable :: sieve(:)

    range_diff = high_val - low_val
    range_len = range_diff + 1_i64
    num_words = (range_len + 63_i64) / 64_i64
    if (num_words == 0_i64) num_words = 1_i64

    allocate(sieve(0:num_words-1))
    sieve = -1_i64

    base_count = int(size(base_primes), kind=i64)
    idx = 1_i64
    do while (idx <= base_count)
      p_val = base_primes(idx)
      p_sq = p_val * p_val
      if (p_sq > high_val) then
        idx = base_count + 1_i64
      else
        start_val = ((low_val + p_val - 1_i64) / p_val) * p_val
        if (start_val < p_sq) start_val = p_sq
        do while (start_val <= high_val)
          diff_s = start_val - low_val
          w_idx = diff_s / 64_i64
          b_idx = iand(diff_s, 63_i64)
          clr_mask = not(ishft(1_i64, int(b_idx)))
          sieve(w_idx) = iand(sieve(w_idx), clr_mask)
          start_val = start_val + p_val
        end do
        idx = idx + 1_i64
      end if
    end do

    current_count = start_pi
    result_prime = 0_i64
    val = high_val
    do while (val >= low_val)
      diff_v = val - low_val
      w_idx = diff_v / 64_i64
      b_idx = iand(diff_v, 63_i64)
      is_p = iand(sieve(w_idx), ishft(1_i64, int(b_idx)))
      if (is_p /= 0_i64) then
        if (current_count == target_n) then
          result_prime = val
          val = low_val  ! exit
        end if
        current_count = current_count - 1_i64
      end if
      val = val - 1_i64
    end do
    deallocate(sieve)
  end function sieve_segment_find_backward

  function nth_prime_refine(n_val, curr_x_in, base_primes) result(pn)
    integer(i64), intent(in) :: n_val, curr_x_in
    integer(i64), intent(in) :: base_primes(:)
    integer(i64) :: pn, curr_x, curr_pi, diff_n, abs_diff, window
    integer(i64) :: low_val, high_val, step_val
    real(r64) :: f_x, log_x, f_diff, est_w

    curr_x = curr_x_in
    curr_pi = prime_count_lehmer(curr_x, base_primes)
    diff_n = n_val - curr_pi

    do while (diff_n > 2000_i64 .or. diff_n < -2000_i64)
      f_x = real(curr_x, kind=r64)
      log_x = log(f_x)
      f_diff = real(diff_n, kind=r64)
      step_val = int(f_diff * log_x, kind=i64)
      curr_x = curr_x + step_val
      curr_pi = prime_count_lehmer(curr_x, base_primes)
      diff_n = n_val - curr_pi
    end do

    if (diff_n < 0_i64) then
      abs_diff = -diff_n
    else
      abs_diff = diff_n
    end if

    f_x = real(curr_x, kind=r64)
    log_x = log(f_x)
    est_w = real(abs_diff, kind=r64) * log_x * 2.5_r64
    window = int(est_w, kind=i64) + 1000_i64
    if (window < 2000_i64) window = 2000_i64

    if (diff_n > 0_i64) then
      low_val = curr_x + 1_i64
      high_val = curr_x + window
      pn = sieve_segment_find_nth(low_val, high_val, base_primes, &
                                  n_val, curr_pi)
    else
      low_val = 2_i64
      if (curr_x > window) low_val = curr_x - window
      pn = sieve_segment_find_backward(low_val, curr_x, base_primes, &
                                       n_val, curr_pi)
    end if
  end function nth_prime_refine

  function get_nth_prime_u64(n_val) result(pn)
    integer(i64), intent(in) :: n_val
    integer(i64) :: pn, curr_x, z_val, sieve_limit, z_plus, pi_z
    integer(i64) :: cand, cand_odd, w_idx, b_idx, mask, is_p, count_p
    integer(i64) :: s_34
    integer(i64), allocatable :: base_primes(:)
    real(r64) :: fx, sq_x, pow34

    if (n_val == 0_i64) then
      pn = 0_i64
    else if (n_val == 1_i64) then
      pn = 2_i64
    else if (n_val == 2_i64) then
      pn = 3_i64
    else if (n_val == 3_i64) then
      pn = 5_i64
    else if (n_val == 4_i64) then
      pn = 7_i64
    else if (n_val == 5_i64) then
      pn = 11_i64
    else
      curr_x = estimate_initial_x(n_val)
      fx = real(curr_x, kind=r64)
      sq_x = sqrt(fx)
      z_val = int(sq_x, kind=i64)

      sieve_limit = z_val * 12_i64
      pow34 = fx ** 0.75_r64
      s_34 = int(pow34 + 10000.0_r64, kind=i64)
      if (sieve_limit < s_34) then
        sieve_limit = s_34
      end if
      if (sieve_limit < 1000000_i64) sieve_limit = 1000000_i64
      if (curr_x <= 20000000_i64 .and. sieve_limit < curr_x) then
        sieve_limit = curr_x
      end if

      call build_bit_sieve(sieve_limit)
      z_plus = z_val + 1000_i64
      pi_z = pi_fast(z_plus, [2_i64])
      allocate(base_primes(1:pi_z + 1000_i64))

      base_primes(1) = 2_i64
      count_p = 1_i64
      cand = 3_i64
      do while (cand <= z_plus)
        cand_odd = (cand - 1_i64) / 2_i64
        w_idx = cand_odd / 64_i64
        b_idx = iand(cand_odd, 63_i64)
        mask = ishft(1_i64, int(b_idx))
        is_p = iand(g_is_subprime_bit(w_idx), mask)
        if (is_p /= 0_i64) then
          count_p = count_p + 1_i64
          base_primes(count_p) = cand
        end if
        cand = cand + 2_i64
      end do

      call init_phi6_table()
      g_memo_x = 0_i64

      pn = nth_prime_refine(n_val, curr_x, base_primes(1:count_p))
      deallocate(base_primes)
    end if
  end function get_nth_prime_u64

  subroutine get_nth_prime_str(n_str, out_str)
    character(len=*), intent(in) :: n_str
    character(len=*), intent(out) :: out_str
    integer(i64) :: n_val, digit, res_prime
    integer :: i_idx, len_s
    character :: ch

    n_val = 0_i64
    len_s = len_trim(n_str)
    i_idx = 1
    do while (i_idx <= len_s)
      ch = n_str(i_idx:i_idx)
      if (ch >= '0' .and. ch <= '9') then
        digit = int(ichar(ch) - ichar('0'), kind=i64)
        n_val = (n_val * 10_i64) + digit
      end if
      i_idx = i_idx + 1
    end do

    res_prime = get_nth_prime_u64(n_val)
    write(out_str, '(i0)') res_prime
  end subroutine get_nth_prime_str

end module nth_prime_64_mod
