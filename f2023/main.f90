program main_program
  use, intrinsic :: iso_fortran_env, only : input_unit, output_unit
  use nth_prime_64_mod, only : get_nth_prime_u64
  use types_mod, only : i64, r128
  implicit none

  integer :: arg_count, ios
  character(len=256) :: arg_str
  integer(i64) :: n_val, result_prime

  arg_count = command_argument_count()

  if (arg_count > 0) then
    call get_command_argument(1, arg_str)
    call parse_input_number(arg_str, n_val, ios)
  else
    read(input_unit, '(a)', iostat=ios) arg_str
    if (ios == 0) then
      call parse_input_number(arg_str, n_val, ios)
    end if
  end if

  if (ios == 0) then
    if (n_val > 0_i64) then
      result_prime = get_nth_prime_u64(n_val)
      write(output_unit, "(i0)") result_prime
    end if
  end if

contains

  subroutine clean_input_string(raw_str, clean_str, len_clean)
    character(len=*), intent(in) :: raw_str
    character(len=256), intent(out) :: clean_str
    integer, intent(out) :: len_clean
    integer :: i, rlen
    character :: ch

    clean_str = ""
    len_clean = 0
    rlen = len(trim(raw_str))
    i = 1
    do while (i <= rlen)
      ch = raw_str(i:i)
      if (ch /= ',' .and. ch /= '_' .and. ch /= ' ') then
        len_clean = len_clean + 1
        clean_str(len_clean:len_clean) = ch
      end if
      i = i + 1
    end do
  end subroutine clean_input_string

  subroutine parse_ascii_digits(str, len_s, n_val, success)
    character(len=*), intent(in) :: str
    integer, intent(in) :: len_s
    integer(i64), intent(out) :: n_val
    logical, intent(out) :: success
    integer :: i
    character :: ch
    integer(i64) :: digit_val

    success = .true.
    n_val = 0_i64
    i = 1
    do while (i <= len_s)
      ch = str(i:i)
      if (ch >= '0' .and. ch <= '9') then
        digit_val = int(ichar(ch) - ichar('0'), i64)
        n_val = n_val * 10_i64
        n_val = n_val + digit_val
      else
        success = .false.
        i = len_s + 1
      end if
      if (success) then
        i = i + 1
      end if
    end do
  end subroutine parse_ascii_digits

  subroutine parse_power_expr(str, len_s, n_val, success)
    character(len=*), intent(in) :: str
    integer, intent(in) :: len_s
    integer(i64), intent(out) :: n_val
    logical, intent(out) :: success
    integer :: p_idx, ios1, ios2
    integer(i64) :: base_val, exp_val

    success = .false.
    n_val = 0_i64

    p_idx = index(str(1:len_s), '**')
    if (p_idx > 0) then
      read(str(1:p_idx-1), *, iostat=ios1) base_val
      read(str(p_idx+2:len_s), *, iostat=ios2) exp_val
      if (ios1 == 0 .and. ios2 == 0) then
        n_val = base_val ** exp_val
        success = .true.
        return
      end if
    end if

    p_idx = index(str(1:len_s), '^')
    if (p_idx > 0) then
      read(str(1:p_idx-1), *, iostat=ios1) base_val
      read(str(p_idx+1:len_s), *, iostat=ios2) exp_val
      if (ios1 == 0 .and. ios2 == 0) then
        n_val = base_val ** exp_val
        success = .true.
        return
      end if
    end if
  end subroutine parse_power_expr

  subroutine parse_float_expr(str, len_s, n_val, success)
    character(len=*), intent(in) :: str
    integer, intent(in) :: len_s
    integer(i64), intent(out) :: n_val
    logical, intent(out) :: success
    integer :: ios
    real(r128) :: rval

    success = .false.
    n_val = 0_i64
    read(str(1:len_s), *, iostat=ios) rval
    if (ios == 0) then
      n_val = int(rval, i64)
      success = .true.
    end if
  end subroutine parse_float_expr

  subroutine parse_input_number(raw_str, n_val, ios)
    character(len=*), intent(in) :: raw_str
    integer(i64), intent(out) :: n_val
    integer, intent(out) :: ios

    character(len=256) :: clean_str
    integer :: len_clean
    logical :: success

    call clean_input_string(raw_str, clean_str, len_clean)
    ios = -1
    n_val = 0_i64

    if (len_clean == 0) return

    call parse_ascii_digits(clean_str, len_clean, n_val, success)
    if (success) then
      ios = 0
      return
    end if

    call parse_power_expr(clean_str, len_clean, n_val, success)
    if (success) then
      ios = 0
      return
    end if

    call parse_float_expr(clean_str, len_clean, n_val, success)
    if (success) then
      ios = 0
    end if
  end subroutine parse_input_number

end program main_program
