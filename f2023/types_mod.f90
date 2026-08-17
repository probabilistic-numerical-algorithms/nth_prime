module types_mod
  use, intrinsic :: iso_fortran_env, only : int64, real64, real128
  implicit none
  private
  public :: i64, r64, r128

  integer, parameter :: i64 = int64
  integer, parameter :: r64 = real64
  integer, parameter :: r128 = real128
end module types_mod
