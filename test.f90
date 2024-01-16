program test
  use stdlib_kinds, functions
  implicit none 

  real(dp), allocatable :: r(:), f(:)
  integer :: i,j,n

  n = 50000
  allocate(r(n), f(n))

  call random_seed(put=666)
  call random_number(r)

  for i=1, n
    f = krc_pot(r,


