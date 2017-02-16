module glob_data
  use prec_def

  implicit none
  public :: p_bond, beta
  integer, parameter :: L =  250, MAX_ITER = 5000, MEAS_STEP = 500 ! lunghezza del reticolo
  integer, parameter :: N = L**2 ! numero di spin
  integer, parameter :: N_BOND = 2*N 
  ! numero di legami = dimensione del sistema per numero di
  ! particelle

  real( prec ), parameter :: J_CONST = 1._prec ! exchange constant
  real( prec ) :: p_bond, beta
end module glob_data
