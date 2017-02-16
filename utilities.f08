module utilities
  use prec_def

contains
  
  subroutine gen_table( NN_table, BOND_table, NN_BOND_table )
    ! crea la tabella dei primi vicini per ogni spin del
    ! reticolo: riga 1 - NN a destra, gli altri in senso
    ! antiorario a partire da quello a destra

    use glob_data

    implicit none
    integer :: i ! indice dello spin
    integer :: x0, y0, x1, x2, y1, y2
    ! (x0,y0) posizione dello spin attuale, (x1,y0) NN a dx, (x2,y0)
    ! NN a sx, (x0,y1) NN sopra, (x0,y2) NN sotto

    integer, dimension( 0 : 3 ) :: nn 
    ! vettore contentene gli indici dei NN

    integer, dimension( 0 : 3 , 0 : N-1 ) :: NN_table 
    integer, dimension( 0 : N_BOND - 1, 1 : 2 ) :: BOND_table
    ! per ogni bond definisce gli spin legati, i bond orizz.
    ! sono definiti pari, quelli verticali sono dispari: se j
    ! è il numero di spin attuale, 2*j è il bond a dx, 2*j+1 è
    ! il bond in alto
    integer, dimension( 0 : 3, 0 : N-1 ) :: NN_BOND_table
    ! per ciascuno spin definisce con quali bond è collegato ai
    ! suoi primi vicini

    do i = 0, N-1

       x0 = mod( i, L )
       y0 = i / L

       x1 = mod( x0 + 1, L )
       x2 = mod( x0 - 1 + L, L )
       y1 = mod( y0 + 1, L )
       y2 = mod( y0 -1 + L, L )

       nn(0) = x1 + y0*L
       nn(1) = y1*L + x0
       nn(2) = x2 + y0*L
       nn(3) = y2*L + x0

       NN_table( : , i ) = nn( : )
       ! costruisco la tabella dei primi vicini

       BOND_table( 2*i, : ) = (/ i, nn(0) /)
       BOND_table( 2*i + 1 , : ) = (/ i, nn(1) /)
       ! costruisco la tabella dei bond

       NN_BOND_table( 0 , i ) = 2*i
       NN_BOND_table( 1 , i ) = 2*i + 1
       NN_BOND_table( 2 , nn(0) ) = 2*i
       NN_BOND_table( 3 , nn(1) ) = 2*i + 1

    end do
  end subroutine gen_table

  subroutine init_rand_seed
    ! inizializza il seme del generatore di numeri casuali

    integer :: clock, n, x
    integer, dimension( : ), allocatable :: s ! seed

    call random_seed( size = n )
    allocate( s( n ) )

    call system_clock( count = clock )

    s = clock + (/(i-1, i = 1, n)/)*23

    call random_seed( put = s )
    deallocate( s )

  end subroutine init_rand_seed

end module utilities
