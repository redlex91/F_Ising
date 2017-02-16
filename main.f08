! Modello di Ising con algoritmo di Swendsen-Wang
! Reticolo quadrato LxL, assumo passo reticolare a=1, pertanto il
! numero di spin è N = L^2
!==============================================================================
!      START OF PROGRAM === START OF PROGRAM === START OF PROGRAM
!==============================================================================

program IsingModel
  ! dichiarazione variabili  
  use glob_data
  use utilities
  use lifo

  implicit none
  integer, dimension( 0 : 3 , 0 : N-1 ) :: NN_table 
  ! tabella dei primi vicini (NN) per ogni spin
  integer, dimension( 0 : N_BOND - 1, 1 : 2 ) :: BOND_table
  ! per ogni bond definisce gli spin legati, i bond orizz.
  ! sono definiti pari, quelli verticali sono dispari: se j
  ! è il numero di spin attuale, 2*j è il bond a dx, 2*j+1 è
  ! il bond in alto
  integer, dimension( 0 : 3, 0 : N-1 ) :: NN_BOND_table

  integer, dimension( 0 : N-1 ) :: spin
  ! vettore di configurazioni degli spin: for each point in the lattice it defines if its spin is up or down

  logical, dimension( 0 : N_BOND - 1 ) :: bond_stat ! if bond_stat is set to .true. the bond is active, otherwise it is inactive
  logical, dimension( 0 : N - 1 ) :: unexplored ! .false. if already visited
  logical :: flip, done, flag = .false. ! if flip is .true. we flip the cluster, flag is .false. for the first measure

  type( nodePtr ) :: clustP

  real( prec ) :: r
  ! numero random per inizializzazione
  ! degli spin

  real( prec ) :: ham, zeta = 0, ave_en = 0, mag = 0._prec, ave_mag = 0._prec, err, mean, stuff
  integer :: s, b, iter, curr_s, nn_s ! indice di spin: range 1-N, bond index: range 0:N_BOND-1

  integer :: i, j, var, bin ! indici per cicli

  ! executable statements

  open( unit = 4, file = 'energyvsbeta.dat', status = 'replace' )
  var = 120
  do   
     beta =  1._prec / ( real( var, prec ) / 100._prec )
     select case( var )
     case( 120 : 200 )
        var = var + 5
     case( 201 : 249 )
        var = var + 1
     case( 250 : 350 )
        var = var + 5
     case default
        exit
     end select

     p_bond = 1._prec - exp( -2._prec * beta )
     flag = .false.
     rewind 2

     call init_rand_seed( )
     ! inizializzo il seme del generatore di numeri casuali

     ! 0. GENERO LE TABELLE CHE DEFINISCONO LA TOPOLOGIA DEL
     ! SISTEMA

     call gen_table( NN_table, BOND_table, NN_BOND_table )

     open( 3, file = 'grid', status = 'replace', action = &
          & 'write' )
     open( 2, file = 'energy.dat', status = 'replace',  position = 'rewind ' )

     ! 1. INIZIALIZZO IL SISTEMA CON SPIN CASUALI
     do  s = 0, N-1
        call random_number( r )
        spin( s ) = (-1)**( nint( r ) ) 
        ! write( 3 , '(i5,a2)' , advance = 'no' ) spin(s), ' '
        ! if ( mod( s+1, L ) .eq. 0 ) then 
        !    write ( 3 , * ) ' '
        ! end if
     end do

     print *, 'Start evolution.'
     ! 2. START OF THE ALGORYTHM
     evolution: do iter = 1, MAX_ITER
        ! compute hamiltonian
        ham = 0._prec; mag = 0._prec
        do i = 0, N - 1
           do j = 0, 3
              if( NN_table( j, i ) > i ) ham = ham + ( - J_CONST * spin( NN_table( j, i ) ) * spin( i ) )
           end do
        end do
        ! compute magnetisation
        mag = sum( spin( 0 : N-1 ) )

        zeta = zeta + exp( - beta * ham )
        ave_en = ave_en + ( 1._prec / real( L, prec )**2._prec ) * ham * exp( - beta * ham )
        ave_mag = ave_mag + ( 1._prec / real( L, prec )**2._prec ) * abs( mag ) * exp( -beta * ham )
        if( mod( iter, MEAS_STEP ) == 0 ) then ! take measures
           if( flag ) then
              ave_en = ave_en / zeta
              ave_mag = ave_mag / zeta
              !           print *, ham, zeta, ave_en, ave_mag
              write( unit = 2, fmt = '(f9.3)' ) ave_en
              ave_en = 0._prec
              ave_mag = 0._prec
              zeta = 0._prec
           else ! discard first measure
              flag = .true.
              ave_en = 0._prec
              ave_mag = 0._prec
              zeta = 0._prec
           end if
        end if


        ! decide bonds' states
        bond_state: do b = 0, N_BOND - 1
           if( spin( BOND_table( b, 1 ) ) == spin( BOND_table( b, 2 ) ) ) then
              if( rand() <= p_bond ) then
                 ! activate bond
                 bond_stat( b ) = .true.
              else
                 bond_stat( b ) = .false.
              end if
           else
              bond_stat( b ) = .false.
           end if
        end do bond_state

        forall( i = 0 : N - 1 )
           unexplored( i ) = .true.
        end forall
        s = 0

        done = .false.
        flip_clusters: do
           if( done ) exit
           done = .true.
           ! the first element in the stack must ppoint to null
           nullify( clustP%Ptr )

           ! we evaluate the probability for flipping the cluster
           if( rand( ) <= .5_prec ) then
              flip = .true.
           else
              flip = .false.
           end if
           ! we mark the spin as visited
           unexplored( s ) = .false.
           if( flip ) spin( s ) = - spin( s ) 
           ! start new cluster from spin s
           call push( clustP, s )

           do
              if( .not.associated( clustP%Ptr ) ) exit
              curr_s = pop( clustP )
              ! check for active bonds among the current spin and its nn
              do i = 0, 3
                 nn_s = NN_table( i, curr_s )
                 if( bond_stat( NN_BOND_table( i, curr_s ) ) .and. unexplored( nn_s ) ) then
                    unexplored( nn_s ) = .false.
                    if( flip ) spin( nn_s ) = - spin( nn_s )
                    call push( clustP, nn_s )
                 end if
              end do
           end do

           ! look for the next spin which has not been visited yet
           do i = s + 1, N - 1
              if( unexplored( i ) ) then
                 s = i
                 done = .false. ! there is at least another spin to visit
                 exit
              end if
           end do
        end do flip_clusters
     end do evolution

     do s = 0, N - 1
        write( 3 , '(i5,a2)' , advance = 'no' ) spin(s), ' '
        if ( mod( s+1, L ) .eq. 0 ) then 
           write ( 3 , * ) ' '
        end if
     end do

     rewind 2
     bin = 0; j = 0
     do
        read( unit = 2, fmt = '(f9.3)', iostat = j ) stuff
        !     print *, j, stuff
        if( j/= 0 ) exit
        bin = bin + 1
     end do

     rewind 2
     mean = 0._prec; err = 0._prec
     do i = 1, bin
        read( unit = 2, fmt = '(f9.3)' ) stuff
        mean = mean + stuff
        err = err + stuff**2._prec
     end do
     mean = mean / real( bin, prec )
     err = err / real( bin, prec )
     err = sqrt( ( err - mean**2._prec ) / real( bin - 1, prec ) )
     write( unit = 4, fmt = '(3f9.3)' ) 1._prec/beta, mean, err
     print *, 1._prec/beta, mean, err
     !  print *, bin

     ! ! questa parte è da eliminare - inizio (1)
     ! open( 2, file = 'tables', status = 'replace', action =&
     !      & 'write', position = 'append' )
     ! ! apro il file in cui memeorizzo la tebella dei primi
     ! ! vicini

     ! print *, L,N,beta


     ! do i = 0, 3
     !    do j = 0, N-1
     !       write( 2, '(i5, a2)', advance = 'no' ), &
     !            & NN_table( i,j ), ' '
     !    end do
     !    write( 2, * ), '\\'
     ! end do

     ! do i = 0, N_BOND-1
     !    do j = 1, 2
     !       write( 2, '(i5, a2)', advance = 'no' ), &
     !            & BOND_table( i, j ), ' '
     !    end do
     !    write( 2, * ), '\\'
     ! end do

     ! do i = 0, 3
     !    do j = 0, N-1
     !       write( 2, '(i5, a2)', advance = 'no' ), &
     !            & NN_BOND_table( i, j ), ' '
     !    end do
     !    write( 2, * ), '//'
     ! end do

     ! close( 2 ) ! chiudo il file della tabella dei primi vicini
     ! ! questa parte è da eliminare - fine (1)
     close( 3 )
     close( 2, status = 'delete' )
  end do
  close( 4 )
end program IsingModel ! FINE DEL PROGRAMMA

