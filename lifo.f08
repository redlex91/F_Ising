module lifo

  implicit none
  type node
     integer :: data
     type( node ), pointer :: nxtPtr
  end type node

  type nodePtr
     type( node ), pointer :: Ptr
  end type nodePtr

contains
  subroutine push( stackP, info )
    ! Purpose: this sbrt creates a new node at the top of the stack.
    implicit none
    integer, intent( in ) :: info
    type( nodePtr ), intent ( inout ) :: stackP 
    type( nodePtr ) :: newP
    integer :: ok

    allocate( newP%Ptr, stat = ok )
    if( ok == 0 ) then ! allocation's been successful
       newP%Ptr%data = info
       newP%Ptr%nxtPtr => stackP%Ptr
       stackP%Ptr => newP%Ptr
    else
       print *, "Allocation error: no memory available!"
       call abort
    end if
  end subroutine push

  integer function pop( stackP )
    ! Purpose: this function deletes the last node at the top of the stack, frees the corresponding memory and return the valued in it.
    implicit none
    type( nodePtr ), intent( inout ) :: stackP
    type( nodePtr ) :: tmpP

    pop = stackP%Ptr%data
    tmpP%Ptr => stackP%Ptr
    stackP%Ptr => stackP%Ptr%nxtPtr
    deallocate( tmpP%Ptr )

  end function pop
  
end module lifo
