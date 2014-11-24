PROGRAM Squaring
  
!Physics is like sex: sure, it may give some practical results,
!but that's not why we do it. --Richard Feynman.

  use math
  USE init
  USE periodicbc

  IMPLICIT NONE
  
  INTEGER(kind = isoluku)                :: blokkeja
  INTEGER(kind = isoluku)                :: NMax
  INTEGER, PARAMETER                :: kanava = 6
  INTEGER                           :: allocstat, i, j, blk
  INTEGER                           :: siemen1
  CHARACTER(LEN=8)                  :: pvm
  CHARACTER(LEN=10)                 :: aloitusaika
  CHARACTER(LEN=30)                 :: teksti
  CHARACTER(LEN=10)                 :: n1,n2
    
  CALL read_INPUT(blokkeja,NMax)
  
  beta = tau_in*real(hiuk(1)%Trotter,kind=rk)
  Temperature = 1.0_rk/beta/kB 
  
  CALL alustamuuttujat()

  write(*,*) ' '
  WRITE(*,*) "--------------------------------------------" 
  write(*,*) 'Begin squaring '
  WRITE(*,*) "--------------------------------------------"
  write(*,*) ' '

  call alustaparipotentiaalit()

  write(*,*) ' '
  WRITE(*,*) "--------------------------------------------" 
  write(*,*) 'Done with the pair potentials '
  WRITE(*,*) "--------------------------------------------"
  write(*,*) ' '
  STOP

END PROGRAM Squaring

