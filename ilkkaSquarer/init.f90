

MODULE init

  use math
  
  IMPLICIT NONE
  
  TYPE:: hi_tiedot
     
     INTEGER(kind = isoluku) :: Trotter
     INTEGER :: lkm ! number of similar particles
     INTEGER :: start_i ! start index for species
     INTEGER :: species ! number of the species
     INTEGER                 :: multi_L
     INTEGER                 :: ParticleKind !0->Boltzmannon, 1->Fermion
     CHARACTER( len=3 )      :: tyyppi
     CHARACTER( len=11 )     :: statistics !Boltzmannon or Fermion
     REAL(KIND=RK)               :: Z
     REAL(KIND=RK)               :: meff
     REAL(KIND=RK)               :: lambda
     REAL(KIND=RK)               :: spin
     REAL(KIND=RK), DIMENSION(:,:), ALLOCATABLE  :: vpaikka
     REAL(KIND=RK), DIMENSION(:,:), ALLOCATABLE  :: upaikka

     INTEGER(kind = isoluku), DIMENSION(:), ALLOCATABLE :: permutation
     INTEGER(kind = isoluku), DIMENSION(:), ALLOCATABLE :: permutationOld
    
     
     !SVA@CSC siirretytpaikat contains a list of beads 
     !which have moved during one MC trial step. Used to
     !selectively copy bead positions between upaikka nad vpaikka
     INTEGER(kind = isoluku), DIMENSION(:), ALLOCATABLE ::siirretytpaikat 
     !SVA@CSC  siirretytpaikatlkm contains the number bead indexes in siirretytpaikat 
     INTEGER:: siirretytpaikatlkm     
     
  END TYPE hi_tiedot
   
  TYPE( hi_tiedot ), DIMENSION(:), ALLOCATABLE,SAVE  :: hiuk


  character(len=11), dimension(3), private :: ParticleKind = &
       (/'Boltzmannon','Fermion    ','Boson      '/)

  
  TYPE:: exchangeMatrix
    
    INTEGER(kind = isoluku)            :: Trotter
    INTEGER            :: lkm ! number of ex particles
    INTEGER            :: nauha
    REAL(KIND=RK)               :: Z
    REAL(KIND=RK)               :: meff
    REAL(KIND=RK)               :: lambda
    REAL(KIND=RK)               :: spin
    REAL(KIND=RK), DIMENSION(:,:), ALLOCATABLE  :: Matrix
    
  END TYPE exchangeMatrix
  
  TYPE( exchangeMatrix ), DIMENSION(:), ALLOCATABLE,SAVE  :: exchange

  
  TYPE:: effpot_tiedot
     
     REAL(KIND=RK)    :: ZZ
     REAL(KIND=RK)    :: m_red
     REAL(KIND=RK)    :: lambda_red
     integer          :: modcusp 
     integer          :: potnro = 0     
     
  END TYPE effpot_tiedot
  
  TYPE( effpot_tiedot ), DIMENSION(:,:), ALLOCATABLE,SAVE  :: effpot_data

  integer, save :: NumPairPots
  integer, save :: pairpot_dim
  integer(kind = isoluku), save :: NumSquares

  logical, save :: GetBlockConf

  ! if DisplaceMoves=0.1 then about 10% 
  ! of the moves are displace moves
  ! and 90% are multilevel bisection moves
  real(kind = rk), save :: DisplaceMoves

  real(kind = rk), dimension(:), allocatable :: rex

  ! Long range parameters
  real(kind = rk), dimension(:), allocatable, private :: kvalues
  real(kind = rk), dimension(:), allocatable, private :: Vk
  real(kind = rk), save :: alpha
  
  real(kind = rk), dimension(:,:), allocatable, save :: karray
  real(kind = rk), save :: kcutoff
  real(kind = rk), save :: kcutoff_in
  integer, save :: karray_dim1  
  integer, save :: PPn_max
  logical, save :: LongRangeOn

  integer, save   :: NumImages
  CHARACTER(LEN=11), SAVE :: NumOfImages

  integer, save :: NumSpecies
  integer, save :: NumSpeciesEQ1

  logical, save :: ReadPairPotFile

  logical, save :: ElectronHoleExOn

  integer, save                  :: PCGridSize
  real(kind=rk), save            :: PCLowLimit
  real(kind=rk), save            :: PCUppLimit

  logical, save :: LevelUpdateOn
  real(kind=rk), save :: AccLow
  real(kind=rk), save :: AccHigh       

  integer, dimension(:), allocatable :: PathSign
  integer, dimension(:), allocatable :: PathSignOld


  real(kind=rk), save            :: beta,tau_in
  REAL(KIND=RK), SAVE            :: Temperature
  INTEGER, SAVE                  :: otosvali
  
  REAL(KIND=RK), PARAMETER       :: a0 = 0.529177208e-10_rk
  REAL(KIND=RK), PARAMETER       :: ev = 3.6749326014e-2_rk
  REAL(KIND=RK), PARAMETER       :: hviiva =  1.0_rk

  INTEGER, SAVE                :: dimensio
  INTEGER, DIMENSION(3), SAVE  :: hi_lkm
  INTEGER, SAVE                :: s_lkm
  REAL(KIND=RK), DIMENSION(2,2), SAVE   :: Uex


  LOGICAL, save                :: alusta_satunnais
  LOGICAL, SAVE                :: bisection_hyv
  LOGICAL, SAVE                :: exchange_on
  LOGICAL, SAVE                :: DisplaceMoveOn
  LOGICAL, SAVE                :: ThermalEstimator
  CHARACTER(LEN=30), SAVE      :: otsikko
  character(len=40), save      :: tiedAlku_e, tiedAlku_y





CONTAINS





!--------------------------------------------------------
!                         ALUSTA
!--------------------------------------------------------
! Tämä tekee kaikki alustamiseen tarvittavat. Siis 
!  -lukee tiedostoista vanhat arvot tai
!  -keksii tyhjästä uudet arvot

  SUBROUTINE alustamuuttujat
    use math
    IMPLICIT NONE
    
    INTEGER               :: i, j
    LOGICAL               :: olemassa
  
    CHARACTER(LEN=20)     :: TrotterChar,NParticlesChar
    CHARACTER(LEN=20)     :: nimi


    
    ! --------------------------
    ! Elektronit (or negative particles)
    ! --------------------------

    !elektroneilla sama Trotterin luku
    WRITE( TrotterChar, '(I10)'  ) hiuk(1)%Trotter 
    WRITE( NParticlesChar, '(I5)'  ) hi_lkm(1)

       
    nimi = 'M' // TRIM(ADJUSTL(TrotterChar)) // 'N' //&
         TRIM(ADJUSTL(NParticlesChar))
    
    INQUIRE ( FILE = TRIM(nimi) // '.electrons', EXIST=olemassa )

    tiedAlku_e = nimi

    IF ( olemassa ) THEN
       
       OPEN( 10, FILE=TRIM(nimi) // '.electrons' )
       
       DO i = 1, hi_lkm(1) 
          DO j = 1, hiuk(i)%Trotter
             
             READ( 10, * ) hiuk(i)%vpaikka(1,j), hiuk(i)%vpaikka(2,j),&
                  hiuk(i)%vpaikka(3,j)
             
          END DO
       END DO
       
       CLOSE( 10 )
       
    ELSE       
       
       !WRITE(*,*) "------------------------------------------------"
       !WRITE(*,*) "Olis syytä olla elektronien paikat jo tiedossa!!"
       !WRITE(*,*) "------------------------------------------------"
       
       DO i = 1, hi_lkm(1) 
          hiuk(i)%vpaikka = 0.0_rk
          hiuk(i)%upaikka = 0.0_rk
          hiuk(i)%siirretytpaikatlkm=0
       END DO
       
    END IF




    ! -------------------
    ! Ytimet (positive particles)
    ! -------------------

    ! ytimilla sama Trotterin luku
    if (hi_lkm(2)==0) then
       WRITE( TrotterChar, '(I10)'  ) 0 
    else
       WRITE( TrotterChar, '(I10)'  ) hiuk(hi_lkm(1)+1)%Trotter 
    end if
    WRITE( NParticlesChar, '(I5)'  ) hi_lkm(2)
        
   !WRITE( TrotterChar, '(I10)'  ) hiuk(hi_lkm(1)+1)%Trotter 
   WRITE( NParticlesChar, '(I5)'  ) hi_lkm(2)

       
   nimi = 'M' // TRIM(ADJUSTL(TrotterChar)) // 'N' //&
        & TRIM(ADJUSTL(NParticlesChar))
    
   INQUIRE ( FILE = TRIM(nimi) // '.nuclei', EXIST=olemassa )

   tiedAlku_y = nimi
    
   IF ( olemassa ) THEN
       
      OPEN( 10, FILE=TRIM(nimi) // '.nuclei' )
       
      DO i = hi_lkm(1)+1,hi_lkm(3) 
         DO j = 1, hiuk(i)%Trotter
            
            READ( 10, * ) hiuk(i)%vpaikka(1,j), hiuk(i)%vpaikka(2,j),&
                 hiuk(i)%vpaikka(3,j)
            
         END DO
      END DO

      CLOSE( 10 )
           
   ELSE

      if (hi_lkm(2) .ne. 0) then
         !WRITE(*,*) "--------------------------------------------"
         !WRITE(*,*) "Olis syytä olla ytimien paikat jo tiedossa!!"
         !WRITE(*,*) "--------------------------------------------"
         
         DO i = hi_lkm(1)+1, hi_lkm(3) 
            hiuk(i)%vpaikka = 0.0_rk
            hiuk(i)%upaikka = 0.0_rk
            hiuk(i)%siirretytpaikatlkm=0
         END DO
      end if
      
   END IF
        
  END SUBROUTINE alustamuuttujat
  

  !-----------------------------------------------
  ! Alustetaan satunnaislukugeneraattori
  !

  SUBROUTINE alustasatunnaisgeneraattori( siemen1 ) 
    
    IMPLICIT NONE
    
    INTEGER, INTENT(OUT)               :: siemen1
    INTEGER, DIMENSION(:), ALLOCATABLE :: siemen
    INTEGER, DIMENSION(8)              :: t
    INTEGER                            :: koko, tila


    !Tehdään vähän satunnaisuutta
    CALL RANDOM_SEED(size=koko)

    ALLOCATE(siemen(koko), stat=tila)
    IF ( tila > 0 ) STOP 'Tilanvaraus epäonnistui'

    IF ( alusta_satunnais ) THEN
       CALL DATE_AND_TIME(values = t)
       !siemen = 100*t(7) + t(8)/10 + 1
       siemen = 1000*t(7) + t(8) + 1

       IF (siemen(1) == 0) THEN
          siemen = 1768
       END IF
    ELSE
       siemen = 1768
    END IF

    CALL RANDOM_SEED( put=siemen )

    siemen1 = siemen(1) 

  END SUBROUTINE alustasatunnaisgeneraattori


  SUBROUTINE read_INPUT(blokkeja,NMax)
    use math
    use periodicbc
    IMPLICIT NONE
    
    INTEGER(kind = isoluku) :: blokkeja,NMax,Trotter
    integer                 :: i,nP1,nP_tot
    INTEGER               :: lkm, lkm2,spin, multi_L,ex_lkm
    REAL(KIND=RK)         :: Z, m
    real(kind = rk)       :: Lx,Ly,Lz, temp
    integer               :: tauORtemp,intParticleKind
    integer               :: intDisplaceMoveOn, ReadFile
    INTEGER               :: allocstat
    !CHARACTER( len=5 )    :: exchange_on_text
    LOGICAL               :: olemassa
    integer               :: ThermEstim,LRon,GetConf,Lupdateon,eh_ex
    CHARACTER( len=20 )   :: estimator,paircors,helptext,LongRange,Conf,Lupdate

    ThermalEstimator = .false.
    
    INQUIRE ( FILE = 'INPUT', EXIST=olemassa )
    
    IF ( olemassa ) THEN
       
       lkm2 = 0
       ex_lkm = 0
       dimensio = 3
       hi_lkm = 0

       OPEN( 10, FILE='INPUT')
       
       READ( 10, *) otsikko
       READ( 10, *) nP1, nP_tot

       NumSpecies = nP1

       hi_lkm(3) = nP_tot

       ! varataan muistia hiukkasille
       ALLOCATE( hiuk( nP_tot ), STAT = allocstat )
       IF ( allocstat /= 0 ) THEN
          WRITE(*,*) 'Virhe muistinvarauksessa'
          STOP
       END IF

       !luetaan hiukkasten datat
       NumSpeciesEQ1=0
       DO i = 1,nP1
          READ( 10, *) lkm, Z, m, spin, Trotter, multi_L,intParticleKind
          if (lkm==1) then
             NumSpeciesEQ1 = NumSpeciesEQ1 + 1
          end if
             
          CALL alustaHiukkaset ( lkm2, lkm, Z, m, &
               spin, Trotter, multi_L,intParticleKind,i)
          lkm2 = lkm2 + lkm
          IF (lkm>1 .and. intParticleKind==1) THEN
             ex_lkm = ex_lkm + 1
          END IF
       END DO
       
       
       READ( 10, *) s_lkm
       READ( 10, *) otosvali
       READ( 10, *) temp, tauORtemp
       READ( 10, *) blokkeja, NMax
       READ( 10, *) helptext,pairpot_dim, NumSquares, ReadFile
       READ( 10, *) helptext,Lx,Ly,Lz
       READ( 10, *) NumOfImages, NumImages
       READ( 10, *) helptext,intDisplaceMoveOn, DisplaceMoves
       READ( 10, *) estimator, ThermEstim
       READ( 10, *) paircors, PCLowLimit, PCUppLimit, PCGridSize
       READ( 10, *) LongRange, LRon, kcutoff_in
       READ( 10, *) Conf, GetConf
       READ( 10, *) Lupdate,Lupdateon,AccLow,AccHigh
       
       if (ReadFile==1) then
          ReadPairPotFile = .true.
       else
          ReadPairPotFile = .false.
       end if

       if (Lupdateon==1) then
          LevelUpdateOn = .true.
       else
          LevelUpdateOn = .false.
       end if       
       if (GetConf==1) then
          GetBlockConf = .true.
       else
          GetBlockConf = .false.
       end if
       if (LRon==1) then
          LongRangeOn=.true.
       else
          LongRangeOn=.false.
       end if
       
       if (ThermEstim==1) then
          ThermalEstimator = .true.
       else
          ThermalEstimator = .false.
       end if

       if (intDisplaceMoveOn==1) then
          DisplaceMoveOn = .true.
       else
          DisplaceMoveOn = .false.
       end if

       if (tauORtemp==0) then
          tau_in = temp
       else
          tau_in = 1.0_rk/kB/Temp/hiuk(1)%Trotter
       end if
       
       periodic_on = .true. ! always true for this version of the code
       call SetBox(Lx,Ly,Lz)       
       
       CLOSE( 10 )
       
       IF (ex_lkm > 0) THEN
          exchange_on = .TRUE.
          ALLOCATE( exchange( ex_lkm ), STAT = allocstat )
          IF ( allocstat /= 0 ) THEN
             WRITE(*,*) 'Virhe muistinvarauksessa'
             STOP
          END IF
          CALL alustaExchange 
       ELSE
          exchange_on = .FALSE.
       END IF
           
    ELSE
       WRITE(*,*) 'You need to have an INPUT file'
    END IF

  END SUBROUTINE read_INPUT

  SUBROUTINE alustaHiukkaset( lkm2, lkm, Z, m, spin, Trotter, multi_L,intParticleKind,SpeciesNum)
    use math
    IMPLICIT NONE
    INTEGER :: hi,SpeciesNum
    INTEGER :: lkm, lkm2, spin, multi_L,intParticleKind
    integer(kind = isoluku) :: Trotter
    REAL(KIND=RK)    :: Z, m
    INTEGER :: allocstat

    DO hi = lkm2 + 1, lkm2 + lkm

       hiuk(hi)%Trotter = Trotter
       hiuk(hi)%lkm     = lkm
       hiuk(hi)%start_i = lkm2+1
       hiuk(hi)%species = SpeciesNum
       hiuk(hi)%meff    = m
       hiuk(hi)%multi_L = multi_L
       hiuk(hi)%spin    = spin
       hiuk(hi)%Z       = Z
       hiuk(hi)%ParticleKind = intParticleKind
       hiuk(hi)%statistics = ParticleKind(intParticleKind+1)
       IF ( hiuk(hi)%Z < 0.) THEN
          hiuk(hi)%tyyppi  = 'neg'
          hi_lkm(1) = hi_lkm(1) + 1
       ELSE
          hiuk(hi)%tyyppi  = 'pos'
          hi_lkm(2) = hi_lkm(2) + 1
       END IF       
       
       if (hiuk(hi)%meff>1.0e10_rk .or. hiuk(hi)%meff<1.0e-10_rk) then
          hiuk(hi)%lambda = 0.0_rk
       else
          hiuk(hi)%lambda = 1.0_rk/hiuk(hi)%meff/2
       end if
       
       ALLOCATE(hiuk(hi)%vpaikka(dimensio,hiuk(hi)%Trotter),&
            STAT=allocstat)
       ALLOCATE(hiuk(hi)%upaikka(dimensio,hiuk(hi)%Trotter),&
            STAT=allocstat)

       hiuk(hi)%vpaikka = 0.0_rk
       hiuk(hi)%upaikka = 0.0_rk

       
       !ALLOCATE(hiuk(hi)%siirretytpaikat(2**hiuk(hi)%multi_L),&
       !    STAT=allocstat)
       ALLOCATE(hiuk(hi)%permutation(hiuk(hi)%Trotter),&
            STAT=allocstat)
       ALLOCATE(hiuk(hi)%permutationOld(hiuk(hi)%Trotter),&
            STAT=allocstat)
       hiuk(hi)%permutation = hi
       hiuk(hi)%permutationOld = hi

       ALLOCATE(hiuk(hi)%siirretytpaikat(hiuk(hi)%Trotter),&
            STAT=allocstat)
       hiuk(hi)%siirretytpaikat = 0
       hiuk(hi)%siirretytpaikatlkm = 0

     END DO
     
    IF ( allocstat /= 0 ) THEN
       WRITE(*,*) 'Virhe muistinvarauksessa'
       STOP
    END IF
    
  END SUBROUTINE alustaHiukkaset

  SUBROUTINE alustaExchange
    use math
    IMPLICIT NONE
    INTEGER :: i, j
    INTEGER :: allocstat

    j = 1
    i = 1
    DO WHILE (i<hi_lkm(3))
       IF (hiuk(i)%lkm > 1 .and. hiuk(i)%ParticleKind==1) THEN
          exchange(j)%Trotter = hiuk(i)%Trotter
          exchange(j)%meff    = hiuk(i)%meff
          exchange(j)%lambda  = hiuk(i)%lambda
          exchange(j)%lkm     = hiuk(i)%lkm
          exchange(j)%nauha   = i
          exchange(j)%Z       = hiuk(i)%Z
          
          ALLOCATE(exchange(j)%Matrix(hiuk(i)%lkm,hiuk(i)%lkm),&
               STAT=allocstat)
          
          exchange(j)%Matrix = 0.0_rk
          j = j + 1
          i = i + hiuk(i)%lkm
       ELSE
          i = i+1
       END IF
    END DO
    
    IF ( allocstat /= 0 ) THEN
       WRITE(*,*) 'Virhe muistinvarauksessa'
       STOP
    END IF
    
  END SUBROUTINE alustaExchange

  
  subroutine alustaparipotentiaalit
    use math
    use periodicbc
    implicit none
    integer :: i,j,m,n,jj,k, r_dim,kk,ind,ind2
    real(kind = rk), dimension(:),allocatable :: r
    real(kind = rk), dimension(3) :: Box
    logical :: loyty

    write(*,*) ' '

    allocate(effpot_data( hi_lkm(3),hi_lkm(3) ) )
    
    do i=1,hi_lkm(3)
       do j=1,hi_lkm(3)
          effpot_data(i,j)%ZZ = hiuk(i)%Z*hiuk(j)%Z
          effpot_data(i,j)%m_red = 1.0_rk/(1.0_rk/hiuk(i)%meff + 1.0_rk/hiuk(j)%meff )
          effpot_data(i,j)%lambda_red = hviiva**2/(2.0_rk*effpot_data(i,j)%m_red)
          effpot_data(i,j)%potnro = 0
          !if (hiuk(i)%species==hiuk(j)%species .and. hiuk(i)%ParticleKind==1) then
          !   effpot_data(i,j)%modcusp = 1
          !else
             effpot_data(i,j)%modcusp = 0
          !end if

       end do
    end do

    jj = 1
    do i=1,hi_lkm(3)-1
       do j=i+1,hi_lkm(3)

          loyty = .false.

          if (i==1 .and. j==2) then
             effpot_data(i,j)%potnro = jj
          else
             do n=1,hi_lkm(3)-1
                do m=n+1,hi_lkm(3)
                   if (n .ne. i .or. m .ne. j) then
                      if (effpot_data(i,j)%ZZ==effpot_data(n,m)%ZZ .and.&
                           effpot_data(i,j)%m_red==effpot_data(n,m)%m_red .and.&
                           !effpot_data(i,j)%modcusp==effpot_data(n,m)%modcusp .and.&
                           effpot_data(n,m)%potnro .ne. 0) then                         
                         effpot_data(i,j)%potnro = effpot_data(n,m)%potnro
                         loyty = .true.
                      end if
                   end if
                end do
             end do
             if (.not. loyty) then
                jj = jj + 1
                effpot_data(i,j)%potnro = jj
             end if                          
          end if

          effpot_data(j,i)%potnro = effpot_data(i,j)%potnro
       end do
    end do
    
    NumPairPots = jj
    
    r_dim = pairpot_dim
    allocate(r(r_dim))
    

    call alustaCSpline2d2(NumPairPots,r_dim,r_dim)
    if (ThermalEstimator) then
       call alusta_dU_dbeta(NumPairPots,r_dim,r_dim)    
    end if
    if (LongRangeOn) then
       allocate(SRcorrection(NumPairPots))
    end if
    
    if (tau_in/2**Numsquares>1.0e-11_rk) then
       open(123,file='NumSquaresModified.dat',status='replace')
       write(123,*) 'Modifying NumSquares = ', NumSquares
       do while (tau_in/2**Numsquares>1.0e-11_rk)
          NumSquares = NumSquares + 1
       end do
       write(123,*) 'Modified NumSquares = ', NumSquares
       write(123,*) 'giving tau_start    = ', tau_in/2**Numsquares
       close(123)
    end if
       
    
    !r = gridIK2(1.0e-4_rk,r_dim)
    
    do i=1,NumPairPots       
       do kk=1,hi_lkm(3)-1
          do jj=kk+1,hi_lkm(3)
             if (effpot_data(kk,jj)%potnro==i) then
                k=kk
                j=jj
             end if
          end do
       end do

       r = gridIK2(sqrt(2.0_rk*effpot_data(k,j)%lambda_red*tau_in/2**NumSquares),r_dim)
       if (1==1) then
          ind2=1
          do ind=1,r_dim
             if (r(ind)<abs(effpot_data(k,j)%ZZ)) then
                ind2=ind
             end if
          end do
          if (r_dim-ind2>2) then
             r(ind2+1:r_dim) = gridExtend(r(ind2+1),r(ind2+2)-r(ind2+1),r_dim-ind2)
          end if
       end if

       if (ThermalEstimator) then
          call GeneratePairPotWithBetaDeriv(effpot_data(k,j)%lambda_red,effpot_data(k,j)%ZZ,r,r_dim,i,effpot_data(k,j)%modcusp)    
       else
          call GeneratePairPot(effpot_data(k,j)%lambda_red,effpot_data(k,j)%ZZ,r,r_dim,i,effpot_data(k,j)%modcusp)    
       end if
    end do
    
    deallocate(r)

  end subroutine alustaparipotentiaalit


  subroutine GeneratePairPot(lambda,ZZ,r,r_dim,PairPotNum,modcusp)
    use math    
    use matrixsquaring
    use periodicbc
    implicit none
    real(kind = rk), intent(in) :: lambda,ZZ
    integer, intent(in) :: r_dim,PairPotNum,modcusp
    real(kind = rk), dimension(r_dim), intent(in) :: r
    real(kind = rk) :: dr,tau,kerroin,u0approx,z,cuspcor
    integer(kind = isoluku) :: M
    integer :: dr_dim,i,j,iovar
    real(kind = rk), dimension(:), allocatable :: r2
    real(kind = rk), dimension(:,:), allocatable :: u0
    real(kind = rk), dimension(r_dim,r_dim) :: u    
    real(kind = rk), dimension(3) :: Box
    character(len=20) :: i1,help_text
    logical :: generate,olemassa
    
    generate = .true.

    !if here then LongRangeOn needs to be false, sorry
    if (LongRangeOn) then
       write(*,*) 'Long Range set to false due to Estimator!'
       write(*,*) ' '
    end if
    LongRangeOn=.false.
    
    !if (lambda<0.01_rk) then 
    !   write(*,'((A),(I2),(A3),(I2))') 'Generating Pair Potential ', PairPotNum, ' / ', NumPairPots
    !   write(*,*) 'Bare Coulomb Potential: lambda = ', lambda
    !   write(*,*) ' '
    !   u=0.0_rk
    !   call initCSpline2d2(r,r,u,r_dim,r_dim,PairPotNum)
    !   return
    !end if
    
    if (ReadPairPotFile) then
       
       write( i1, '(I6)'  ) PairPotNum
       inquire(file='us.' // trim(adjustl(i1)) // '.txt', EXIST=olemassa )
       if (olemassa) then
          write(*,'((A),(I2),(A3),(I2))') 'Reading Pair Potential ', PairPotNum, ' / ', NumPairPots
          write(*,*) ' '

          allocate(r2(r_dim))
          r2=0.0_rk
          open(123,file='us.' // trim(adjustl(i1)) // '.txt')
          do i=1,r_dim
             do j=1,r_dim
                read(123,'(3(ES20.10))') r2(i),r2(j),u(i,j)
             end do
             read(123,'(A)',iostat=iovar) help_text
          end do
          close(123)
          
          call initCSpline2d2(r2,r2,u,r_dim,r_dim,PairPotNum)
          generate = .false.
       else
          write(*,*) 'Cannot find pair potential... thus, generating'
          write(*,*) ' '
          generate = .true.
       end if
    end if

    if (generate) then
       M = 2**NumSquares
       tau = tau_in/M    
       
       lambdaMS = lambda
       ZZMS = ZZ
       
       ! addition to the grid wanted
       dr = 20.0_rk*sqrt(2.0_rk*lambda*tau_in)
       
       dr_dim = ceiling(dr/(r(r_dim)-r(r_dim-1)))
       
       allocate(r2(r_dim+dr_dim))
       allocate(u0(r_dim+dr_dim,r_dim+dr_dim))
       call alustaCSpline2d(1,r_dim+dr_dim,r_dim+dr_dim)
       
       r2 = 0.0_rk
       r2(2:r_dim+1) = r
       dr = (r(r_dim)-r(r_dim-1))
       do i=r_dim+2,r_dim+dr_dim
          r2(i) = r2(i-1)+dr
       end do
       
       !r2(1:r_dim) = r
       !dr = (r(r_dim)-r(r_dim-1))
       !do i=r_dim+1,r_dim+dr_dim
       !   r2(i) = r2(i-1)+dr
       !end do

       if (modcusp==0) then
          cuspcor=1.0_rk
       else
          cuspcor=2.0_rk
       end if
       
       ! make the high temp. approx of the potential
       call u0_HighTempApprox(r2,u0,tau,r_dim+dr_dim)
       
       call initCSpline2dmod(r2,r2,u0,r_dim+dr_dim,r_dim+dr_dim,1,ZZ,lambda*cuspcor,tau) 
       !call initCSpline2d(r2,r2,u0,r_dim+dr_dim,r_dim+dr_dim,1) 
       
       write(*,'((A),(I2),(A3),(I2))') 'Generating Pair Potential ', PairPotNum, ' / ', NumPairPots
       
       write(*,*) 'Initial tau ', tau
       
       ! Matrix squaring
       do i = 1,NumSquares       
          ! calculate u_L for L=0     
          
          call MatrixSquaringIntegration(r2,u0,tau,r_dim+dr_dim)
          
          call initCSpline2dmod(r2,r2,u0,r_dim+dr_dim,r_dim+dr_dim,1,ZZ,lambda*cuspcor,tau)      
          !call initCSpline2d(r2,r2,u0,r_dim+dr_dim,r_dim+dr_dim,1) 
          
          tau = tau*2.0_rk
          write(*,*) '        tau ', tau
       end do
       
       write(*,*) ' '
       !write(*,*) 'At the Origin 1'
       !write(*,*) NearOrigin(0.0_rk,0.0_rk,tau),u0(1,1)
       !write(*,*) u0(1,1),u0(2,2),u0(3,3)
       write(*,*) ' '
       
       !! print diagonal of u_0, i.e. u0(r,r)

       if (1==2) then              
       write( i1, '(I6)'  ) PairPotNum
       open(123,file='u0.' // trim(adjustl(i1)) // '.txt', status='replace')
       do i=1,r_dim
          write(123,*) r2(i),u0(i,i)
       end do
       close(123)
       
       
       !! print u_0 surface ... gnuplot format (splot 'file' w l)
       write( i1, '(I6)'  ) PairPotNum
       open(123,file='u0s.' // trim(adjustl(i1)) // '.txt', status='replace')
       do i=1,r_dim
          do j=1,r_dim
             write(123,*) r2(i),r2(j),u0(i,j)
          end do
          write(123,*) ' '
       end do
       close(123)
       end if

       ! calculate pair action u(x,y) from u0
       u0 = PairPotential(r2,u0,tau,r_dim+dr_dim)
       
       write( i1, '(I6)'  ) PairPotNum
       open(123,file='ud.' // trim(adjustl(i1)) // '.txt', status='replace')
       do i=1,r_dim
          write(123,*) r2(i),u0(i,i)
       end do
       close(123)
       
       !! print u surface ... gnuplot format (splot 'file' w l)
       write( i1, '(I6)'  ) PairPotNum
       open(123,file='us.' // trim(adjustl(i1)) // '.txt', status='replace')
       do i=1,r_dim
          do j=1,r_dim
             write(123,'(3(ES20.10))') r2(i),r2(j),u0(i,j)
          end do
          write(123,*) ' '
       end do
       close(123)
       
    
       call initCSpline2d2(r2(1:r_dim),r2(1:r_dim),u0(1:r_dim,1:r_dim),r_dim,r_dim,PairPotNum)  

       call lopetaCSpline2d
       deallocate(u0)
       deallocate(r2)
      
    end if

  end subroutine GeneratePairPot





  
  subroutine GeneratePairPotWithBetaDeriv(lambda,ZZ,r,r_dim,PairPotNum,modcusp)
    use math    
    use matrixsquaring
    use periodicbc
    implicit none
    real(kind = rk), intent(in) :: lambda,ZZ
    integer, intent(in) :: r_dim,PairPotNum,modcusp
    real(kind = rk), dimension(r_dim), intent(in) :: r
    real(kind = rk) :: dr,tau,kerroin,z,u0approx,deltaU,cuspcor
    integer(kind = isoluku) :: M
    integer :: dr_dim,i,j,iovar
    real(kind = rk), dimension(:), allocatable :: r2
    real(kind = rk), dimension(:,:), allocatable :: u0,du0,ut
    real(kind = rk), dimension(r_dim,r_dim) :: u
    real(kind = rk), dimension(3) :: Box
    character(len=20) :: i1,help_text    
    logical :: generate,olemassa


    generate = .true.
    
    !if (lambda<0.01_rk) then 
    !   write(*,'((A),(I2),(A3),(I2))') 'Generating Pair Potential ', PairPotNum, ' / ', NumPairPots
    !   write(*,*) 'Bare Coulomb Potential: lambda = ', lambda
    !   write(*,*) ' '
    !   u=0.0_rk
    !   call initCSpline2d2(r,r,u,r_dim,r_dim,PairPotNum)
    !   call init_dU_dbeta(r,r,u,r_dim,r_dim,PairPotNum)
    !   return
    !end if

    if (ReadPairPotFile) then
       
       write( i1, '(I6)'  ) PairPotNum
       inquire(file='us.' // trim(adjustl(i1)) // '.txt', EXIST=olemassa )
       if (olemassa) then
          allocate(r2(r_dim))
          open(123,file='us.' // trim(adjustl(i1)) // '.txt')
          do i=1,r_dim
             do j=1,r_dim
                read(123,'(3(ES20.10))') r2(i),r2(j),u(i,j)
             end do
             read(123,'(A)',iostat=iovar) help_text
          end do
          close(123)
          
          call initCSpline2d2(r2,r2,u,r_dim,r_dim,PairPotNum)

          open(123,file='dus.' // trim(adjustl(i1)) // '.txt')
          do i=1,r_dim
             do j=1,r_dim
                read(123,'(3(ES20.10))') r2(i),r2(j),u(i,j)
             end do
             read(123,'(A)',iostat=iovar) help_text
          end do
          close(123)

          call init_dU_dbeta(r2,r2,u,r_dim,r_dim,PairPotNum)
          
          if (LongRangeON) then
             call GetTheBox(Box)
             call ReadLongRangeSpline(tau,ZZ,PairPotNum)
             if (PairPotNum==1) then
                call make_karray()
             end if
          end if

          generate = .false.
       else
          write(*,*) 'Cannot find pair potentials... thus, generating'
          write(*,*) ' '
          generate = .true.
       end if
    end if
    
    if (generate) then
          
       M = 2**NumSquares
       tau = tau_in/M    
       
       lambdaMS = lambda
       ZZMS = ZZ
       
       ! addition to the grid wanted
       dr = 20.0_rk*sqrt(2.0_rk*lambda*tau_in)
       
       dr_dim = ceiling(dr/(r(r_dim)-r(r_dim-1)))
       
       allocate(r2(r_dim+dr_dim))
       allocate(u0(r_dim+dr_dim,r_dim+dr_dim))
       allocate(du0(r_dim+dr_dim,r_dim+dr_dim))
       allocate(ut(r_dim+dr_dim,r_dim+dr_dim))
       call alustaCSpline2d(2,r_dim+dr_dim,r_dim+dr_dim)
       
       r2 = 0.0_rk
       r2(2:r_dim+1) = r
       dr = (r(r_dim)-r(r_dim-1))
       do i=r_dim+2,r_dim+dr_dim
          r2(i) = r2(i-1)+dr
       end do
       
       
       if (modcusp==0) then
          cuspcor=1.0_rk
       else
          cuspcor=2.0_rk
       end if
       
       ! make the high temp. approx of the potential
       call u0_HighTempApprox(r2,u0,tau,r_dim+dr_dim)
       
       call initCSpline2dmod(r2,r2,u0,r_dim+dr_dim,r_dim+dr_dim,1,ZZ,cuspcor*lambda,tau) 
       du0 = u0/tau
       call initCSpline2dmod(r2,r2,du0,r_dim+dr_dim,r_dim+dr_dim,2,ZZ,cuspcor*lambda,tau)     
       
       write(*,'((A),(I2),(A3),(I2))') 'Generating Pair Potential and Beta Derivative', PairPotNum, ' / ', NumPairPots
       
       write(*,*) 'Initial tau ', tau
       
       ! Matrix squaring
       do i = 1,NumSquares       
          ! calculate u_L for L=0     
          
          call MatrixSquaringIntegrationWithBetaDeriv(r2,u0,du0,tau,r_dim+dr_dim)
          tau = tau*2.0_rk
          
          call initCSpline2dmod(r2,r2,u0,r_dim+dr_dim,r_dim+dr_dim,1,ZZ,cuspcor*lambda,tau) 
          call initCSpline2dmod(r2,r2,du0,r_dim+dr_dim,r_dim+dr_dim,2,ZZ,cuspcor*lambda,tau)     
          
          if (i==3) then
             ut=u0
          else if (i==4) then
             du0=(u0-ut)/(tau-tau/2)
          end if
          
          write(*,*) '        tau ', tau
       end do
       
       
       write(*,*) ' '
       !write(*,*) 'At the Origin 1'
       !write(*,*)  u0(1,1)
       !write(*,*) u0(1,1),u0(2,2),u0(3,3)
       write(*,*) ' '
       
       !! print diagonal of u_0, i.e. u0(r,r)                                          
       if (1==2) then
       write( i1, '(I6)'  ) PairPotNum
       open(123,file='u0.' // trim(adjustl(i1)) // '.txt', status='replace')
       open(223,file='du0.' // trim(adjustl(i1)) // '.txt', status='replace')
       do i=1,r_dim
          write(123,*) r2(i),u0(i,i)
          write(223,*) r2(i),du0(i,i)
       end do
       close(123)
       close(223)
       end if

       if (1==2) then
       !print u0 and du0 surfaces ... gnuplot format (splot 'file' w l)
       write( i1, '(I6)'  ) PairPotNum
       open(123,file='u0s.' // trim(adjustl(i1)) // '.txt', status='replace')
       open(223,file='du0s.' // trim(adjustl(i1)) // '.txt', status='replace')
       do i=1,r_dim
          do j=1,r_dim
             write(123,*) r2(i),r2(j),u0(i,j)
             write(223,*) r2(i),r2(j),du0(i,j)
          end do
          write(123,*) ' '
          write(223,*) ' '
       end do
       close(123)
       close(223)
       end if
       
       if (LongRangeON) then
          call GetTheBox(Box)
          call makeLongRangeSpline(tau,ZZ,PairPotNum)
          if (PairPotNum==1) then
             call make_karray()
          end if          
       end if
       
       
       ! open(123,file='u0srdiag.dat')
       ! open(213,file='du0srdiag.dat')
       ! do i=1,r_dim
       !    write(123,*) r2(i),u0(i,i)-evalLR(r2(i),1)
       !    write(213,*) r2(i),du0(i,i)-evalLR(r2(i),2)
       ! end do
       ! close(123)
       ! close(213)
       
       
       !if (LongRangeON) then
       !   call PairPotentials2SRC(r2,u0,du0,tau,r_dim+dr_dim,PairPotNum)
       !else
          call PairPotentials2(r2,u0,du0,tau,r_dim+dr_dim)
       !end if


       !print u and du diagonals
       write( i1, '(I6)'  ) PairPotNum
       open(123,file='ud.' // trim(adjustl(i1)) // '.txt', status='replace')
       open(223,file='dud.' // trim(adjustl(i1)) // '.txt', status='replace')
       do i=1,r_dim
          write(123,'(2(ES20.10))') r2(i),u0(i,i)
          write(223,'(2(ES20.10))') r2(i),du0(i,i)
       end do
       close(123)
       close(223)

       !print u and du surfaces ... gnuplot format (splot 'file' w l)
       write( i1, '(I6)'  ) PairPotNum
       open(123,file='us.' // trim(adjustl(i1)) // '.txt', status='replace')
       open(223,file='dus.' // trim(adjustl(i1)) // '.txt', status='replace')
       do i=1,r_dim
          do j=1,r_dim
             write(123,'(3(ES20.10))') r2(i),r2(j),u0(i,j)
             write(223,'(3(ES20.10))') r2(i),r2(j),du0(i,j)
          end do
          write(123,*) ' '
          write(223,*) ' '
       end do
       close(123)
       close(223)

       
       ! potential
       call initCSpline2d2(r2(1:r_dim),r2(1:r_dim),u0(1:r_dim,1:r_dim),r_dim,r_dim,PairPotNum)
       
       ! beta derivative of potential
       call init_dU_dbeta(r2(1:r_dim),r2(1:r_dim),du0(1:r_dim,1:r_dim),r_dim,r_dim,PairPotNum)
       
       
       ! open(123,file='usrdiag.dat')
       ! open(213,file='dusrdiag.dat')
       ! do i=1,r_dim
       !    write(123,*) r2(i),u0(i,i)-evalLR(r2(i),1)
       !    write(213,*) r2(i),du0(i,i)-evalLR(r2(i),2)
       ! end do
       ! close(123)
       ! close(213)
       
       call lopetaCSpline2d
       deallocate(u0)
       deallocate(ut)
       deallocate(du0)
       deallocate(r2)
    end if
    

  end subroutine GeneratePairPotWithBetaDeriv

  function primApproximation(r,rp,ZZ,tau) result(V)
    use math
    implicit none
    real(kind = rk) :: r,rp,tau,ZZ,V,x,y
    
    if (r<1.0e-7) then
       x = 1.0e-7_rk
    else
       x = r
    end if
    if (rp<1.0e-7) then
       y = 1.0e-7_rk
    else
       y = rp
    end if

    V = 0.5_rk*tau*ZZ*(1.0_rk/x + 1.0_rk/y)
    
  end function primApproximation

  subroutine makeLongRangeSpline(tau,ZZ,PairPotNum)
    use periodicbc
    implicit none
    real(kind =rk), dimension(10000) :: x,u0erf
    real(kind = rk), dimension(3) :: Box
    real(kind = rk) :: kerroin,tau,ZZ,ee,maxabsZZ
    integer :: i,PairPotNum
    character(len=20) :: i1

    call GetTheBox(Box)

    !if (alpha .eq. 0.0_rk) then
    ee = 1.0e-45_rk
    maxabsZZ = 1.0_rk
    alpha = 2.0_rk/Box(1)*sqrt(-log(sqrt(pi)/2*ee/maxabsZZ))
    !end if

    x=linspace(0.0_rk,sqrt(3.0_rk)*Box(1)/2+5.0_rk,10000)
    
    u0erf=0.0_rk
    do i=1,10000
       u0erf(i) = eval_u0at(x(i),x(i),1,tau,ZZ)*erff(alpha*x(i))
    end do
    
    !! print long range correction
    write( i1, '(I6)'  ) PairPotNum
    open(123,file='uL.' // trim(adjustl(i1)) // '.txt', status='replace')
    do i=1,10000
       write(123,'(2(ES20.10))') x(i),u0erf(i)
    end do
    close(123)
    
    call alustaLRcspline1d(PairpotNum,2,10000)
    
    call initLRcspline1d(x,u0erf,10000,PairpotNum,1)            
    
    u0erf=0.0_rk
    do i=1,10000
       u0erf(i) = eval_u0at(x(i),x(i),2,tau,ZZ)*erff(alpha*x(i))
    end do
    
    call initLRcspline1d(x,u0erf,10000,PairpotNum,2)            
    
    !! print long range correction
    write( i1, '(I6)'  ) PairPotNum
    open(223,file='duL.' // trim(adjustl(i1)) // '.txt', status='replace')
    do i=1,10000
       write(223,'(2(ES20.10))') x(i),u0erf(i)
    end do
    close(223)    
    
  end subroutine makeLongRangeSpline

  subroutine ReadLongRangeSpline(tau,ZZ,PairPotNum)
    use periodicbc
    implicit none
    real(kind =rk), dimension(10000) :: x,u0erf
    real(kind = rk), dimension(3) :: Box
    real(kind = rk) :: kerroin,tau,ZZ,ee,maxabsZZ
    integer :: i,PairPotNum
    character(len=20) :: i1

    call GetTheBox(Box)

    !if (alpha .eq. 0.0_rk) then
    ee = 1.0e-45_rk
    maxabsZZ = 1.0_rk
    alpha = 2.0_rk/Box(1)*sqrt(-log(sqrt(pi)/2*ee/maxabsZZ))
    !end if

    !x=linspace(0.0_rk,sqrt(3.0_rk)*Box(1)/2+5.0_rk,10000)
    x=0.0_rk
    u0erf = 0.0_rk
        
    !! read long range correction
    write( i1, '(I6)'  ) PairPotNum
    open(123,file='uL.' // trim(adjustl(i1)) // '.txt')
    do i=1,10000
       read(123,'(2(ES20.10))') x(i),u0erf(i)
    end do
    close(123)
    
    call alustaLRcspline1d(PairpotNum,2,10000)
    
    call initLRcspline1d(x,u0erf,10000,PairpotNum,1)            
    
    u0erf=0.0_rk

    !! read long range correction
    write( i1, '(I6)'  ) PairPotNum
    open(223,file='duL.' // trim(adjustl(i1)) // '.txt')
    do i=1,10000
       read(223,'(2(ES20.10))') x(i),u0erf(i)
    end do
    close(223)    
    
    call initLRcspline1d(x,u0erf,10000,PairpotNum,2)            
    
    
  end subroutine ReadLongRangeSpline



  function eval_u0at(x,y,pot,tau,ZZ) result(U)
    implicit none
    real(kind=rk), intent(in) :: x,y
    real(kind=rk) :: U,tau,ZZ,rx,ry,Vlr
    integer,intent(in) :: pot
    
    U = 0.0_rk
    rx = x
    ry = y        

    if (rx<=Ucspline(pot)%x(Ucspline(pot)%x_dim) &
         .and. ry<=Ucspline(pot)%y(Ucspline(pot)%y_dim)) then
       
       if (rx<Ucspline(pot)%x(1)) then 
          rx = Ucspline(pot)%x(1)
       end if
       
       if (ry<Ucspline(pot)%y(1)) then
          ry = Ucspline(pot)%y(1)
       end if
       
       U = eval2dcspline(rx,ry,pot)
       
    else       
       
       if (rx>Ucspline(pot)%x(Ucspline(pot)%x_dim) &
            .and. ry>Ucspline(pot)%y(Ucspline(pot)%y_dim)) then
          U = (eval2dcspline(Ucspline(pot)%x(Ucspline(pot)%x_dim),&
               Ucspline(pot)%y(Ucspline(pot)%y_dim),pot))/&
               (primApproximation(Ucspline(pot)%x(Ucspline(pot)%x_dim),&
               Ucspline(pot)%y(Ucspline(pot)%y_dim),tau,ZZ))&
               *(primApproximation(x,y,tau,ZZ))
       elseif (rx<=Ucspline(pot)%x(Ucspline(pot)%x_dim) &
            .and. ry>Ucspline(pot)%y(Ucspline(pot)%y_dim)) then
       
          if (rx<Ucspline(pot)%x(1)) then 
             rx = Ucspline(pot)%x(1)
          end if
          U = (eval2dcspline(rx,Ucspline(pot)%y(Ucspline(pot)%y_dim),pot))/&
               (primApproximation(rx,Ucspline(pot)%y(Ucspline(pot)%y_dim),tau,ZZ))&
               *(primApproximation(rx,ry,tau,ZZ))
       else
          if (ry<Ucspline(pot)%y(1)) then
             ry = Ucspline(pot)%y(1)
          end if
          U = (eval2dcspline(Ucspline(pot)%x(Ucspline(pot)%x_dim),ry,pot))/&
               (primApproximation(Ucspline(pot)%x(Ucspline(pot)%x_dim),ry,tau,ZZ))&
               *(primApproximation(rx,ry,tau,ZZ))
       end if
    end if
    
    
  end function eval_u0at


  function LRpotential(r,rp,pot,spline_indx) result(LRpot)
    implicit none
    real(kind = rk), intent(in) :: r,rp
    !real(kind = rk), dimension(3), intent(in) :: rv,rpv
    integer, intent(in) :: spline_indx,pot
    real(kind =rk) :: LRpot
    integer :: i
    
    ! short range correction    

    if (LongRangeOn) then
       LRpot = -(evalLR(r,pot,spline_indx)+evalLR(rp,pot,spline_indx))/2
    else
       LRpot = 0.0_rk
    end if

  end function LRpotential

  function LRcos(ZZ,rv,rpv) result(LRpot)
    implicit none
    real(kind = rk), intent(in) :: ZZ
    real(kind = rk), dimension(3), intent(in) :: rv,rpv
    real(kind =rk) :: LRpot
    integer :: i    
    
    LRpot = 0.0_rk

    if (LongRangeOn) then
       ! long range part
       
       LRpot = 0.0_rk
       do i=1,karray_dim1
          !if (karray(i,5)<kcutoff) then ! no need for this anymore
             LRpot = LRpot + karray(i,4)*(cos(dot_product(karray(i,1:3),rv))+cos(dot_product(karray(i,1:3),rpv)))
          !end if
       end do
       LRpot = LRpot*ZZ
    end if
    
  end function LRcos

  

  subroutine make_karray
    use periodicbc
    implicit none
    real(kind = rk), dimension(3) :: Box,kv
    real(kind = rk) :: kerroin,k2,k
    integer :: i,lukum,loppu,coef,lukum2
    integer :: n_max,n_x,n_y,n_z

    call GetTheBox(Box)
    
    !coef = floor(20.0_rk/pi)+1
    !kcutoff = real(coef,kind=rk)*2.0_rk*pi/Box(1)
    !kcutoff = 14.0_rk*pi/Box(1)    

    if (kcutoff_in< 1.0e-10_rk) then
       kcutoff = 12.0_rk*pi/Box(1)
       n_max = 6
    else
       kcutoff = kcutoff_in
       n_max = floor(kcutoff/2/pi*Box(1))+1
    end if
    !n_max = coef

    lukum = 0
    lukum2 = 0
    do n_x=-n_max,n_max
       do n_y=-n_max,n_max
          do n_z=-n_max,n_max
             
             kv=(/real(n_x,kind=rk)/Box(1),real(n_y,kind=rk)/Box(2),real(n_z,kind=rk)/Box(3)/)*2.0_rk*pi
             k = sqrt(sum(kv**2))
             if (k<kcutoff) then
                lukum = lukum+1
             end if
             lukum2 = lukum2 + 1
          end do
       end do
    end do

    allocate(karray((lukum-1)/2,5))
   

    lukum = (lukum-1)/2
    lukum2 = (lukum2-1)/2
    karray = 0.0_rk
    karray_dim1 = lukum


    loppu = 1
    i=1
    do n_x=-n_max,n_max
       do n_y=-n_max,n_max
          do n_z=-n_max,n_max

             if (i<= lukum2) then
                kv=(/real(n_x,kind=rk)/Box(1),real(n_y,kind=rk)/Box(2),real(n_z,kind=rk)/Box(3)/)*2.0_rk*pi
                k = sqrt(sum(kv**2))
                if (k<kcutoff) then
                   if (loppu <= lukum) then
                      karray(loppu,1:3) = kv
                      karray(loppu,5) = k
                   end if
                   loppu=loppu+1
                end if
             end if
             i=i+1
          end do
       end do
    end do
    
    !write(*,*) loppu-1,lukum
    
    kerroin = 4.0_rk*pi/Box(1)/Box(2)/Box(3)
    
    do i=1,lukum
       k2 = karray(i,5)**2
       karray(i,4) = kerroin*exp(-k2/4/alpha**2)/k2
    end do
               
  end subroutine make_karray


  

  


  subroutine make_karrayOld
    use periodicbc
    implicit none
    real(kind = rk), dimension(3) :: Box
    real(kind = rk) :: kerroin,k2
    integer :: i,lukum,loppu,coef
    integer :: n_max,n_x,n_y,n_z
    logical, save :: firsttime=.true.

    
    if (firsttime) then
       call GetTheBox(Box)
       
       !coef = floor(10.0_rk/pi)+1
       !kcutoff = real(coef,kind=rk)*2.0_rk*pi/Box(1)
       kcutoff = 12.0_rk*pi/Box(1)    
       
       n_max = 6
       !n_max = coef
       
       lukum = 0
       do n_x=-n_max,n_max
          do n_y=-n_max,n_max
             do n_z=-n_max,n_max
                lukum = lukum+1
             end do
          end do
       end do
       
       allocate(karray((lukum-1)/2,5))
       
       
       lukum = (lukum-1)/2
       
       karray = 0.0_rk
       karray_dim1 = lukum
       
       
       loppu = 1
       do n_x=-n_max,n_max
          do n_y=-n_max,n_max
             do n_z=-n_max,n_max
                
                if (loppu <= lukum) then
                   karray(loppu,1:3) = (/real(n_x,kind=rk)/Box(1),real(n_y,kind=rk)/Box(2),real(n_z,kind=rk)/Box(3)/)*2.0_rk*pi
                   karray(loppu,5) = sqrt(sum(karray(loppu,1:3)**2))
                end if
                loppu=loppu+1
             end do
          end do
       end do
       
       kerroin = 4.0_rk*pi/Box(1)/Box(2)/Box(3)
       
       do i=1,lukum
          k2 = karray(i,5)**2
          karray(i,4) = kerroin*exp(-k2/4/alpha**2)/k2
       end do
       
    end if


  end subroutine make_karrayOld


  

END MODULE init
