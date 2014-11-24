module math
  
  integer, parameter :: r4k = KIND(1.0)                  !single precision
  integer, parameter :: r8k = SELECTED_REAL_KIND(2*PRECISION(1.0_r4k))
  integer, parameter :: rk=r8k !this is the real kinddd we are using 

  INTEGER, PARAMETER :: isoluku  = SELECTED_INT_KIND(12)
  
  !pi = 4.0_rk*atan(1.0_rk)
  real(kind = rk), parameter :: pi = 3.14159265358979_rk 
  real(kind = rk), parameter :: kB = 3.1668152e-6_rk
  !real(kind = rk), save :: dielectric_constant
  !real(kind = rk), save :: e_h_exchange_scale

  TYPE:: CSpline1d
     
     integer :: x_dim
     real(kind=rk), dimension(:), allocatable :: F
     real(kind=rk), dimension(:), allocatable :: Fx
     real(kind=rk), dimension(:), allocatable :: x
     real(kind=rk), dimension(:), allocatable :: hx
     
  END TYPE CSpline1d

  type :: ShortRangeCorrection
     type( CSpline1d ), dimension(:), allocatable  :: LRcspline
  end type ShortRangeCorrection
  type( ShortRangeCorrection ), dimension(:), allocatable  :: SRcorrection
  


  TYPE:: CSpline2d
     
     integer :: x_dim
     integer :: y_dim
     real(kind=rk), dimension(:,:), allocatable :: F
     real(kind=rk), dimension(:,:), allocatable :: Fx
     real(kind=rk), dimension(:,:), allocatable :: Fy
     real(kind=rk), dimension(:,:), allocatable :: Fxy
     real(kind=rk), dimension(:), allocatable :: x
     real(kind=rk), dimension(:), allocatable :: y
     real(kind=rk), dimension(:), allocatable :: hx
     real(kind=rk), dimension(:), allocatable :: hy
     
  END TYPE CSpline2d
  
  type( CSpline2d ), dimension(:), allocatable  :: Ucspline
  type( CSpline2d ), dimension(:), allocatable  :: Ucspline2
  type( CSpline2d ), dimension(:), allocatable  :: dU_dbeta

  real(kind = rk), private :: GaussKronrodValue
  real(kind = rk), private :: GaussKronrodMinStep
  real(kind = rk), save :: lambdaMS,ZZMS

contains


  subroutine alustaLRcspline1d(pot,f_dim,x_dim)
    implicit none
    integer :: i,x_dim,f_dim,pot
    
    allocate(SRcorrection(pot)%LRcspline(f_dim))
    
    do i = 1,f_dim
       allocate(SRcorrection(pot)%LRcspline(i)%F(x_dim))
       allocate(SRcorrection(pot)%LRcspline(i)%Fx(x_dim))
       allocate(SRcorrection(pot)%LRcspline(i)%x(x_dim))
       allocate(SRcorrection(pot)%LRcspline(i)%hx(x_dim-1))

       SRcorrection(pot)%LRcspline(i)%x_dim = x_dim
       SRcorrection(pot)%LRcspline(i)%F = 0.0_rk
       SRcorrection(pot)%LRcspline(i)%Fx = 0.0_rk
       SRcorrection(pot)%LRcspline(i)%x = 0.0_rk
       SRcorrection(pot)%LRcspline(i)%hx = 0.0_rk

    end do
    
  end subroutine alustaLRcspline1d

  subroutine lopetaLRcspline1d
    implicit none
    
    deallocate(SRcorrection)
    
  end subroutine lopetaLRcspline1d


  
  subroutine alustaCSpline2d(f_dim,x_dim,y_dim)
    implicit none
    integer :: i,x_dim,y_dim,f_dim
    
    allocate(Ucspline(f_dim))
    
    do i  = 1,f_dim
       allocate(Ucspline(i)%F(x_dim,y_dim))
       allocate(Ucspline(i)%Fx(x_dim,y_dim))
       allocate(Ucspline(i)%Fy(x_dim,y_dim))
       allocate(Ucspline(i)%Fxy(x_dim,y_dim))
       allocate(Ucspline(i)%x(x_dim))
       allocate(Ucspline(i)%y(y_dim))
       allocate(Ucspline(i)%hx(x_dim-1))
       allocate(Ucspline(i)%hy(y_dim-1))

       Ucspline(i)%x_dim = x_dim
       Ucspline(i)%y_dim = y_dim
       Ucspline(i)%F = 0.0_rk
       Ucspline(i)%Fx = 0.0_rk
       Ucspline(i)%Fy = 0.0_rk
       Ucspline(i)%Fxy = 0.0_rk
       Ucspline(i)%x = 0.0_rk
       Ucspline(i)%y = 0.0_rk
       Ucspline(i)%hx = 0.0_rk
       Ucspline(i)%hy = 0.0_rk
    end do
    
  end subroutine alustaCSpline2d

  subroutine lopetaCSpline2d
    implicit none
    
    deallocate(Ucspline)
    
  end subroutine lopetaCSpline2d

  
  subroutine alustaCSpline2d2(f_dim,x_dim,y_dim)
    implicit none
    integer :: i,x_dim,y_dim,f_dim
    
    allocate(Ucspline2(f_dim))
    
    do i  = 1,f_dim
       allocate(Ucspline2(i)%F(x_dim,y_dim))
       allocate(Ucspline2(i)%Fx(x_dim,y_dim))
       allocate(Ucspline2(i)%Fy(x_dim,y_dim))
       allocate(Ucspline2(i)%Fxy(x_dim,y_dim))
       allocate(Ucspline2(i)%x(x_dim))
       allocate(Ucspline2(i)%y(y_dim))
       allocate(Ucspline2(i)%hx(x_dim-1))
       allocate(Ucspline2(i)%hy(y_dim-1))

       Ucspline2(i)%x_dim = x_dim
       Ucspline2(i)%y_dim = y_dim
       Ucspline2(i)%F = 0.0_rk
       Ucspline2(i)%Fx = 0.0_rk
       Ucspline2(i)%Fy = 0.0_rk
       Ucspline2(i)%Fxy = 0.0_rk
       Ucspline2(i)%x = 0.0_rk
       Ucspline2(i)%y = 0.0_rk
       Ucspline2(i)%hx = 0.0_rk
       Ucspline2(i)%hy = 0.0_rk
    end do
    
  end subroutine alustaCSpline2d2

  
  subroutine alusta_dU_dbeta(f_dim,x_dim,y_dim)
    implicit none
    integer :: i,x_dim,y_dim,f_dim
    
    allocate(dU_dbeta(f_dim))
    
    do i  = 1,f_dim
       allocate(dU_dbeta(i)%F(x_dim,y_dim))
       allocate(dU_dbeta(i)%Fx(x_dim,y_dim))
       allocate(dU_dbeta(i)%Fy(x_dim,y_dim))
       allocate(dU_dbeta(i)%Fxy(x_dim,y_dim))
       allocate(dU_dbeta(i)%x(x_dim))
       allocate(dU_dbeta(i)%y(y_dim))
       allocate(dU_dbeta(i)%hx(x_dim-1))
       allocate(dU_dbeta(i)%hy(y_dim-1))

       dU_dbeta(i)%x_dim = x_dim
       dU_dbeta(i)%y_dim = y_dim
       dU_dbeta(i)%F = 0.0_rk
       dU_dbeta(i)%Fx = 0.0_rk
       dU_dbeta(i)%Fy = 0.0_rk
       dU_dbeta(i)%Fxy = 0.0_rk
       dU_dbeta(i)%x = 0.0_rk
       dU_dbeta(i)%y = 0.0_rk
       dU_dbeta(i)%hx = 0.0_rk
       dU_dbeta(i)%hy = 0.0_rk
    end do
    
  end subroutine alusta_dU_dbeta


  subroutine lopetaCSpline2d2
    implicit none
    
    deallocate(Ucspline2)
    
  end subroutine lopetaCSpline2d2

  subroutine lopeta_dU_dbeta
    implicit none
    
    deallocate(dU_dbeta)
    
  end subroutine lopeta_dU_dbeta



!   subroutine initCSpline1d(r,f,h,r_dim,fp) 
!     implicit none
!     integer, intent(in) :: r_dim
!     integer, dimension(r_dim) :: IPIV
!     integer :: INFO
!     real(kind = rk), dimension(r_dim), intent(in)   :: r,f
!     real(kind = rk), dimension(r_dim-1),intent(in)  :: h
!     real(kind = rk), dimension(r_dim), intent(out)  :: fp
!     !real(kind = rk), dimension(r_dim)   :: d
!     real(kind = rk), dimension(r_dim,r_dim)   :: A
!     integer :: i    
    
!     A = 0.0_rk
!     fp(1) = (f(2)-f(1))/(r(2)-r(1))
!     A(1,1) = 1.0_rk
!     do i = 2,r_dim-1  
!        fp(i) = 6.0_rk*(h(i)/h(i-1)-h(i-1)/h(i))*f(i) - 6.0_rk*h(i)/h(i-1)*f(i-1) &
!             + 6.0_rk*h(i-1)/h(i)*f(i+1)
!        A(i,i-1:i+1) = (/2.0_rk*h(i), 4.0_rk*(h(i)+h(i-1)), 2.0_rk*h(i-1)/)
!     end do
!     A(r_dim,r_dim) = 1.0_rk
!     fp(r_dim) = (f(r_dim)-f(r_dim-1))/(r(r_dim)-r(r_dim-1))
        
!     call DGESV( r_dim, 1, A, r_dim, IPIV, fp, r_dim, INFO )


!   end subroutine initCSpline1d


  subroutine initCSpline1d(r,f,h,r_dim,fp) 
    implicit none
    integer, intent(in) :: r_dim
    real(kind = rk), dimension(r_dim), intent(in)   :: r,f
    real(kind = rk), dimension(r_dim-1),intent(in)  :: h
    real(kind = rk), dimension(r_dim), intent(out)  :: fp
    real(kind = rk), dimension(r_dim)   :: a,b,c,d
    real(kind = rk) :: id
    integer :: i    
    real(kind = rk) :: h1,h2

    h1 = r(2)-r(1)
    h2 = r(3)-r(2)
    
    
    !d(1) = (f(2)-f(1))/(r(2)-r(1))
    a(1) = 1.0_rk ! not used
    b(1) = 1.0_rk

    ! for complete boundary conditions
    c(1) = 0.0_rk 
    d(1) = (f(2)-f(1))/(r(2)-r(1))
    !d(1) = (-h1**2*f(3)+(h1+h2)**2*f(2)-(2.0_rk*h1*h2+h2**2)*f(1))/&
    !     (h1**2*h2+h1*h2**2)
    ! ---------------
    
    ! for natural boundary conditions    
    !c(1) = 0.5_rk 
    !d(1) = 1.5_rk*(f(2)-f(1))/(r(2)-r(1))
    ! ---------------


    do i = 2,r_dim-1
       d(i) = 6.0_rk*(h(i)/h(i-1)-h(i-1)/h(i))*f(i) - 6.0_rk*h(i)/h(i-1)*f(i-1) &
            + 6.0_rk*h(i-1)/h(i)*f(i+1)
       a(i) = 2.0_rk*h(i)
       b(i) = 4.0_rk*(h(i)+h(i-1))
       c(i) = 2.0_rk*h(i-1)
    end do
    
    ! for complete boundary conditions
    a(r_dim) = 0.0_rk
    d(r_dim) = (f(r_dim)-f(r_dim-1))/(r(r_dim)-r(r_dim-1))
    ! ------------- 

    ! for natural boundary conditions
    !a(r_dim) = 0.5_rk
    !d(r_dim) = 1.5_rk*(f(r_dim)-f(r_dim-1))/(r(r_dim)-r(r_dim-1))
    ! --------------


    b(r_dim) = 1.0_rk
    c(r_dim) = 1.0_rk ! not used    



    ! Solve tridiagonal eq. A*fp = d
    c(1) = c(1)/b(1)
    d(1) = d(1)/b(1)
    do i = 2,r_dim-1
       id = b(i) - c(i-1)*a(i)
       c(i) = c(i)/id
       d(i) = (d(i) - d(i-1)*a(i))/id
    end do

    fp(r_dim) = d(r_dim)
    do i = r_dim-1,1,-1
       fp(i) = d(i)-c(i)*fp(i+1)
    end do

  end subroutine initCSpline1d

  subroutine initLRcspline1d(x,F,x_dim,pot,spline_indx) 
    implicit none
    integer, intent(in) :: x_dim,spline_indx
    real(kind = rk), dimension(x_dim), intent(in)   :: F
    real(kind = rk), dimension(x_dim), intent(in)   :: x
    integer :: i,pot
    
    SRcorrection(pot)%LRcspline(spline_indx)%x = x
    
    SRcorrection(pot)%LRcspline(spline_indx)%hx = x(2:x_dim)-x(1:x_dim-1)

    SRcorrection(pot)%LRcspline(spline_indx)%F = F 
    
    call initCSpline1d(x,F,SRcorrection(pot)%LRcspline(spline_indx)%hx,x_dim,SRcorrection(pot)%LRcspline(spline_indx)%Fx)     
    
  end subroutine initLRcspline1d

  subroutine initCSpline1dtest22(r,f,h,r_dim,fp,ZZ,lambda,ind,tau) 
    implicit none
    integer, intent(in) :: r_dim,ind
    real(kind = rk), dimension(r_dim), intent(in)   :: r
    real(kind = rk), dimension(r_dim), intent(inout)   :: f
    real(kind = rk), dimension(r_dim-1),intent(in)  :: h
    real(kind = rk), dimension(r_dim), intent(out)  :: fp
    real(kind = rk), dimension(r_dim)   :: a,b,c,d
    real(kind = rk) :: id,ZZ,lambda,delta,d2,uusi,der,der2,tau
    integer :: i,count,k
    real(kind = rk) :: h1,h2

    h1 = r(2)-r(1)
    h2 = r(3)-r(2)

    if (1==2) then
    do i = r_dim-3,1,-1
       
       do k=1,5
          if (abs(f(i))<abs(f(i+1))) then
             uusi=extrapolatePP(r(i:i+3),f(i:i+3),r(i),4)
             f(i) = uusi
          end if
          if (ind==1) then
             der = (f(i)-f(i+1))/(r(i)-r(i+1))
             if (i>1) then
                der2 = (f(i)-f(i-1))/(r(i)-r(i-1))
                if (abs(der)>abs(der2)) then
                   uusi=extrapolatePP(r(i-1:i+3),f(i-1:i+3),r(i),5)
                   f(i) = uusi             
                end if
             end if
             
             
             !             if (tau>abs(r(i)-r(i+1))) then
             der = (f(i)-f(i+1))/(r(i)-r(i+1))
             if (abs(der)>abs(ZZ/lambda/2)) then
                f(i) = f(i+1)-ZZ/lambda/2*(r(i)-r(i+1))
             end if
             !             end if
          end if
       end do       
    end do
    end if

    a(1) = 1.0_rk ! not used
    b(1) = 1.0_rk

    ! for complete boundary conditions
    c(1) = 0.0_rk 
    !do k=1,10
    f(1) = extrapolatePP(r(2:4),f(2:4),r(1),3)
    !end do
    d(1) = (f(2)-f(1))/(r(2)-r(1))
    ! if (ind==1) then
    !    if (abs(d(1)-ZZ/lambda/2)>1.0e-3_rk) then
    !       !d(1) = -ZZ/lambda/2
    !       do k=1,5
    !          f(1) = extrapolatePP(r(1:4),f(1:4),r(1),4)
    !       end do
    !       d(1) = (f(2)-f(1))/(r(2)-r(1))
    !       !f(1) = d(1)*(r(1)-r(2))+f(2)
    !    end if
    ! else
    !    if (abs(d(1))>1.0e-3_rk) then
    !       d(1) = 0.0_rk
    !       f(1) = d(1)*(r(1)-r(2))+f(2)
    !    end if
    ! end if

    !d(1) = (-h1**2*f(3)+(h1+h2)**2*f(2)-(2.0_rk*h1*h2+h2**2)*f(1))/&
    !     (h1**2*h2+h1*h2**2)
    ! ---------------
    
    ! for natural boundary conditions    
    !c(1) = 0.5_rk 
    !d(1) = 1.5_rk*(f(2)-f(1))/(r(2)-r(1))
    ! ---------------
    
    do i = 2,r_dim-1
       d(i) = 6.0_rk*(h(i)/h(i-1)-h(i-1)/h(i))*f(i) - 6.0_rk*h(i)/h(i-1)*f(i-1) &
            + 6.0_rk*h(i-1)/h(i)*f(i+1)
       a(i) = 2.0_rk*h(i)
       b(i) = 4.0_rk*(h(i)+h(i-1))
       c(i) = 2.0_rk*h(i-1)
    end do
    
    ! for complete boundary conditions
    a(r_dim) = 0.0_rk
    d(r_dim) = (f(r_dim)-f(r_dim-1))/(r(r_dim)-r(r_dim-1))
    ! ------------- 

    ! for natural boundary conditions
    !a(r_dim) = 0.5_rk
    !d(r_dim) = 1.5_rk*(f(r_dim)-f(r_dim-1))/(r(r_dim)-r(r_dim-1))
    ! --------------


    b(r_dim) = 1.0_rk
    c(r_dim) = 1.0_rk ! not used    



    ! Solve tridiagonal eq. A*fp = d
    c(1) = c(1)/b(1)
    d(1) = d(1)/b(1)
    do i = 2,r_dim-1
       id = b(i) - c(i-1)*a(i)
       c(i) = c(i)/id
       d(i) = (d(i) - d(i-1)*a(i))/id
    end do


    fp(r_dim) = d(r_dim)
    do i = r_dim-1,1,-1
       fp(i) = d(i)-c(i)*fp(i+1)
    end do


  end subroutine initCSpline1dtest22



  
  subroutine initCSpline2d(x,y,F,x_dim,y_dim,spline_indx) 
    implicit none
    integer, intent(in) :: x_dim,y_dim,spline_indx
    real(kind = rk), dimension(x_dim,y_dim), intent(in)   :: F
    real(kind = rk), dimension(x_dim), intent(in)   :: x
    real(kind = rk), dimension(y_dim), intent(in)   :: y
    integer :: i
    
    Ucspline(spline_indx)%x = x
    Ucspline(spline_indx)%y = y
    
    Ucspline(spline_indx)%hx = x(2:x_dim)-x(1:x_dim-1)
    Ucspline(spline_indx)%hy = y(2:y_dim)-y(1:y_dim-1)

    Ucspline(spline_indx)%F = F    
    
    do i = 1,maxval((/x_dim,y_dim/))
       if (i<=y_dim) then
          call initCSpline1d(x,F(:,i),Ucspline(spline_indx)%hx,x_dim,Ucspline(spline_indx)%Fx(:,i)) 
       end if
       if (i<=x_dim) then
          call initCSpline1d(y,F(i,:),Ucspline(spline_indx)%hy,y_dim,Ucspline(spline_indx)%Fy(i,:)) 
       end if
    end do

    do i=1,y_dim
       call initCSpline1d(x,Ucspline(spline_indx)%Fy(:,i),Ucspline(spline_indx)%hx,x_dim,Ucspline(spline_indx)%Fxy(:,i))
    end do
    
  end subroutine initCSpline2d

  
  subroutine initCSpline2d2(x,y,F,x_dim,y_dim,spline_indx) 
    implicit none
    integer, intent(in) :: x_dim,y_dim,spline_indx
    real(kind = rk), dimension(x_dim,y_dim), intent(in)   :: F
    real(kind = rk), dimension(x_dim), intent(in)   :: x
    real(kind = rk), dimension(y_dim), intent(in)   :: y
    integer :: i
    
    Ucspline2(spline_indx)%x = x
    Ucspline2(spline_indx)%y = y
    
    Ucspline2(spline_indx)%hx = x(2:x_dim)-x(1:x_dim-1)
    Ucspline2(spline_indx)%hy = y(2:y_dim)-y(1:y_dim-1)

    Ucspline2(spline_indx)%F = F    
    
    do i = 1,maxval((/x_dim,y_dim/))
       if (i<=y_dim) then
          call initCSpline1d(x,F(:,i),Ucspline2(spline_indx)%hx,x_dim,Ucspline2(spline_indx)%Fx(:,i)) 
       end if
       if (i<=x_dim) then
          call initCSpline1d(y,F(i,:),Ucspline2(spline_indx)%hy,y_dim,Ucspline2(spline_indx)%Fy(i,:)) 
       end if
    end do

    do i=1,y_dim
       call initCSpline1d(x,Ucspline2(spline_indx)%Fy(:,i),Ucspline2(spline_indx)%hx,x_dim,Ucspline2(spline_indx)%Fxy(:,i))
    end do
    
  end subroutine initCSpline2d2

    
  subroutine initCSpline2dmod(x,y,F,x_dim,y_dim,spline_indx,ZZ,lambda,tau) 
    implicit none
    integer, intent(in) :: x_dim,y_dim,spline_indx
    real(kind = rk), dimension(x_dim,y_dim), intent(inout)   :: F
    !real(kind = rk), dimension(x_dim,y_dim)   :: temp
    real(kind = rk), dimension(x_dim), intent(in)   :: x
    real(kind = rk), dimension(y_dim), intent(in)   :: y
    real(kind =rk), intent(in) :: ZZ,lambda,tau
    integer :: i,j,k,ii,ii1,ii2,count
    
    Ucspline(spline_indx)%x = x
    Ucspline(spline_indx)%y = y
    
    Ucspline(spline_indx)%hx = x(2:x_dim)-x(1:x_dim-1)
    Ucspline(spline_indx)%hy = y(2:y_dim)-y(1:y_dim-1)

    !temp = F
    do i = 1,x_dim
       call initCSpline1dtest22(x,F(:,i),Ucspline(spline_indx)%hx,x_dim,Ucspline(spline_indx)%Fx(:,i),ZZ,lambda,spline_indx,tau) 
    end do
    
    do i = 1,y_dim
       call initCSpline1dtest22(y,F(i,:),Ucspline(spline_indx)%hy,y_dim,Ucspline(spline_indx)%Fy(i,:),ZZ,lambda,spline_indx,tau) 
    end do
        
    do i = 1,maxval((/x_dim,y_dim/))
       if (i<=y_dim) then
          call initCSpline1d(x,F(:,i),Ucspline(spline_indx)%hx,x_dim,Ucspline(spline_indx)%Fx(:,i)) 
       end if
       if (i<=x_dim) then
          call initCSpline1d(y,F(i,:),Ucspline(spline_indx)%hy,y_dim,Ucspline(spline_indx)%Fy(i,:)) 
       end if
    end do
    
    do i=1,y_dim
       call initCSpline1d(x,Ucspline(spline_indx)%Fy(:,i),Ucspline(spline_indx)%hx,x_dim,Ucspline(spline_indx)%Fxy(:,i))
    end do
    Ucspline(spline_indx)%F = F
    
  end subroutine initCSpline2dmod

  

  
  subroutine init_dU_dbeta(x,y,F,x_dim,y_dim,spline_indx) 
    implicit none
    integer, intent(in) :: x_dim,y_dim,spline_indx
    real(kind = rk), dimension(x_dim,y_dim), intent(in)   :: F
    real(kind = rk), dimension(x_dim), intent(in)   :: x
    real(kind = rk), dimension(y_dim), intent(in)   :: y
    integer :: i
    
    dU_dbeta(spline_indx)%x = x
    dU_dbeta(spline_indx)%y = y
    
    dU_dbeta(spline_indx)%hx = x(2:x_dim)-x(1:x_dim-1)
    dU_dbeta(spline_indx)%hy = y(2:y_dim)-y(1:y_dim-1)

    dU_dbeta(spline_indx)%F = F    
    
    do i = 1,maxval((/x_dim,y_dim/))
       if (i<=y_dim) then
          call initCSpline1d(x,F(:,i),dU_dbeta(spline_indx)%hx,x_dim,dU_dbeta(spline_indx)%Fx(:,i)) 
       end if
       if (i<=x_dim) then
          call initCSpline1d(y,F(i,:),dU_dbeta(spline_indx)%hy,y_dim,dU_dbeta(spline_indx)%Fy(i,:)) 
       end if
    end do

    do i=1,y_dim
       call initCSpline1d(x,dU_dbeta(spline_indx)%Fy(:,i),dU_dbeta(spline_indx)%hx,x_dim,dU_dbeta(spline_indx)%Fxy(:,i))
    end do
    
  end subroutine init_dU_dbeta

  
  function extrapolatePP(r,u,r0,r_dim) result(u0)
    implicit none
    integer :: r_dim
    real(kind = rk), dimension(r_dim), intent(in) :: r,u
    real(kind = rk), intent(in) :: r0
    real(kind = rk), dimension(3,3) :: a,ainv
    real(kind = rk), dimension(3) :: coef,y
    real(kind = rk), dimension(r_dim) :: rs
    real(kind = rk) :: u0,det_a,a1,a2,a3,a4,a5
    integer :: i

    rs = r/r(r_dim)
    a = 0.0_rk
    ainv = 0.0_rk
    a1 = sum(rs**4)
    a2 = sum(rs**3)
    a3 = sum(rs**2)
    a4 = sum(rs)
    a5 = r_dim

    a(1,1) = a1
    a(1,2) = a2
    a(1,3) = a3
    a(2,1) = a2
    a(2,2) = a3
    a(2,3) = a4
    a(3,1) = a3
    a(3,2) = a4
    a(3,3) = a5
    
    y(1) = dot_product(rs**2,u)
    y(2) = dot_product(rs,u)
    y(3) = sum(u)

    ainv(1,1) = a(3,3)*a(2,2)-a(3,2)*a(2,3)
    ainv(1,2) = -(a(3,3)*a(1,2)-a(3,2)*a(1,3))
    ainv(1,3) = a(2,3)*a(1,2)-a(2,2)*a(1,3)

    ainv(2,1) = -(a(3,3)*a(2,1)-a(3,1)*a(2,3))
    ainv(2,2) = a(3,3)*a(1,1)-a(3,1)*a(1,3)
    ainv(2,3) = -(a(2,3)*a(1,1)-a(2,1)*a(1,3))

    ainv(3,1) = a(3,2)*a(2,1)-a(3,1)*a(2,2)
    ainv(3,2) = -(a(3,2)*a(1,1)-a(3,1)*a(1,2))
    ainv(3,3) = a(2,2)*a(1,1)-a(2,1)*a(1,2)

    det_a = a(1,1)*(a(3,3)*a(2,2)-a(3,2)*a(2,3)) &
         -a(2,1)*(a(3,3)*a(1,2)-a(3,2)*a(1,3)) &
         +a(3,1)*(a(2,3)*a(1,2)-a(2,2)*a(1,3))

    ainv = ainv/det_a

    coef = matmul(ainv,y)
    
    u0 = coef(1)*(r0/r(r_dim))**2+coef(2)*r0/r(r_dim)+coef(3)    

  end function extrapolatePP




  function evalLR(x,pot,spline_indx) result(value)
    implicit none
    integer :: spline_indx,pot
    real(kind = rk) :: x,y
    real(kind = rk) :: u,v
    real(kind = rk) :: p1u,p2u,q1u,q2u
    real(kind = rk) :: value
    real(kind = rk), dimension(4,4) :: CSmatrix
    real(kind = rk), dimension(4) :: pqu,pqv
    integer :: i,j
    
    !i = get_i(LRcspline(spline_indx)%x,LRcspline(spline_indx)%x_dim,x)
    
    !uniform grid (should extend over the box limits)
    i = floor(x/SRcorrection(pot)%LRcspline(spline_indx)%hx(1)) + 1
    if (i>=SRcorrection(pot)%LRcspline(spline_indx)%x_dim) then
       i = i-1
    end if

    u = (x-SRcorrection(pot)%LRcspline(spline_indx)%x(i))/SRcorrection(pot)%LRcspline(spline_indx)%hx(i)
    
    p1u = (1.0_rk+2.0_rk*u)*(u-1.0_rk)**2
    
    p2u = u**2*(3.0_rk-2.0_rk*u)
    
    q1u = u*(u-1.0_rk)**2
    
    q2u = u**2*(u-1.0_rk)

    
    value = SRcorrection(pot)%LRcspline(spline_indx)%F(i)*p1u &
         +SRcorrection(pot)%LRcspline(spline_indx)%F(i+1)*p2u &
         +SRcorrection(pot)%LRcspline(spline_indx)%hx(i)*&
         (SRcorrection(pot)%LRcspline(spline_indx)%Fx(i)*q1u &
         +SRcorrection(pot)%LRcspline(spline_indx)%Fx(i+1)*q2u)    
   
    
  end function evalLR


  function evaldLR(x,pot,spline_indx) result(value)
    implicit none
    integer :: spline_indx,pot
    real(kind = rk) :: x,y
    real(kind = rk) :: u,v
    real(kind = rk) :: p1u,p2u,q1u,q2u
    real(kind = rk) :: value
    real(kind = rk), dimension(4,4) :: CSmatrix
    real(kind = rk), dimension(4) :: pqu,pqv
    integer :: i,j
    
    !i = get_i(LRcspline(spline_indx)%x,LRcspline(spline_indx)%x_dim,x)
    
    !uniform grid (should extend over the box limits)
    i = floor(x/SRcorrection(pot)%LRcspline(spline_indx)%hx(1)) + 1
    if (i>=SRcorrection(pot)%LRcspline(spline_indx)%x_dim) then
       i = i-1
    end if

    u = (x-SRcorrection(pot)%LRcspline(spline_indx)%x(i))/SRcorrection(pot)%LRcspline(spline_indx)%hx(i)
    
    p1u = dp1(u)
    
    p2u = dp2(u)
    
    q1u = dq1(u)
    
    q2u = dq2(u)
 
    value = (SRcorrection(pot)%LRcspline(spline_indx)%F(i)*p1u &
         +SRcorrection(pot)%LRcspline(spline_indx)%F(i+1)*p2u)/&
         SRcorrection(pot)%LRcspline(spline_indx)%hx(i) &
         +SRcorrection(pot)%LRcspline(spline_indx)%Fx(i)*q1u &
         +SRcorrection(pot)%LRcspline(spline_indx)%Fx(i+1)*q2u
   
    
  end function evaldLR


  function evald2LR(x,pot,spline_indx) result(value)
    implicit none
    integer :: spline_indx,pot
    real(kind = rk) :: x,y
    real(kind = rk) :: u,v
    real(kind = rk) :: p1u,p2u,q1u,q2u
    real(kind = rk) :: value
    real(kind = rk), dimension(4,4) :: CSmatrix
    real(kind = rk), dimension(4) :: pqu,pqv
    integer :: i,j
    
    !i = get_i(LRcspline(spline_indx)%x,LRcspline(spline_indx)%x_dim,x)
    
    !uniform grid (should extend over the box limits)
    i = floor(x/SRcorrection(pot)%LRcspline(spline_indx)%hx(1)) + 1
    if (i>=SRcorrection(pot)%LRcspline(spline_indx)%x_dim) then
       i = i-1
    end if

    u = (x-SRcorrection(pot)%LRcspline(spline_indx)%x(i))/SRcorrection(pot)%LRcspline(spline_indx)%hx(i)
    
    p1u = d2p1(u)
    
    p2u = d2p2(u)
    
    q1u = d2q1(u)
    
    q2u = d2q2(u)
 
    value = (SRcorrection(pot)%LRcspline(spline_indx)%F(i)*p1u &
         +SRcorrection(pot)%LRcspline(spline_indx)%F(i+1)*p2u)/&
         SRcorrection(pot)%LRcspline(spline_indx)%hx(i)**2 &
         +(SRcorrection(pot)%LRcspline(spline_indx)%Fx(i)*q1u &
         +SRcorrection(pot)%LRcspline(spline_indx)%Fx(i+1)*q2u)/&
         SRcorrection(pot)%LRcspline(spline_indx)%hx(i)
   
    
  end function evald2LR




  function eval2dcspline(x,y,spline_indx) result(value)
    implicit none
    integer :: spline_indx
    real(kind = rk) :: x,y
    real(kind = rk) :: u,v
    real(kind = rk) :: p1v,p2v,q1v,q2v
    real(kind = rk) :: p1u,p2u,q1u,q2u
    real(kind = rk) :: value
    real(kind = rk), dimension(4,4) :: CSmatrix
    real(kind = rk), dimension(4) :: pqu,pqv
    integer :: i,j
    
    i = get_i(Ucspline(spline_indx)%x,Ucspline(spline_indx)%x_dim,x)
    if (i==Ucspline(spline_indx)%x_dim) then
       i = i-1
    end if

    j = get_i(Ucspline(spline_indx)%y,Ucspline(spline_indx)%y_dim,y)
    if (j==Ucspline(spline_indx)%y_dim) then
       j = j-1
    end if

    u = (x-Ucspline(spline_indx)%x(i))/Ucspline(spline_indx)%hx(i)
    v = (y-Ucspline(spline_indx)%y(j))/Ucspline(spline_indx)%hy(j)
    
    p1u = (1.0_rk+2.0_rk*u)*(u-1.0_rk)**2
    p1v = (1.0_rk+2.0_rk*v)*(v-1.0_rk)**2
    
    p2u = u**2*(3.0_rk-2.0_rk*u)
    p2v = v**2*(3.0_rk-2.0_rk*v)
    
    q1u = u*(u-1.0_rk)**2
    q1v = v*(v-1.0_rk)**2
    
    q2u = u**2*(u-1.0_rk)
    q2v = v**2*(v-1.0_rk)

    CSmatrix(1,1) = Ucspline(spline_indx)%F(i,j)
    CSmatrix(1,2) = Ucspline(spline_indx)%F(i,j+1)
    CSmatrix(1,3) = Ucspline(spline_indx)%Fy(i,j)
    CSmatrix(1,4) = Ucspline(spline_indx)%Fy(i,j+1)

    CSmatrix(2,1) = Ucspline(spline_indx)%F(i+1,j)
    CSmatrix(2,2) = Ucspline(spline_indx)%F(i+1,j+1)
    CSmatrix(2,3) = Ucspline(spline_indx)%Fy(i+1,j)
    CSmatrix(2,4) = Ucspline(spline_indx)%Fy(i+1,j+1)
    
    CSmatrix(3,1) = Ucspline(spline_indx)%Fx(i,j)
    CSmatrix(3,2) = Ucspline(spline_indx)%Fx(i,j+1)
    CSmatrix(3,3) = Ucspline(spline_indx)%Fxy(i,j)
    CSmatrix(3,4) = Ucspline(spline_indx)%Fxy(i,j+1)
    
    CSmatrix(4,1) = Ucspline(spline_indx)%Fx(i+1,j)
    CSmatrix(4,2) = Ucspline(spline_indx)%Fx(i+1,j+1)
    CSmatrix(4,3) = Ucspline(spline_indx)%Fxy(i+1,j)
    CSmatrix(4,4) = Ucspline(spline_indx)%Fxy(i+1,j+1)

    pqu = (/p1u, p2u, Ucspline(spline_indx)%hx(i)*q1u, &
         Ucspline(spline_indx)%hx(i)*q2u/)

    pqv = (/p1v, p2v, Ucspline(spline_indx)%hy(j)*q1v, &
         Ucspline(spline_indx)%hy(j)*q2v/)


    

    value = 0.0_rk
    do i=1,4
       !write(*,'(4(F10.5))') CSmatrix(i,1:4)       
       value = value + dot_product(CSmatrix(i,1:4),pqv)*pqu(i)
    end do
    !pause


    !value = dot_product(pqu,matmul(CSmatrix,pqv))

    
  end function eval2dcspline


  function eval2dcspline2(x,y,spline_indx) result(value)
    implicit none
    integer :: spline_indx
    real(kind = rk) :: x,y
    real(kind = rk) :: u,v
    real(kind = rk) :: p1v,p2v,q1v,q2v
    real(kind = rk) :: p1u,p2u,q1u,q2u
    real(kind = rk) :: value
    real(kind = rk), dimension(4,4) :: CSmatrix
    real(kind = rk), dimension(4) :: pqu,pqv
    integer :: i,j
    
    i = get_i(Ucspline2(spline_indx)%x,Ucspline2(spline_indx)%x_dim,x)
    if (i==Ucspline2(spline_indx)%x_dim) then
       i = i-1
    end if

    j = get_i(Ucspline2(spline_indx)%y,Ucspline2(spline_indx)%y_dim,y)
    if (j==Ucspline2(spline_indx)%y_dim) then
       j = j-1
    end if

    u = (x-Ucspline2(spline_indx)%x(i))/Ucspline2(spline_indx)%hx(i)
    v = (y-Ucspline2(spline_indx)%y(j))/Ucspline2(spline_indx)%hy(j)
    
    p1u = (1.0_rk+2.0_rk*u)*(u-1.0_rk)**2
    p1v = (1.0_rk+2.0_rk*v)*(v-1.0_rk)**2
    
    p2u = u**2*(3.0_rk-2.0_rk*u)
    p2v = v**2*(3.0_rk-2.0_rk*v)
    
    q1u = u*(u-1.0_rk)**2
    q1v = v*(v-1.0_rk)**2
    
    q2u = u**2*(u-1.0_rk)
    q2v = v**2*(v-1.0_rk)

    CSmatrix(1,1) = Ucspline2(spline_indx)%F(i,j)
    CSmatrix(1,2) = Ucspline2(spline_indx)%F(i,j+1)
    CSmatrix(1,3) = Ucspline2(spline_indx)%Fy(i,j)
    CSmatrix(1,4) = Ucspline2(spline_indx)%Fy(i,j+1)

    CSmatrix(2,1) = Ucspline2(spline_indx)%F(i+1,j)
    CSmatrix(2,2) = Ucspline2(spline_indx)%F(i+1,j+1)
    CSmatrix(2,3) = Ucspline2(spline_indx)%Fy(i+1,j)
    CSmatrix(2,4) = Ucspline2(spline_indx)%Fy(i+1,j+1)
    
    CSmatrix(3,1) = Ucspline2(spline_indx)%Fx(i,j)
    CSmatrix(3,2) = Ucspline2(spline_indx)%Fx(i,j+1)
    CSmatrix(3,3) = Ucspline2(spline_indx)%Fxy(i,j)
    CSmatrix(3,4) = Ucspline2(spline_indx)%Fxy(i,j+1)
    
    CSmatrix(4,1) = Ucspline2(spline_indx)%Fx(i+1,j)
    CSmatrix(4,2) = Ucspline2(spline_indx)%Fx(i+1,j+1)
    CSmatrix(4,3) = Ucspline2(spline_indx)%Fxy(i+1,j)
    CSmatrix(4,4) = Ucspline2(spline_indx)%Fxy(i+1,j+1)

    pqu = (/p1u, p2u, Ucspline2(spline_indx)%hx(i)*q1u, &
         Ucspline2(spline_indx)%hx(i)*q2u/)

    pqv = (/p1v, p2v, Ucspline2(spline_indx)%hy(j)*q1v, &
         Ucspline2(spline_indx)%hy(j)*q2v/)


    

    value = 0.0_rk
    do i=1,4
       !write(*,'(4(F10.5))') CSmatrix(i,1:4)       
       value = value + dot_product(CSmatrix(i,1:4),pqv)*pqu(i)
    end do
    !pause


    !value = dot_product(pqu,matmul(CSmatrix,pqv))

    
  end function eval2dcspline2


  function eval2d_dU_dbeta(x,y,spline_indx) result(value)
    implicit none
    integer :: spline_indx
    real(kind = rk) :: x,y
    real(kind = rk) :: u,v
    real(kind = rk) :: p1v,p2v,q1v,q2v
    real(kind = rk) :: p1u,p2u,q1u,q2u
    real(kind = rk) :: value
    real(kind = rk), dimension(4,4) :: CSmatrix
    real(kind = rk), dimension(4) :: pqu,pqv
    integer :: i,j
    
    i = get_i(dU_dbeta(spline_indx)%x,dU_dbeta(spline_indx)%x_dim,x)
    if (i==dU_dbeta(spline_indx)%x_dim) then
       i = i-1
    end if

    j = get_i(dU_dbeta(spline_indx)%y,dU_dbeta(spline_indx)%y_dim,y)
    if (j==dU_dbeta(spline_indx)%y_dim) then
       j = j-1
    end if

    u = (x-dU_dbeta(spline_indx)%x(i))/dU_dbeta(spline_indx)%hx(i)
    v = (y-dU_dbeta(spline_indx)%y(j))/dU_dbeta(spline_indx)%hy(j)
    
    p1u = (1.0_rk+2.0_rk*u)*(u-1.0_rk)**2
    p1v = (1.0_rk+2.0_rk*v)*(v-1.0_rk)**2
    
    p2u = u**2*(3.0_rk-2.0_rk*u)
    p2v = v**2*(3.0_rk-2.0_rk*v)
    
    q1u = u*(u-1.0_rk)**2
    q1v = v*(v-1.0_rk)**2
    
    q2u = u**2*(u-1.0_rk)
    q2v = v**2*(v-1.0_rk)

    CSmatrix(1,1) = dU_dbeta(spline_indx)%F(i,j)
    CSmatrix(1,2) = dU_dbeta(spline_indx)%F(i,j+1)
    CSmatrix(1,3) = dU_dbeta(spline_indx)%Fy(i,j)
    CSmatrix(1,4) = dU_dbeta(spline_indx)%Fy(i,j+1)

    CSmatrix(2,1) = dU_dbeta(spline_indx)%F(i+1,j)
    CSmatrix(2,2) = dU_dbeta(spline_indx)%F(i+1,j+1)
    CSmatrix(2,3) = dU_dbeta(spline_indx)%Fy(i+1,j)
    CSmatrix(2,4) = dU_dbeta(spline_indx)%Fy(i+1,j+1)
    
    CSmatrix(3,1) = dU_dbeta(spline_indx)%Fx(i,j)
    CSmatrix(3,2) = dU_dbeta(spline_indx)%Fx(i,j+1)
    CSmatrix(3,3) = dU_dbeta(spline_indx)%Fxy(i,j)
    CSmatrix(3,4) = dU_dbeta(spline_indx)%Fxy(i,j+1)
    
    CSmatrix(4,1) = dU_dbeta(spline_indx)%Fx(i+1,j)
    CSmatrix(4,2) = dU_dbeta(spline_indx)%Fx(i+1,j+1)
    CSmatrix(4,3) = dU_dbeta(spline_indx)%Fxy(i+1,j)
    CSmatrix(4,4) = dU_dbeta(spline_indx)%Fxy(i+1,j+1)

    pqu = (/p1u, p2u, dU_dbeta(spline_indx)%hx(i)*q1u, &
         dU_dbeta(spline_indx)%hx(i)*q2u/)

    pqv = (/p1v, p2v, dU_dbeta(spline_indx)%hy(j)*q1v, &
         dU_dbeta(spline_indx)%hy(j)*q2v/)


    

    value = 0.0_rk
    do i=1,4
       !write(*,'(4(F10.5))') CSmatrix(i,1:4)       
       value = value + dot_product(CSmatrix(i,1:4),pqv)*pqu(i)
    end do
    !pause


    !value = dot_product(pqu,matmul(CSmatrix,pqv))

    
  end function eval2d_dU_dbeta



  
  function Uij(i,j,spline_indx) result(value)
    implicit none
    integer :: i,j,spline_indx
    real(kind = rk) :: value
    
    value = Ucspline(spline_indx)%F(i,j)
    
  end function Uij


  function dU_dBetaij(i,j,spline_indx) result(value)
    implicit none
    integer :: i,j,spline_indx
    real(kind = rk) :: value

    value = dU_dbeta(spline_indx)%F(i,j)

  end function dU_dBetaij




  function dUx(rx,ry,spline_indx) result(value)
    implicit none
    integer :: i,j,spline_indx
    real(kind = rk) :: value,x,y,u,v,rx,ry
    real(kind = rk), dimension(4) :: a,b 
    real(kind = rk), dimension(4,4) :: CSmatrix

    x = rx
    y = ry
    
    i = get_i(Ucspline(spline_indx)%x,Ucspline(spline_indx)%x_dim,x)
    if (i==Ucspline(spline_indx)%x_dim) then
       i = i-1
    end if
    
    j = get_i(Ucspline(spline_indx)%y,Ucspline(spline_indx)%y_dim,y)
    if (j==Ucspline(spline_indx)%y_dim) then
       j = j-1
    end if

    u = (x-Ucspline(spline_indx)%x(i))/Ucspline(spline_indx)%hx(i)
    v = (y-Ucspline(spline_indx)%y(j))/Ucspline(spline_indx)%hy(j)
    
    a(1) = dp1(u)/Ucspline(spline_indx)%hx(i)
    a(2) = dp2(u)/Ucspline(spline_indx)%hx(i)
    a(3) = dq1(u)
    a(4) = dq2(u)

    b(1) = p1(v)
    b(2) = p2(v)
    b(3) = q1(v)*Ucspline(spline_indx)%hy(j)
    b(4) = q2(v)*Ucspline(spline_indx)%hy(j)
    
    CSmatrix(1,1) = Ucspline(spline_indx)%F(i,j)
    CSmatrix(1,2) = Ucspline(spline_indx)%F(i,j+1)
    CSmatrix(1,3) = Ucspline(spline_indx)%Fy(i,j)
    CSmatrix(1,4) = Ucspline(spline_indx)%Fy(i,j+1)
    
    CSmatrix(2,1) = Ucspline(spline_indx)%F(i+1,j)
    CSmatrix(2,2) = Ucspline(spline_indx)%F(i+1,j+1)
    CSmatrix(2,3) = Ucspline(spline_indx)%Fy(i+1,j)
    CSmatrix(2,4) = Ucspline(spline_indx)%Fy(i+1,j+1)
    
    CSmatrix(3,1) = Ucspline(spline_indx)%Fx(i,j)
    CSmatrix(3,2) = Ucspline(spline_indx)%Fx(i,j+1)
    CSmatrix(3,3) = Ucspline(spline_indx)%Fxy(i,j)
    CSmatrix(3,4) = Ucspline(spline_indx)%Fxy(i,j+1)
    
    CSmatrix(4,1) = Ucspline(spline_indx)%Fx(i+1,j)
    CSmatrix(4,2) = Ucspline(spline_indx)%Fx(i+1,j+1)
    CSmatrix(4,3) = Ucspline(spline_indx)%Fxy(i+1,j)
    CSmatrix(4,4) = Ucspline(spline_indx)%Fxy(i+1,j+1)


    value = 0.0_rk
    do i=1,4
       !write(*,'(4(F10.5))') CSmatrix(i,1:4)       
       value = value + dot_product(CSmatrix(i,1:4),b)*a(i)
    end do

  end function dUx


  function dUy(rx,ry,spline_indx) result(value)
    implicit none
    integer :: i,j,spline_indx
    real(kind = rk) :: value,x,y,u,v,rx,ry
    real(kind = rk), dimension(4) :: a,b 
    real(kind = rk), dimension(4,4) :: CSmatrix

    x = rx
    y = ry


    i = get_i(Ucspline(spline_indx)%x,Ucspline(spline_indx)%x_dim,x)
    if (i==Ucspline(spline_indx)%x_dim) then
       i = i-1
    end if
    
    j = get_i(Ucspline(spline_indx)%y,Ucspline(spline_indx)%y_dim,y)
    if (j==Ucspline(spline_indx)%y_dim) then
       j = j-1
    end if

    u = (x-Ucspline(spline_indx)%x(i))/Ucspline(spline_indx)%hx(i)
    v = (y-Ucspline(spline_indx)%y(j))/Ucspline(spline_indx)%hy(j)
    
    a(1) = p1(u)
    a(2) = p2(u)
    a(3) = q1(u)*Ucspline(spline_indx)%hx(i)
    a(4) = q2(u)*Ucspline(spline_indx)%hx(i)

    b(1) = dp1(v)/Ucspline(spline_indx)%hy(j)
    b(2) = dp2(v)/Ucspline(spline_indx)%hy(j)
    b(3) = dq1(v)
    b(4) = dq2(v)
    
    CSmatrix(1,1) = Ucspline(spline_indx)%F(i,j)
    CSmatrix(1,2) = Ucspline(spline_indx)%F(i,j+1)
    CSmatrix(1,3) = Ucspline(spline_indx)%Fy(i,j)
    CSmatrix(1,4) = Ucspline(spline_indx)%Fy(i,j+1)
    
    CSmatrix(2,1) = Ucspline(spline_indx)%F(i+1,j)
    CSmatrix(2,2) = Ucspline(spline_indx)%F(i+1,j+1)
    CSmatrix(2,3) = Ucspline(spline_indx)%Fy(i+1,j)
    CSmatrix(2,4) = Ucspline(spline_indx)%Fy(i+1,j+1)
    
    CSmatrix(3,1) = Ucspline(spline_indx)%Fx(i,j)
    CSmatrix(3,2) = Ucspline(spline_indx)%Fx(i,j+1)
    CSmatrix(3,3) = Ucspline(spline_indx)%Fxy(i,j)
    CSmatrix(3,4) = Ucspline(spline_indx)%Fxy(i,j+1)
    
    CSmatrix(4,1) = Ucspline(spline_indx)%Fx(i+1,j)
    CSmatrix(4,2) = Ucspline(spline_indx)%Fx(i+1,j+1)
    CSmatrix(4,3) = Ucspline(spline_indx)%Fxy(i+1,j)
    CSmatrix(4,4) = Ucspline(spline_indx)%Fxy(i+1,j+1)


    value = 0.0_rk
    do i=1,4
       !write(*,'(4(F10.5))') CSmatrix(i,1:4)       
       value = value + dot_product(CSmatrix(i,1:4),b)*a(i)
    end do
    

  end function dUy



  function d2Ux(rx,ry,spline_indx) result(value)
    implicit none
    integer :: i,j,spline_indx
    real(kind = rk) :: value,x,y,u,v,rx,ry
    real(kind = rk), dimension(4) :: a,b 
    real(kind = rk), dimension(4,4) :: CSmatrix

    x = rx
    y = ry
    
    i = get_i(Ucspline(spline_indx)%x,Ucspline(spline_indx)%x_dim,x)
    if (i==Ucspline(spline_indx)%x_dim) then
       i = i-1
    end if
    
    j = get_i(Ucspline(spline_indx)%y,Ucspline(spline_indx)%y_dim,y)
    if (j==Ucspline(spline_indx)%y_dim) then
       j = j-1
    end if

    u = (x-Ucspline(spline_indx)%x(i))/Ucspline(spline_indx)%hx(i)
    v = (y-Ucspline(spline_indx)%y(j))/Ucspline(spline_indx)%hy(j)
    
    a(1) = d2p1(u)/Ucspline(spline_indx)%hx(i)**2
    a(2) = d2p2(u)/Ucspline(spline_indx)%hx(i)**2
    a(3) = d2q1(u)/Ucspline(spline_indx)%hx(i)
    a(4) = d2q2(u)/Ucspline(spline_indx)%hx(i)
    
    b(1) = p1(v)
    b(2) = p2(v)
    b(3) = q1(v)*Ucspline(spline_indx)%hy(j)
    b(4) = q2(v)*Ucspline(spline_indx)%hy(j)
    
    CSmatrix(1,1) = Ucspline(spline_indx)%F(i,j)
    CSmatrix(1,2) = Ucspline(spline_indx)%F(i,j+1)
    CSmatrix(1,3) = Ucspline(spline_indx)%Fy(i,j)
    CSmatrix(1,4) = Ucspline(spline_indx)%Fy(i,j+1)
    
    CSmatrix(2,1) = Ucspline(spline_indx)%F(i+1,j)
    CSmatrix(2,2) = Ucspline(spline_indx)%F(i+1,j+1)
    CSmatrix(2,3) = Ucspline(spline_indx)%Fy(i+1,j)
    CSmatrix(2,4) = Ucspline(spline_indx)%Fy(i+1,j+1)
    
    CSmatrix(3,1) = Ucspline(spline_indx)%Fx(i,j)
    CSmatrix(3,2) = Ucspline(spline_indx)%Fx(i,j+1)
    CSmatrix(3,3) = Ucspline(spline_indx)%Fxy(i,j)
    CSmatrix(3,4) = Ucspline(spline_indx)%Fxy(i,j+1)
    
    CSmatrix(4,1) = Ucspline(spline_indx)%Fx(i+1,j)
    CSmatrix(4,2) = Ucspline(spline_indx)%Fx(i+1,j+1)
    CSmatrix(4,3) = Ucspline(spline_indx)%Fxy(i+1,j)
    CSmatrix(4,4) = Ucspline(spline_indx)%Fxy(i+1,j+1)


    value = 0.0_rk
    do i=1,4
       !write(*,'(4(F10.5))') CSmatrix(i,1:4)       
       value = value + dot_product(CSmatrix(i,1:4),b)*a(i)
    end do


  end function d2Ux


  function d2Uy(rx,ry,spline_indx) result(value)
    implicit none
    integer :: i,j,spline_indx
    real(kind = rk) :: value,x,y,u,v,rx,ry
    real(kind = rk), dimension(4) :: a,b 
    real(kind = rk), dimension(4,4) :: CSmatrix

    x = rx
    y = ry


    i = get_i(Ucspline(spline_indx)%x,Ucspline(spline_indx)%x_dim,x)
    if (i==Ucspline(spline_indx)%x_dim) then
       i = i-1
    end if
    
    j = get_i(Ucspline(spline_indx)%y,Ucspline(spline_indx)%y_dim,y)
    if (j==Ucspline(spline_indx)%y_dim) then
       j = j-1
    end if

    u = (x-Ucspline(spline_indx)%x(i))/Ucspline(spline_indx)%hx(i)
    v = (y-Ucspline(spline_indx)%y(j))/Ucspline(spline_indx)%hy(j)
    
    a(1) = p1(u)
    a(2) = p2(u)
    a(3) = q1(u)*Ucspline(spline_indx)%hx(i)
    a(4) = q2(u)*Ucspline(spline_indx)%hx(i)
    
    b(1) = d2p1(v)/Ucspline(spline_indx)%hy(j)**2
    b(2) = d2p2(v)/Ucspline(spline_indx)%hy(j)**2
    b(3) = d2q1(v)/Ucspline(spline_indx)%hy(j)
    b(4) = d2q2(v)/Ucspline(spline_indx)%hy(j)
    
    CSmatrix(1,1) = Ucspline(spline_indx)%F(i,j)
    CSmatrix(1,2) = Ucspline(spline_indx)%F(i,j+1)
    CSmatrix(1,3) = Ucspline(spline_indx)%Fy(i,j)
    CSmatrix(1,4) = Ucspline(spline_indx)%Fy(i,j+1)
    
    CSmatrix(2,1) = Ucspline(spline_indx)%F(i+1,j)
    CSmatrix(2,2) = Ucspline(spline_indx)%F(i+1,j+1)
    CSmatrix(2,3) = Ucspline(spline_indx)%Fy(i+1,j)
    CSmatrix(2,4) = Ucspline(spline_indx)%Fy(i+1,j+1)
    
    CSmatrix(3,1) = Ucspline(spline_indx)%Fx(i,j)
    CSmatrix(3,2) = Ucspline(spline_indx)%Fx(i,j+1)
    CSmatrix(3,3) = Ucspline(spline_indx)%Fxy(i,j)
    CSmatrix(3,4) = Ucspline(spline_indx)%Fxy(i,j+1)
    
    CSmatrix(4,1) = Ucspline(spline_indx)%Fx(i+1,j)
    CSmatrix(4,2) = Ucspline(spline_indx)%Fx(i+1,j+1)
    CSmatrix(4,3) = Ucspline(spline_indx)%Fxy(i+1,j)
    CSmatrix(4,4) = Ucspline(spline_indx)%Fxy(i+1,j+1)


    value = 0.0_rk
    do i=1,4
       !write(*,'(4(F10.5))') CSmatrix(i,1:4)       
       value = value + dot_product(CSmatrix(i,1:4),b)*a(i)
    end do

  end function d2Uy



  function d2Uxy(rx,ry,spline_indx) result(value)
    implicit none
    integer :: i,j,spline_indx
    real(kind = rk) :: value,x,y,u,v,rx,ry
    real(kind = rk), dimension(4) :: a,b 
    real(kind = rk), dimension(4,4) :: CSmatrix

    x = rx
    y = ry
    
    i = get_i(Ucspline(spline_indx)%x,Ucspline(spline_indx)%x_dim,x)
    if (i==Ucspline(spline_indx)%x_dim) then
       i = i-1
    end if
    
    j = get_i(Ucspline(spline_indx)%y,Ucspline(spline_indx)%y_dim,y)
    if (j==Ucspline(spline_indx)%y_dim) then
       j = j-1
    end if

    u = (x-Ucspline(spline_indx)%x(i))/Ucspline(spline_indx)%hx(i)
    v = (y-Ucspline(spline_indx)%y(j))/Ucspline(spline_indx)%hy(j)
    
    a(1) = dp1(u)/Ucspline(spline_indx)%hx(i)
    a(2) = dp2(u)/Ucspline(spline_indx)%hx(i)
    a(3) = dq1(u)
    a(4) = dq2(u)

    b(1) = dp1(v)/Ucspline(spline_indx)%hy(j)
    b(2) = dp2(v)/Ucspline(spline_indx)%hy(j)
    b(3) = dq1(v)
    b(4) = dq2(v)
    
    CSmatrix(1,1) = Ucspline(spline_indx)%F(i,j)
    CSmatrix(1,2) = Ucspline(spline_indx)%F(i,j+1)
    CSmatrix(1,3) = Ucspline(spline_indx)%Fy(i,j)
    CSmatrix(1,4) = Ucspline(spline_indx)%Fy(i,j+1)
    
    CSmatrix(2,1) = Ucspline(spline_indx)%F(i+1,j)
    CSmatrix(2,2) = Ucspline(spline_indx)%F(i+1,j+1)
    CSmatrix(2,3) = Ucspline(spline_indx)%Fy(i+1,j)
    CSmatrix(2,4) = Ucspline(spline_indx)%Fy(i+1,j+1)
    
    CSmatrix(3,1) = Ucspline(spline_indx)%Fx(i,j)
    CSmatrix(3,2) = Ucspline(spline_indx)%Fx(i,j+1)
    CSmatrix(3,3) = Ucspline(spline_indx)%Fxy(i,j)
    CSmatrix(3,4) = Ucspline(spline_indx)%Fxy(i,j+1)
    
    CSmatrix(4,1) = Ucspline(spline_indx)%Fx(i+1,j)
    CSmatrix(4,2) = Ucspline(spline_indx)%Fx(i+1,j+1)
    CSmatrix(4,3) = Ucspline(spline_indx)%Fxy(i+1,j)
    CSmatrix(4,4) = Ucspline(spline_indx)%Fxy(i+1,j+1)


    value = 0.0_rk
    do i=1,4
       !write(*,'(4(F10.5))') CSmatrix(i,1:4)       
       value = value + dot_product(CSmatrix(i,1:4),b)*a(i)
    end do

  end function d2Uxy


  function get_i(x,x_dim,r) result(ind_x)
    implicit none
    integer :: x_dim,i
    real(kind =rk),dimension(x_dim) ::  x
    integer ::  ind_x
    real(kind =rk) :: r
    
    do i = 1,x_dim-1
       if (r>=x(i) .and. r<=x(i+1)) then
          ind_x = i
          return
       end if
    end do
    ind_x = x_dim

  end function get_i



  function p1(t) result(value)
    implicit none
    real(kind = rk) :: t,value
    
    value = (t-1.0_rk)*(t-1.0_rk)*(1.0_rk+2.0_rk*t)
    
  end function p1

  
  function p2(t) result(value)
    implicit none
    real(kind = rk) :: t,value

    value = t*t*(3.0_rk-2.0_rk*t)

  end function p2

  function q1(t) result(value)
    implicit none
    real(kind = rk) :: t,value

    value = t*(t-1.0_rk)*(t-1.0_rk)

  end function q1

  function q2(t) result(value)
    implicit none
    real(kind = rk) :: t,value

    value = t*t*(t-1.0_rk)

  end function q2

  function dp1(t) result(value)
    implicit none
    real(kind = rk) :: t,value

    value = 6.0_rk*t*(t-1.0_rk)

  end function dp1

  function dq1(t) result(value)
    implicit none
    real(kind = rk) :: t,value

    value = (t-1.0_rk)*(3.0_rk*t-1.0_rk)

  end function dq1
  
  function dp2(t) result(value)
    implicit none
    real(kind = rk) :: t,value
    
    value = -dp1(t)
    
  end function dp2

  function dq2(t) result(value)
    implicit none
    real(kind = rk) :: t,value

    value = 3.0_rk*t*t - 2.0_rk*t

  end function dq2


  function d2p1(t) result(value)
    implicit none
    real(kind = rk) :: t,value

    value = 12.0_rk*t-6.0_rk

  end function d2p1


  function d2q1(t) result(value)
    implicit none
    real(kind = rk) :: t,value

    value = 6.0_rk*t - 4.0_rk

  end function d2q1

  function d2p2(t) result(value)
    implicit none
    real(kind = rk) :: t,value
    
    value = -d2p1(t)

  end function d2p2


  function d2q2(t) result(value)
    implicit none
    real(kind = rk) :: t,value

    value = 6.0_rk*t - 2.0_rk

  end function d2q2


  
  FUNCTION gaussinen() RESULT( Gauss_taul )
    !--------------------------------------------------------
    ! Tama generoi normaalisti jakautuneita satunnaislukuja
    !
    
    !Box-Muller algoritmi
    
    
    REAL(kind = rk), DIMENSION(3)      :: Gauss_taul
    REAL(kind = rk), DIMENSION(2)      :: ksi
    
    CALL RANDOM_NUMBER( ksi )
    
    gauss_taul(1) = SQRT( -2.0_rk*LOG(ksi(1)) )*COS( 2.0_rk*pi*ksi(2) )
    gauss_taul(2) = SQRT( -2.0_rk*LOG(ksi(1)) )*SIN( 2.0_rk*pi*ksi(2) )
    
    !Toinen kerta toden sanoo
    CALL RANDOM_NUMBER( ksi )
    
    gauss_taul(3) = SQRT( -2.0_rk*LOG(ksi(1)) )*COS( 2.0_rk*pi*ksi(2) )
    
  END FUNCTION gaussinen


  
  FUNCTION erfc(x)
    !Numerical Recipes:
    !Returns the complementary error function erfc(x) with
    !fractional error everywhere less than 1.2 function erfcc(x)
    !Numerical Recipes:
    IMPLICIT NONE
    REAL(kind=rk) :: x,z,t,ans,erfc
    z=ABS(x);
    t=1.0_rk/(1.0_rk+z/2);
    ans=t*EXP(-z*z-1.26551223_rk+t*(1.00002368_rk+t*(0.37409196_rk+t*(0.09678418_rk+&
        t*(-0.18628806_rk+t*(0.27886807_rk+t*(-1.13520398_rk+t*(1.48851587_rk+&
        t*(-0.82215223_rk+t*0.17087277_rk)))))))));
    IF (x>=0.0_rk) THEN
       erfc = ans
    ELSE
       erfc = 2.0_rk - ans
    END IF
  END FUNCTION erfc



  FUNCTION erff(t)
    IMPLICIT NONE
    REAL(kind=rk) :: t,erff

    erff = 1.0_rk - erfc(t)

  END FUNCTION erff

    
  ! function log1p(z) result(y)
  !   ! does not include any check-ups, e.g. if z is NaN
  !   real(kind = rk) :: z,y,u
    
  !   u = 1.0_rk+z
  !   y = z
  !   if (u <= 0.0_rk) then
  !      y = log(u)
  !   !else if (u .ne. 1.0_rk) then
  !   else if (abs(u-1.0_rk)> 1.0e-20_rk) then
  !      y = log(u)*(z/(u-1.0_rk))
  !   end if

  ! end function log1p

  
  function log1p(z) result(y)
    ! does not include any check-ups, e.g. if z is NaN
    real(kind = rk) :: z,y,u
    
    if (z<1.0e-10_rk .and. z>-1.0e-10_rk) then
       y = z-z**2/2+z**3/6-z**4/24
    else
       y = log(1.0_rk+z)
    end if

  end function log1p


  function expm1(z) result(y)
    ! does not include any check-ups, e.g. if z is NaN
    real(kind = rk) :: z,y,u
    
    u = exp(z)
    y = z
    !if (u .eq. 0.0_rk) then
    !if (u < 1.0e-20_rk) then
    if (u < 1.0e-40_rk) then
       y = u-1.0_rk
       return
    end if
    !if (u .ne. 1.0_rk) then
    !if (abs(u-1.0_rk)>1.0e-20_rk) then
    if (abs(u-1.0_rk)>1.0e-40_rk) then
       y = (u-1.0_rk)*(z/log(u))
    end if

    
  end function expm1


!   function expm1(z) result(y)
!     ! does not include any check-ups, e.g. if z is NaN
!     real(kind = rk) :: z,y,u
            
!     if (z<1.0e-10_rk .and. z>-1.0e-10_rk) then
!        y = z+z**2/2+z**3/6+z**4/24
!     else
!        y = exp(z)-1.0_rk
!     end if
    
!   end function expm1


  function sinc(z,z_dim) result(y)
    implicit none
    integer, intent(in) :: z_dim
    real(kind = rk), dimension(z_dim), intent(in) :: z
    real(kind = rk), dimension(z_dim) :: y
    integer :: i
    
    y=0.0_rk
    do i=1,z_dim
       ! error term is of order z**6
       if (z(i)<1.0e-6_rk) then
          y(i) = 1.0_rk-z(i)**2/6+z(i)**4/120
       else
          y(i) = sin(z(i))/z(i)
       end if
    end do
       
  end function sinc

  function sinc2(z) result(y)
    implicit none
    real(kind = rk), intent(in):: z
    real(kind = rk) :: y
    integer :: i

    y=0.0_rk

    ! error term is of order z**6                                                           
    if (z<1.0E-4_rk) then
       y = 1.0_rk-z**2/6+z**4/120
    else
       y=sin(z)/z
    end if

  end function sinc2


  function dxsinc(z,z_dim) result(y)
    implicit none
    integer, intent(in) :: z_dim
    real(kind = rk), dimension(z_dim), intent(in) :: z
    real(kind = rk), dimension(z_dim) :: y
    integer :: i

    y=0.0_rk
    do i=1,z_dim
       ! error term is of order z**7                                                        
       if (z(i)<1.0E-4_rk) then
          y(i) = -(z(i)/3-z(i)**3/30+z(i)**5/840)
       else
          y(i) = cos(z(i))/z(i)-sin(z(i))/z(i)**2
       end if
    end do

  end function dxsinc

  function dxsinc2(z) result(y)
    implicit none
    real(kind = rk), intent(in) :: z
    real(kind = rk) :: y
    integer :: i

    y=0.0_rk
    ! error term is of order z**7                                                           
    if (z<1.0E-4_rk) then
       y = -(z/3-z**3/30+z**5/840)
    else
       y = cos(z)/z-sin(z)/z**2
    end if


  end function dxsinc2



  function d2xsinc(z,z_dim) result(y)
    implicit none
    integer, intent(in) :: z_dim
    real(kind = rk), dimension(z_dim), intent(in) :: z
    real(kind = rk), dimension(z_dim) :: y
    real(kind = rk) :: v1,v2
    integer :: i

    y=0.0_rk
    do i=1,z_dim
       ! error term is of order z**6                                                        
       v1 = sinc2(z(i))
       if (z(i)<1.0E-4_rk) then
          y(i) = -v1-2.0_rk*(1.0_rk/3-z(i)**2/30+z(i)**4/840)
       else
          v2 = dxsinc2(z(i))
          y(i) = -v1-2.0_rk/z(i)*v2
       end if
    end do

  end function d2xsinc





  function log22(a) result(log2)
    implicit none
    real(kind = rk) :: a,log2
    
    log2 = log10(a)/log10(2.0_rk)
    
  end function log22

  function m0(z) result(value)
    ! m0 = i_0 * exp(-z), where z>=0 and 
    ! i_0 is modified spherical bessel function
    implicit none
    real(kind = rk) :: z,value

    if (z<0.000001_rk) then
       value = exp(-z)*(1.0_rk+z**2/6+z**4/120)
    elseif (z>1000.0_rk) then
       value = 1.0_rk/z/2
    else
       value = (1.0_rk-exp(-2.0_rk*z))/z/2  
    end if

  end function m0


  function dm0(z) result(value)
    ! dm0 = exp(-z) d/dz i_0, where z>=0 and 
    ! i_0 is modified spherical bessel function
    implicit none
    real(kind = rk) :: z,value

    if (z<0.000001_rk) then
       value = exp(-z)*(z/3+z**3/30)
    elseif (z>1000.0_rk) then
       value = 1.0_rk/z/2-1.0_rk/z**2/2
    else
       value = ((z-1.0_rk)+exp(-2.0_rk*z)*(z+1.0_rk))/z**2/2
    end if
    
  end function dm0


  
  function linspace(x,y,dim) result(r)
    implicit none
    integer :: dim, i
    real(kind = rk) :: x,y, dr
    real(kind = rk), dimension(dim) :: r

    
    dr = (y-x)/(dim-1)
    
    do i = 1,dim
       r(i) = x + (i-1)*dr
    end do
    
  end function linspace

  function gridIK2(x,dim) result(r)
    implicit none
    integer :: dim, i
    real(kind = rk) :: x, dr
    real(kind = rk), dimension(dim) :: r,a

    !dr = 0.1_rk
    !dr = 0.2_rk
    dr=0.2_rk
    !dr = 0.05_rk
    a = linspace(-10.0_rk,10.0_rk,dim)
    !a = linspace(-7.0_rk,7.0_rk,dim)
    r(1) = x
    do i = 2,dim
       r(i) = r(i-1) + (-1.0_rk/(exp(a(i))+1.0_rk)+1.0_rk)*dr
    end do
    
  end function gridIK2

  function gridIK2_2(x,y,dim) result(r)
    implicit none
    integer :: dim, i
    real(kind = rk) :: x,y,dx_min,dx_max,a_start,a_end,normi
    real(kind = rk), dimension(dim) :: r
    real(kind = rk), dimension(dim-1) :: a,w
    
    ! preliminary values (will change with respect to the grid size)
    dx_max = 0.2_rk
    dx_min = 1.0e-4_rk
    
    a_start = log( 1.0_rk/(1.0_rk-dx_min/dx_max) - 1.0_rk )
    a_end = -a_start
    
    a = linspace(a_start,a_end,dim-1)
    do i = 1,dim-1
       w(i) = 1.0_rk-1.0_rk/(exp(a(i))+1.0_rk)
    end do
    normi = sum(w)*dx_max/abs(y-x)
    w = w/normi
    
    r(1) = x
    do i = 2,dim
       r(i) = r(i-1) + w(i-1)*dx_max
    end do
    
  end function gridIK2_2

  function gridIK(x,y,dim) result(r)
    implicit none
    integer :: dim, i
    real(kind = rk) :: x,y,dx_max,dx_min,norm
    real(kind = rk), dimension(dim) :: r,coef

    dx_max = 0.2_rk
    dx_min = 1.0e-4_rk

    coef = 0.0_rk
    coef(2:dim) = linspace(1.0_rk,dx_max/dx_min,dim-1)

    norm = dx_min/abs(x-y)*sum(coef)

    coef = coef/norm

    r(1) = x
    do i = 2,dim
       r(i) = r(i-1) + coef(i)*dx_min
    end do


  end function gridIK


  function gridEsler(x,y,dim) result(r)
    implicit none
    integer :: dim, i
    real(kind = rk) :: x,y,a,b
    real(kind = rk), dimension(dim) :: r

    b = log(400.0_rk)/ (dim-2)
    !b = log(10.0_rk)/ (dim-2)
    a = (y-x) / expm1(b*real(dim-1,kind=rk))

    do i = 1,dim
       r(i) = x + a*expm1(b*real((i-1),kind=rk))
    end do


  end function gridEsler

  function gridExtend(x,dr,dim) result(r)
    implicit none
    integer :: dim, i
    real(kind = rk) :: x,y,a,b,dr
    real(kind = rk), dimension(dim) :: r

    b = log(200.0_rk)/ (dim-2)
    a = dr/expm1(b)
    
    do i = 1,dim
       r(i) = x + a*expm1(b*real((i-1),kind=rk))
    end do
    
  end function gridExtend


  function logspace(x,y,dim) result(r)
    implicit none
    integer :: dim, i
    real(kind = rk) :: x,y,a
    real(kind = rk), dimension(dim) :: r
    
    a = (y/x)**(1.0_rk/(dim-1))
    
    r(1) = x
    do i = 2,dim
       r(i) = a*r(i-1)
    end do    
    
  end function logspace

  
  function interpolate2dc(r,rp,ind) result(z)
    implicit none
    real(kind = rk) :: r,rp,z
    integer :: ind
    
    !if (r>rp) then
    !   z = eval2dcspline(r,rp,1)
    !else
    !   z = eval2dcspline(rp,r,1)
    !end if
    
    z = (eval2dcspline(r,rp,ind) + eval2dcspline(rp,r,ind))/2
    
  end function interpolate2dc

  
  subroutine GaussHermite30(x,wt)
    implicit none
    real(kind=rk), dimension(30), intent(out) :: x,wt
    
    
    x = (/-6.8633452935298915810611083575550266_rk,&
         -6.1382792201239346203949923785375795_rk,&
         -5.5331471515674957251183335555803967_rk,&
         -4.9889189685899439444864971063309543_rk,&
         -4.4830553570925183418870376197091052_rk,&
         -4.0039086038612288152278760133218181_rk,&
         -3.5444438731553498869254009021683636_rk,&
         -3.0999705295864417486887333223746390_rk,&
         -2.6671321245356172005711064642208749_rk,&
         -2.2433914677615040724729799948250614_rk,&
         -1.8267411436036880388358804835061281_rk,&
         -1.4155278001981885119407251055475798_rk,&
         -1.0083382710467234618049896086964179_rk,&
         -0.6039210586255523077781556787573418_rk,&
         -0.2011285765488714855457630132436922_rk,&
         0.2011285765488714855457630132436922_rk,&
         0.6039210586255523077781556787573418_rk,&
         1.0083382710467234618049896086964179_rk,&
         1.4155278001981885119407251055475798_rk,&
         1.8267411436036880388358804835061281_rk,&
         2.2433914677615040724729799948250614_rk,&
         2.6671321245356172005711064642208749_rk,&
         3.0999705295864417486887333223746390_rk,&
         3.5444438731553498869254009021683636_rk,&
         4.0039086038612288152278760133218181_rk,&
         4.4830553570925183418870376197091052_rk,&
         4.9889189685899439444864971063309543_rk,&
         5.5331471515674957251183335555803967_rk,&
         6.1382792201239346203949923785375795_rk,&
         6.8633452935298915810611083575550266_rk/)


    wt = (/0.8342474710127617953407203967_rk,& 
         0.649097981554266700709961136_rk,& 
         0.569402691949640503966094892_rk,& 
         0.522525689331354549642402498_rk,& 
         0.49105799583288269650554983_rk,& 
         0.468374812564728816774690513_rk,& 
         0.451321035991188621287464588_rk,& 
         0.438177022652683703695367032_rk,& 
         0.4279180629327437485827730260_rk,&
         0.4198950037368240886418132650_rk,& 
         0.41367936361113893718433910583_rk,& 
         0.408981575003531602497229317388_rk,& 
         0.405605123325684436312140238733_rk,& 
         0.403419816924804022552760123219_rk,& 
         0.402346066701902927115350128334_rk,& 
         0.402346066701902927115350128334_rk,& 
         0.403419816924804022552760123219_rk,& 
         0.405605123325684436312140238733_rk,& 
         0.408981575003531602497229317388_rk,& 
         0.41367936361113893718433910583_rk,& 
         0.4198950037368240886418132650_rk,& 
         0.4279180629327437485827730260_rk,&
         0.438177022652683703695367032_rk,& 
         0.451321035991188621287464588_rk,& 
         0.468374812564728816774690513_rk,& 
         0.49105799583288269650554983_rk,& 
         0.522525689331354549642402498_rk,& 
         0.569402691949640503966094892_rk,& 
         0.649097981554266700709961136_rk,& 
         0.8342474710127617953407203967_rk/)
    
!     wt = (/0.83424747101276179534072039670392718_rk,&
!          0.64909798155426670070996113574601663_rk,&
!          0.56940269194964050396609489150119153_rk,&
!          0.52252568933135454964240249788468650_rk,&
!          0.49105799583288269650554982988781939_rk,&
!          0.46837481256472881677469051268093625_rk,&
!          0.45132103599118862128746458760617011_rk,&
!          0.43817702265268370369536703154242711_rk,&
!          0.42791806293274374858277302602122547_rk,&
!          0.41989500373682408864181326503301818_rk,&
!          0.41367936361113893718433910583495662_rk,&
!          0.40898157500353160249722931738840658_rk,&
!          0.40560512332568443631214023873335870_rk,&
!          0.40341981692480402255276012321945584_rk,&
!          0.40234606670190292711535012833385637_rk,&
!          0.40234606670190292711535012833385637_rk,&
!          0.40341981692480402255276012321945584_rk,&
!          0.40560512332568443631214023873335870_rk,&
!          0.40898157500353160249722931738840658_rk,&
!          0.41367936361113893718433910583495662_rk,&
!          0.41989500373682408864181326503301818_rk,&
!          0.42791806293274374858277302602122547_rk,&
!          0.43817702265268370369536703154242711_rk,&
!          0.45132103599118862128746458760617011_rk,&
!          0.46837481256472881677469051268093625_rk,&
!          0.49105799583288269650554982988781939_rk,&
!          0.52252568933135454964240249788468650_rk,&
!          0.56940269194964050396609489150119153_rk,&
!          0.64909798155426670070996113574601663_rk,&
!          0.83424747101276179534072039670392718_rk/)


  end subroutine GaussHermite30


  
  function GaussKronrodTot(f,a,b,r,rp,xbar,tau,x_grid,x_dim) result(I15)
    implicit none
    real(kind = rk), external :: f
    real(kind = rk) :: a,b,I15
    integer :: x_dim
    real(kind=rk) :: r, rp, xbar,tau,err_limit
    real(kind=rk), dimension(x_dim) :: x_grid

    GaussKronrodMinStep = (b-a)/2**20
    GaussKronrodValue = 0.0_rk
    !err_limit = 1.0e-9_rk
    !err_limit = minval((/ 1.0e-11_rk,tau/10000 /))    
    if (tau<1.0e-3) then
       err_limit = 1.0e-11_rk
    else
       err_limit =  1.0e-10_rk
    end if

    call GKintegrateTot15_7(f,a,b,err_limit,r,rp,xbar,tau,x_grid,x_dim)
    !call GKintegrateTot31_15t(f,a,b,err_limit,r,rp,xbar,tau,x_grid,x_dim)
    
    I15 = GaussKronrodValue
    !write(*,*) I15
    !pause

  end function GaussKronrodTot




  recursive subroutine GKintegrateTot15_7(f,a,b,err_limit,r,rp,xbar,tau,x_grid,x_dim)
    implicit none
    real(kind = rk), external :: f
    real(kind = rk), intent(in) :: a,b,err_limit
    real(kind = rk) :: I15,I7,an,bn,dx,ff,ker
    integer :: j,k
    real(kind = rk), dimension(15) :: r15,w15,r15n
    real(kind = rk), dimension(7) :: w7
    integer, intent(in) :: x_dim
    real(kind=rk), intent(in) :: r, rp, xbar,tau
    real(kind=rk), dimension(x_dim), intent(in) :: x_grid
    real(kind=rk) :: errgk,mid_point
    
    dx = abs(b-a)
    an = a
    bn = b
    I7 = 0.0_rk
    I15 = 0.0_rk
    
    call GK(w7,w15,r15)
    
    ker = (bn-an)/2
    
    r15n = ker*r15+(an+bn)/2
    
    k = 1    
    do j=1,15
       ff = ker*f(r15n(j),r,rp,xbar,tau,x_grid,x_dim)
       I15 = I15 + ff*w15(j)
       if (modulo(j,2) == 0) then
          I7 = I7 + ff*w7(k)
          k = k+1
       end if
    end do
    errgk = abs(I7-I15)
    !errgk = (200.0_rk*abs(I7-I15))**1.5_rk
    
    if (errgk>maxval((/err_limit,abs(I15)*1.0e-11_rk/)) .and. dx>GaussKronrodMinStep) then
    !if (errgk>err_limit .and. dx>GaussKronrodMinStep) then
       mid_point = (bn+an)/2
       call GKintegrateTot15_7(f,an,mid_point,err_limit/2,r,rp,xbar,tau,x_grid,x_dim)
       
       call GKintegrateTot15_7(f,mid_point,bn,err_limit/2,r,rp,xbar,tau,x_grid,x_dim)
    else
       GaussKronrodValue = GaussKronrodValue + I15
    end if
       
    
    
  end subroutine GKintegrateTot15_7


  recursive subroutine GKintegrateTot31_15(f,a,b,err_limit,r,rp,xbar,tau,x_grid,x_dim)
    implicit none
    real(kind = rk), external :: f
    real(kind = rk), intent(in) :: a,b,err_limit
    real(kind = rk) :: I15,I7,an,bn,dx,ff,ker
    integer :: j,k
    real(kind = rk), dimension(31) :: r15,w15,r15n
    real(kind = rk), dimension(15) :: w7
    integer, intent(in) :: x_dim
    real(kind=rk), intent(in) :: r, rp, xbar,tau
    real(kind=rk), dimension(x_dim), intent(in) :: x_grid
    real(kind=rk) :: errgk,mid_point
    

    !dx = (abs(b)+abs(a))
    dx = abs(b-a)
    an = a
    bn = b
    I7 = 0.0_rk
    I15 = 0.0_rk
    
    call GK31_15(w7,w15,r15)
    
    ker = (bn-an)/2
    
    r15n = ker*r15+(an+bn)/2
    
    k = 1    
    do j=1,31
       ff = ker*f(r15n(j),r,rp,xbar,tau,x_grid,x_dim)
       I15 = I15 + ff*w15(j)
       if (modulo(j,2) == 0) then
          I7 = I7 + ff*w7(k)
          k = k+1
       end if
    end do
    errgk = abs(I7-I15)
    !errgk = (200.0_rk*abs(I7-I15))**1.5_rk
    
    if (errgk>err_limit .and. dx>GaussKronrodMinStep) then
       mid_point = (bn+an)/2
       call GKintegrateTot31_15(f,an,mid_point,err_limit/2,r,rp,xbar,tau,x_grid,x_dim)
       
       call GKintegrateTot31_15(f,mid_point,bn,err_limit/2,r,rp,xbar,tau,x_grid,x_dim)
    else
       GaussKronrodValue = GaussKronrodValue + I15
    end if
       
    
    
  end subroutine GKintegrateTot31_15



  recursive subroutine GKintegrateTot31_15t(f,a,b,err_limit,r,rp,xbar,tau,x_grid,x_dim)
    implicit none
    real(kind = rk), external :: f
    real(kind = rk), intent(in) :: a,b,err_limit
    real(kind = rk) :: I15,I7,an,bn,dx,ff,ker,ff2
    integer :: j,k
    real(kind = rk), dimension(31) :: r15,w15,r15n
    real(kind = rk), dimension(15) :: w7
    integer, intent(in) :: x_dim
    real(kind=rk), intent(in) :: r, rp, xbar,tau
    real(kind=rk), dimension(x_dim), intent(in) :: x_grid
    real(kind=rk) :: errgk,mid_point
    

    !dx = (abs(b)+abs(a))
    dx = abs(b-a)
    an = a
    bn = b
    I7 = 0.0_rk
    I15 = 0.0_rk
    
    call GK31_15(w7,w15,r15)
    
    ker = (bn-an)/2
    
    r15n = ker*r15+(an+bn)/2
    
    k = 1    
    do j=1,15
       ff = ker*f(r15n(2*j-1),r,rp,xbar,tau,x_grid,x_dim)
       ff2 = ker*f(r15n(2*j),r,rp,xbar,tau,x_grid,x_dim)
       I15 = I15 + ff*w15(2*j-1)+ff2*w15(2*j)
       I7 = I7 + ff2*w7(j)
    end do
    I15 = I15 + ker*f(r15n(31),r,rp,xbar,tau,x_grid,x_dim)*w15(31)
    errgk = abs(I7-I15)
    !errgk = (200.0_rk*abs(I7-I15))**1.5_rk    

    if (errgk>maxval((/err_limit,abs(I15)*1.0e-8_rk/)) .and. dx>GaussKronrodMinStep) then
       mid_point = (bn+an)/2
       call GKintegrateTot31_15t(f,an,mid_point,err_limit/2,r,rp,xbar,tau,x_grid,x_dim)
       
       call GKintegrateTot31_15t(f,mid_point,bn,err_limit/2,r,rp,xbar,tau,x_grid,x_dim)
    else
       GaussKronrodValue = GaussKronrodValue + I15
    end if
       
    
    
  end subroutine GKintegrateTot31_15t



  
  function GaussKronrodFree(f,a,b,r,rp,xbar,tau) result(I15)
    implicit none
    real(kind = rk), external :: f
    real(kind = rk) :: a,b,I15
    real(kind=rk) :: r, rp, xbar,tau,err_limit
    
    GaussKronrodValue = 0.0_rk
    err_limit = minval((/ 1.0e-11_rk,tau/10000 /))
    !write(*,*) 'in'
    call GKintegrateFree(f,a,b,err_limit,r,rp,xbar,tau)
    !write(*,*) 'out'
    
    I15 = GaussKronrodValue

  end function GaussKronrodFree




  recursive subroutine GKintegrateFree(f,a,b,err_limit,r,rp,xbar,tau)
    implicit none
    real(kind = rk), external :: f
    real(kind = rk), intent(in) :: a,b,err_limit
    real(kind = rk) :: I15,I7,an,bn,dx,ff,ker
    integer :: j,k
    real(kind = rk), dimension(15) :: r15,w15,r15n
    real(kind = rk), dimension(7) :: w7
    real(kind=rk), intent(in) :: r, rp, xbar,tau
    real(kind=rk) :: errgk,mid_point
    

    !dx = (abs(b)+abs(a))
    dx = abs(b-a)
    an = a
    bn = b
    I7 = 0.0_rk
    I15 = 0.0_rk
    
    call GK(w7,w15,r15)
    
    ker = (bn-an)/2
    
    r15n = ker*r15+(an+bn)/2
    
    k = 1    
    do j=1,15
       ff = ker*f(r15n(j),r,rp,xbar,tau)
       I15 = I15 + ff*w15(j)
       if (modulo(j,2) == 0) then
          I7 = I7 + ff*w7(k)
          k = k+1
       end if
    end do
    errgk = abs(I7-I15)
    !errgk = (200.0_rk*abs(I7-I15))**1.5_rk
    
    if (errgk>err_limit .and. dx>1.0e-6_rk) then
       mid_point = (bn+an)/2
       call GKintegrateFree(f,an,mid_point,err_limit/2,r,rp,xbar,tau)
       
       call GKintegrateFree(f,mid_point,bn,err_limit/2,r,rp,xbar,tau)
    else
       GaussKronrodValue = GaussKronrodValue + I15
    end if
       
    
    
  end subroutine GKintegrateFree


  subroutine GK(w7,w15,r15)
    implicit none
    real(kind = rk), dimension(7), intent(out) :: w7
    real(kind = rk), dimension(15), intent(out) :: w15,r15

    !r7 = (/-0.949107912342758524526189684047851_rk,&
    !     -0.741531185599394439863864773280788_rk,&
    !     -0.405845151377397166906606412076961_rk,&
    !     0.000000000000000000000000000000000_rk,&
    !     0.405845151377397166906606412076961_rk,&
    !     0.741531185599394439863864773280788_rk,&
    !     0.949107912342758524526189684047851_rk/)


    w7 = (/0.129484966168869693270611432679082_rk,&
         0.279705391489276667901467771423780_rk,&
         0.381830050505118944950369775488975_rk,&
         0.417959183673469387755102040816327_rk,&
         0.381830050505118944950369775488975_rk,&
         0.279705391489276667901467771423780_rk,&
         0.129484966168869693270611432679082_rk/)

    r15 = (/-0.991455371120812639206854697526329_rk,&
         -0.949107912342758524526189684047851_rk,&
         -0.864864423359769072789712788640926_rk,&
         -0.741531185599394439863864773280788_rk,&
         -0.586087235467691130294144838258730_rk,&
         -0.405845151377397166906606412076961_rk,&
         -0.207784955007898467600689403773245_rk,&
         0.000000000000000000000000000000000_rk,&
         0.207784955007898467600689403773245_rk,&
         0.405845151377397166906606412076961_rk,&
         0.586087235467691130294144838258730_rk,&
         0.741531185599394439863864773280788_rk,&
         0.864864423359769072789712788640926_rk,&
         0.949107912342758524526189684047851_rk,&
         0.991455371120812639206854697526329_rk/)
    

    w15 = (/0.022935322010529224963732008058970_rk,&
         0.063092092629978553290700663189204_rk,&
         0.104790010322250183839876322541518_rk,&
         0.140653259715525918745189590510238_rk,&
         0.169004726639267902826583426598550_rk,&
         0.190350578064785409913256402421014_rk,&
         0.204432940075298892414161999234649_rk,&
         0.209482141084727828012999174891714_rk,&
         0.204432940075298892414161999234649_rk,&
         0.190350578064785409913256402421014_rk,&
         0.169004726639267902826583426598550_rk,&
         0.140653259715525918745189590510238_rk,&
         0.104790010322250183839876322541518_rk,&
         0.063092092629978553290700663189204_rk,&
         0.022935322010529224963732008058970_rk/)

  end subroutine GK



  subroutine GK31_15(w15,w31,r31)
    implicit none
    real(kind = rk), dimension(15), intent(out) :: w15
    real(kind = rk), dimension(31), intent(out) :: w31,r31

    r31 = (/  0.998002298693397060285172840152271_rk,&
         0.987992518020485428489565718586613_rk,&
         0.967739075679139134257347978784337_rk,&
         0.937273392400705904307758947710209_rk,&
         0.897264532344081900882509656454496_rk,&
         0.848206583410427216200648320774217_rk,&
         0.790418501442465932967649294817947_rk,&
         0.724417731360170047416186054613938_rk,&
         0.650996741297416970533735895313275_rk,&
         0.570972172608538847537226737253911_rk,&
         0.485081863640239680693655740232351_rk,&
         0.394151347077563369897207370981045_rk,&
         0.299180007153168812166780024266389_rk,&
         0.201194093997434522300628303394596_rk,&
         0.101142066918717499027074231447392_rk,&
         0.000000000000000000000000000000000_rk,&
         0.101142066918717499027074231447392_rk,&
         0.201194093997434522300628303394596_rk,&
         0.299180007153168812166780024266389_rk,&
         0.394151347077563369897207370981045_rk,&
         0.485081863640239680693655740232351_rk,&
         0.570972172608538847537226737253911_rk,&
         0.650996741297416970533735895313275_rk,&
         0.724417731360170047416186054613938_rk,&
         0.790418501442465932967649294817947_rk,&
         0.848206583410427216200648320774217_rk,&
         0.897264532344081900882509656454496_rk,&
         0.937273392400705904307758947710209_rk,&
         0.967739075679139134257347978784337_rk,&
         0.987992518020485428489565718586613_rk,&
         0.998002298693397060285172840152271_rk /)


    w31 = (/  0.005377479872923348987792051430128_rk,&
         0.015007947329316122538374763075807_rk,&
         0.025460847326715320186874001019653_rk,&
         0.035346360791375846222037948478360_rk,&
         0.044589751324764876608227299373280_rk,&
         0.053481524690928087265343147239430_rk,&
         0.062009567800670640285139230960803_rk,&
         0.069854121318728258709520077099147_rk,&
         0.076849680757720378894432777482659_rk,&
         0.083080502823133021038289247286104_rk,&
         0.088564443056211770647275443693774_rk,&
         0.093126598170825321225486872747346_rk,&
         0.096642726983623678505179907627589_rk,&
         0.099173598721791959332393173484603_rk,&
         0.100769845523875595044946662617570_rk,&
         0.101330007014791549017374792767493_rk,&
         0.100769845523875595044946662617570_rk,&
         0.099173598721791959332393173484603_rk,&
         0.096642726983623678505179907627589_rk,&
         0.093126598170825321225486872747346_rk,&
         0.088564443056211770647275443693774_rk,&
         0.083080502823133021038289247286104_rk,&
         0.076849680757720378894432777482659_rk,&
         0.069854121318728258709520077099147_rk,&
         0.062009567800670640285139230960803_rk,&
         0.053481524690928087265343147239430_rk,&
         0.044589751324764876608227299373280_rk,&
         0.035346360791375846222037948478360_rk,&
         0.025460847326715320186874001019653_rk,&
         0.015007947329316122538374763075807_rk,&
         0.005377479872923348987792051430128_rk/)


    w15 = (/  0.030753241996117268354628393577204_rk,&
         0.070366047488108124709267416450667_rk,&
         0.107159220467171935011869546685869_rk,&
         0.139570677926154314447804794511028_rk,&
         0.166269205816993933553200860481209_rk,&
         0.186161000015562211026800561866423_rk,&
         0.198431485327111576456118326443839_rk,&
         0.202578241925561272880620199967519_rk,&
         0.198431485327111576456118326443839_rk,&
         0.186161000015562211026800561866423_rk,&
         0.166269205816993933553200860481209_rk,&
         0.139570677926154314447804794511028_rk,&
         0.107159220467171935011869546685869_rk,&
         0.070366047488108124709267416450667_rk,&
         0.030753241996117268354628393577204_rk /)

  end subroutine GK31_15

  
  subroutine LUdecomp(A,dim)
    implicit none
    integer, intent(in) :: dim
    real(kind = rk), dimension(dim,dim),intent(inout) :: A
    real(kind = rk), dimension(dim,dim) :: L,U
    integer :: i,j,k

    L = 0.0_rk
    U = 0.0_rk

    do i=1,dim
  
       ! non-zero diagonal assumed
  
       L(i,i) = 1.0_rk
  
       do j=i,dim
          U(i,j) = A(i,j)
          do k=1,i-1
             U(i,j) = U(i,j) - L(i,k)*U(k,j)
          end do
       end do
  
       do j=i+1,dim
          L(j,i) = A(j,i)
          do k=1,i-1
             L(j,i) = L(j,i) - L(j,k)*U(k,i)
          end do
          L(j,i) = L(j,i)/U(i,i)
       end do
    end do
    
    do i=1,dim
       do j=i,dim
          A(j,i) = L(j,i)
          A(i,j) = U(i,j)
       end do
    end do

  end subroutine LUdecomp


  
  FUNCTION det(A,M) RESULT( detvalue )
    IMPLICIT NONE
    !EXTERNAL :: SGETRF
    REAL(KIND=RK), DIMENSION(M,M), INTENT(out) :: A
    REAL(KIND=RK)                             :: detvalue
    INTEGER, INTENT(in)              :: M
    !INTEGER                          :: INFO,i,j
    integer :: i,j
    !INTEGER, DIMENSION(M)            :: IPIV

    !DO i = 1,M
    !  IPIV(i) = i
    !END DO
    !INFO = 0

    ! LU decomposition
    !CALL SGETRF(M,M,A,M,IPIV,INFO)
    
    call Pivoting(A,M,j)
    call LUdecomp(A,M)
    
    ! value of the determinant
    !j = 0
    detvalue = 1.0_rk
    DO i = 1,M
       !IF (IPIV(i) .NE. i) THEN
       !  j = j + 1
       !END IF
       
       detvalue = detvalue*A(i,i)
       
    END DO
    
    detvalue = ((-1.0)**j)*detvalue
    
  END FUNCTION det

  subroutine Pivoting(A,M,j)
    implicit none
    integer, intent(in) :: M
    real(kind=rk), dimension(M,M), intent(inout) :: A
    real(kind=rk), dimension(1,M) :: row
    real(kind=rk), dimension(M,1) :: column
    integer, intent(out) :: j
    integer :: i,k,mm,nn,ind
    real(kind=rk) :: largest
    j=0
    mm=1
    nn=1
    
    do ind=1,M-1
       largest=0.0_rk 
       do i=ind,M
          do k=ind,M
             if (abs(A(i,k))>largest) then
                largest=abs(A(i,k))
                mm=i
                nn=k
             end if
          end do
       end do
    
       if (mm .ne. ind) then
          row(1,:)=A(mm,:)
          A(mm,:)=A(ind,:)
          A(ind,:)=row(1,:)
          j=j+1
       end if
       
       if (nn .ne. ind) then
          column(:,1)=A(:,nn)
          A(:,nn)=A(:,ind)
          A(:,ind)=column(:,1)
          j=j+1
       end if
    end do
    
  end subroutine Pivoting


  function index_check(ind,limit_down,limit_up) result(inbounds)
    implicit none
    integer :: ind, limit_down, limit_up
    logical :: inbounds
    
    if (ind >=1 .and. ind >= limit_down .and. ind <= limit_up) then
       inbounds = .true.
    else
       inbounds = .false.
    end if
    
  end function index_check

  function GetNucleiBeadLabel(i_e,M_e,M_n) result(i_n)
    implicit none
    integer(kind = isoluku),intent(in) :: i_e,M_e,M_n
    integer(kind = isoluku) :: i_n,jakso
    
    ! M_e should be >= M_n
    ! and jakso should be exactly an integer not (0.9 ~ 1)
    jakso = M_e/M_n
    
    if (jakso == 1) then
       i_n = i_e
    elseif (jakso == M_e) then
       i_n = 1
    else
       i_n = ceiling(real(i_e,kind=rk)/jakso)
    end if
    
  end function GetNucleiBeadLabel
  
  
end module math
