module periodicbc
  
  use math
  implicit none
  real(kind =rk), dimension(3), private :: PBCbox    
  logical, save :: periodic_on, PBCcontinue
  
contains

  subroutine SetBox(Lx,Ly,Lz)
    implicit none
    real(kind = rk), intent(in) :: Lx,Ly,Lz
    
    PBCbox = (/Lx,Ly,Lz/)
    
  end subroutine SetBox

  subroutine GetTheBox(L)
    implicit none
    real(kind = rk), dimension(3):: L
    
    L = PBCbox
    
  end subroutine GetTheBox

  subroutine PutInBox(r)
    implicit none
    real(kind = rk), dimension(3), intent(inout) :: r
    
    r = modulo(r,PBCbox)
    
  end subroutine PutInBox
  
  function periodic_distance(r1,r2) result(dr)
    implicit none
    real(kind=rk), dimension(3) :: r1,r2
    real(kind=rk) :: dx,dr
    integer :: i

    dr = 0.0_rk
    do i=1,3
       dx = r2(i)-r1(i)
       if ( abs(dx) > PBCbox(i)/2 ) then
          if (dx>0.0_rk) then 
             !r2(i) = r2(i) - PBCbox(i)
             dx = dx -  PBCbox(i)
          else
             !r2(i) = r2(i) + PBCbox(i)
             dx = dx +  PBCbox(i)
          end if
       end if
       dr = dr + dx**2
    end do
    dr = sqrt( dr )

  end function periodic_distance

  function periodic_distance2(r1,r2) result(dr)
    implicit none
    real(kind=rk), dimension(3) :: r1,r2
    real(kind=rk) :: dx,dr
    integer :: i
    
    dr = 0.0_rk
    do i=1,3
       dx = r2(i)-r1(i)
       if ( abs(dx) > PBCbox(i)/2 ) then
          if (dx>0.0_rk) then 
             !r2(i) = r2(i) - PBCbox(i)
             dx = dx -  PBCbox(i)
          else
             !r2(i) = r2(i) + PBCbox(i)
             dx = dx +  PBCbox(i)
          end if
       end if
       dr = dr + dx**2
    end do

  end function periodic_distance2

  function periodic_mean(r1,r2) result(r)
    implicit none
    real(kind=rk), dimension(3) :: r1,r2,r,rt
    real(kind=rk) :: dx
    integer :: i

    rt = r2

    do i=1,3
       dx = r2(i)-r1(i)
       if ( abs(dx) > PBCbox(i)/2 ) then
          if (dx>0.0_rk) then 
             rt(i) = r2(i) - PBCbox(i)             
          else
             rt(i) = r2(i) + PBCbox(i)
          end if
       end if
    end do
    
    r = ( r1 + rt )/2

  end function periodic_mean

  function periodic_difference(r1,r2) result(r)
    implicit none
    real(kind=rk), dimension(3) :: r1,r2,r
    real(kind=rk) :: dx
    integer :: i
    
    r = 0.0_rk

    do i=1,3
       dx = r2(i)-r1(i)
       if ( abs(dx) > PBCbox(i)/2 ) then
          if (dx>0.0_rk) then 
             !r2(i) = r2(i) - PBCbox(i)
             dx = dx -  PBCbox(i)
          else
             !r2(i) = r2(i) + PBCbox(i)
             dx = dx +  PBCbox(i)
          end if
       end if
       r(i) = dx
    end do
    

  end function periodic_difference


  
  function periodicImagedistance(r1,r2,n_x,n_y,n_z) result(dr)
    implicit none
    real(kind=rk), dimension(3) :: r1,r2
    real(kind=rk) :: dx,dr
    integer :: i
    integer :: n_x,n_y,n_z
    real(kind = rk), dimension(3) :: n

    n = (/ real(n_x,kind=rk),real(n_y,kind=rk),real(n_z,kind=rk)/)
    dr = 0.0_rk
    do i=1,3
       dx = r2(i)-r1(i)
       if ( abs(dx) > PBCbox(i)/2 ) then
          if (dx>0.0_rk) then 
             !r2(i) = r2(i) - PBCbox(i)
             dx = dx -  PBCbox(i)
          else
             !r2(i) = r2(i) + PBCbox(i)
             dx = dx +  PBCbox(i)
          end if
       end if
       dx = abs(dx) + n(i)*PBCbox(i)
       dr = dr + dx**2
    end do
    dr = sqrt( dr )

  end function periodicImagedistance

  function periodicImagedistance2(r1,r2,n_x,n_y,n_z) result(dr)
    implicit none
    real(kind=rk), dimension(3) :: r1,r2
    real(kind=rk) :: dx,dr
    integer :: i
    integer :: n_x,n_y,n_z
    real(kind = rk), dimension(3) :: n

    n = (/ real(n_x,kind=rk),real(n_y,kind=rk),real(n_z,kind=rk)/)    
    dr = 0.0_rk
    do i=1,3
       dx = r2(i)-r1(i)
       if ( abs(dx) > PBCbox(i)/2 ) then
          if (dx>0.0_rk) then 
             !r2(i) = r2(i) - PBCbox(i)
             dx = dx -  PBCbox(i)
          else
             !r2(i) = r2(i) + PBCbox(i)
             dx = dx +  PBCbox(i)
          end if
       end if
       dr = dr + (abs(dx)+n(i)*PBCbox(i))**2
    end do

  end function periodicImagedistance2


end module periodicbc
