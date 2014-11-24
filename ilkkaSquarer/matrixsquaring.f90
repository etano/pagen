module matrixsquaring

  implicit none

contains
  
  subroutine MatrixSquaringIntegration(r,u0,tau,r_dim)
    use math
    implicit none
    integer, intent(in) :: r_dim
    real(kind = rk), intent(in) :: tau
    real(kind = rk), dimension(r_dim), intent(in) :: r
    real(kind = rk), dimension(r_dim,r_dim), intent(out) :: u0
    real(kind = rk) :: sqrt2ltau,xbar
    real(kind = rk), dimension(30) :: x30,wt30,gs
    !real(kind = rk), dimension(64) :: x64,wt64,gs
    real(kind = rk) :: Ival
    real(kind = rk) :: raja,ylaraja,alaraja
    integer :: i,j
    
    sqrt2ltau = sqrt(2.0_rk*lambdaMS*tau)        
    
    call GaussHermite30(x30,wt30)
    !call GaussHermite64(x64,wt64)
    ! x30 is in the range [-6.863*sqrt2ltau,6.863*sqrt2ltau]
    x30 = sqrt2ltau*x30 
    
    !raja = the limit after which the Gauss-Hermite integration is used
    !raja = 9.0_rk*sqrt2ltau
    raja = 12.0_rk*sqrt2ltau

    
    ! set the upper bound (ylaraja) for the Gauss-Kronrod integration
    !  -- (upper bound = xbar + ylaraja)
    !ylaraja = 13.0_rk*sqrt2ltau
    ylaraja = 20.0_rk*sqrt2ltau
    !ylaraja = x30(30)
    !ylaraja = x64(64)


    do i = 1,r_dim
       do j = 1,r_dim
          
          xbar = (r(i)+r(j))/2                    
          
          if (tau<1.0e-03) then
             
             if (xbar > raja) then
                ! Gauss-Hermite integration
                !gs = FreeIntegrand(x30,r(i),r(j),xbar,tau,30)
                !gs = FreeIntegrand(x64,r(i),r(j),xbar,tau,64)
                
                !Check that the norm is one
                !write(*,*) 'Hermite ', sqrt2ltau*dot_product(gs,wt30)
                !pause

                !gs = gs*expm1UIntegrand(x30,r(i),r(j),xbar,tau,30,r,r_dim)
                !gs = gs*expm1UIntegrand(x64,r(i),r(j),xbar,tau,64,r,r_dim)


                gs = SmallTau(x30,r(i),r(j),xbar,tau,30,r,r_dim)
                Ival = sqrt2ltau*dot_product(gs,wt30)

                !Ival = sqrt2ltau*sum(gs*wt30)
             else if (xbar<1.0e-40_rk) then ! i.e. xbar==0.0
                gs = SmallTau(abs(x30),r(i),r(j),xbar,tau,30,r,r_dim)
                Ival = sqrt2ltau*dot_product(gs,wt30)/2
             else
                
                !alaraja = -minval((/xbar,x30(30)/))                
                !alaraja = -minval((/xbar,ylaraja/))                
                alaraja = -xbar
                
                ! Gauss-Kronrod integration
                Ival = GaussKronrodTot(TotalIntegrandSmallTau,alaraja,ylaraja,r(i),r(j),xbar,tau,r,r_dim)
                
                !write(*,*) 'Kronrod ',GaussKronrodFree(FreeIntegrandgk,alaraja,ylaraja,r(i),r(j),xbar,tau)
                !pause

                
             end if
             
             u0(i,j) = -log1p(Ival)            
          else
             
             if (xbar > raja) then
                ! Gauss-Hermite integration
                !gs = FreeIntegrand(x30,r(i),r(j),xbar,tau,30)
                
                !Check that the norm is one
                !write(*,*) 'Hermite ', sqrt2ltau*dot_product(gs,wt30)
                !pause
                
                !gs = gs*expUIntegrand(x30,r(i),r(j),xbar,tau,30,r,r_dim)

                gs = LargeTau(x30,r(i),r(j),xbar,tau,30,r,r_dim)
                Ival = sqrt2ltau*dot_product(gs,wt30)

                !Ival = sqrt2ltau*sum(gs*wt30)
             else if (xbar<1.0e-40_rk) then ! i.e. xbar==0.0
                gs = LargeTau(abs(x30),r(i),r(j),xbar,tau,30,r,r_dim)
                Ival = sqrt2ltau*dot_product(gs,wt30)/2
             else

                !alaraja = -minval((/xbar,ylaraja/))
                alaraja = -xbar
                !alaraja = -minval((/xbar,x30(30)/))                

                ! Gauss-Kronrod integration                
                Ival = GaussKronrodTot(TotalIntegrandLargerTau,alaraja,ylaraja,r(i),r(j),xbar,tau,r,r_dim)
                
                !write(*,*) 'Kronrod ',GaussKronrodFree(FreeIntegrandgk,alaraja,ylaraja,r(i),r(j),xbar,tau)
                !pause
                
             end if
             
             u0(i,j) = -log(Ival)
          end if
          

          u0(j,i) = u0(i,j)
          
       end do
    end do
    
    ! symmetrize
    do j = 1,r_dim
       do i = j,r_dim
          Ival = (u0(i,j)+u0(j,i))/2
          u0(i,j) = Ival
          u0(j,i) = Ival          
       end do
    end do
    
    
  end subroutine MatrixSquaringIntegration


  subroutine MatrixSquaringIntegrationWithBetaDeriv(r,u0,du0,tau,r_dim)
    use math
    implicit none
    integer, intent(in) :: r_dim
    real(kind = rk), intent(in) :: tau
    real(kind = rk), dimension(r_dim), intent(in) :: r
    real(kind = rk), dimension(r_dim,r_dim), intent(out) :: u0,du0
    real(kind = rk) :: sqrt2ltau,xbar
    real(kind = rk), dimension(30) :: x30,wt30,gs
    !real(kind = rk), dimension(64) :: x64,wt64,gs
    real(kind = rk) :: Ival,Ival2
    real(kind = rk) :: raja,ylaraja,alaraja
    integer :: i,j
    
    sqrt2ltau = sqrt(2.0_rk*lambdaMS*tau)        
    
    call GaussHermite30(x30,wt30)
    !call GaussHermite64(x64,wt64)
    ! x30 is in the range [-6.863*sqrt2ltau,6.863*sqrt2ltau]
    x30 = sqrt2ltau*x30 
    
    !raja = the limit after which the Gauss-Hermite integration is used
    raja = 9.0_rk*sqrt2ltau

    
    ! set the upper bound (ylaraja) for the Gauss-Kronrod integration
    !  -- (upper bound = xbar + ylaraja)
    ylaraja = 13.0_rk*sqrt2ltau
    !ylaraja = x30(30)
    !ylaraja = x64(64)


    do i = 1,r_dim
       do j = 1,r_dim
          
          xbar = (r(i)+r(j))/2                    
          
          if (tau<1.0e-03) then
             
             if (xbar > raja) then
                ! Gauss-Hermite integration
                !gs = FreeIntegrand(x30,r(i),r(j),xbar,tau,30)
                !gs = FreeIntegrand(x64,r(i),r(j),xbar,tau,64)
                
                !Check that the norm is one
                !write(*,*) 'Hermite ', sqrt2ltau*dot_product(gs,wt30)
                !pause

                !gs = gs*expm1UIntegrand(x30,r(i),r(j),xbar,tau,30,r,r_dim)
                !gs = gs*expm1UIntegrand(x64,r(i),r(j),xbar,tau,64,r,r_dim)


                gs = SmallTau(x30,r(i),r(j),xbar,tau,30,r,r_dim)
                Ival = sqrt2ltau*dot_product(gs,wt30)
                gs = LargeTau2(x30,r(i),r(j),xbar,tau,30,r,r_dim)
                Ival2 = sqrt2ltau*dot_product(gs,wt30)

                !Ival = sqrt2ltau*sum(gs*wt30)
             else if (xbar<1.0e-40_rk) then ! i.e. xbar==0.0
                gs = SmallTau(abs(x30),r(i),r(j),xbar,tau,30,r,r_dim)
                Ival = sqrt2ltau*dot_product(gs,wt30)/2
                gs = LargeTau2(x30,r(i),r(j),xbar,tau,30,r,r_dim)
                Ival2 = sqrt2ltau*dot_product(gs,wt30)/2
             else
                
                !alaraja = -minval((/xbar,x30(30)/))                
                !alaraja = -minval((/xbar,ylaraja/))                
                alaraja = -xbar
                
                ! Gauss-Kronrod integration
                Ival = GaussKronrodTot(TotalIntegrandSmallTau,alaraja,ylaraja,r(i),r(j),xbar,tau,r,r_dim)

                Ival2 = GaussKronrodTot(TotalIntegrandLargerTau2,alaraja,ylaraja,r(i),r(j),xbar,tau,r,r_dim)
                
                !write(*,*) 'Kronrod ',GaussKronrodFree(FreeIntegrandgk,alaraja,ylaraja,r(i),r(j),xbar,tau)
                !pause

                
             end if
             
             u0(i,j) = -log1p(Ival)
             du0(i,j) = -exp(u0(i,j))*Ival2/2
          else
             
             if (xbar > raja) then
                ! Gauss-Hermite integration
                !gs = FreeIntegrand(x30,r(i),r(j),xbar,tau,30)
                
                !Check that the norm is one
                !write(*,*) 'Hermite ', sqrt2ltau*dot_product(gs,wt30)
                !pause
                
                !gs = gs*expUIntegrand(x30,r(i),r(j),xbar,tau,30,r,r_dim)

                gs = LargeTau(x30,r(i),r(j),xbar,tau,30,r,r_dim)
                Ival = sqrt2ltau*dot_product(gs,wt30)
                gs = LargeTau2(x30,r(i),r(j),xbar,tau,30,r,r_dim)
                Ival2 = sqrt2ltau*dot_product(gs,wt30)

                !Ival = sqrt2ltau*sum(gs*wt30)
             else if (xbar<1.0e-40_rk) then ! i.e. xbar==0.0
                gs = LargeTau(abs(x30),r(i),r(j),xbar,tau,30,r,r_dim)
                Ival = sqrt2ltau*dot_product(gs,wt30)/2
                gs = LargeTau2(x30,r(i),r(j),xbar,tau,30,r,r_dim)
                Ival2 = sqrt2ltau*dot_product(gs,wt30)/2
             else

                !alaraja = -minval((/xbar,ylaraja/))
                alaraja = -xbar
                !alaraja = -minval((/xbar,x30(30)/))                

                ! Gauss-Kronrod integration                
                Ival = GaussKronrodTot(TotalIntegrandLargerTau,alaraja,ylaraja,r(i),r(j),xbar,tau,r,r_dim)

                Ival2 = GaussKronrodTot(TotalIntegrandLargerTau2,alaraja,ylaraja,r(i),r(j),xbar,tau,r,r_dim)
                
                !write(*,*) 'Kronrod ',GaussKronrodFree(FreeIntegrandgk,alaraja,ylaraja,r(i),r(j),xbar,tau)
                !pause
                
             end if
             
             u0(i,j) = -log(Ival)
             du0(i,j) = -exp(u0(i,j))*Ival2/2
             
          end if
       
          !write(*,*) exp(u0(i,j)),Ival2
          !pause
          

          !u0(j,i) = u0(i,j)
          !du0(j,i) = du0(i,j)
          
       end do
    end do

    ! symmetrize
    do j = 1,r_dim
       do i = j,r_dim
          Ival = (u0(i,j)+u0(j,i))/2
          u0(i,j) = Ival
          u0(j,i) = Ival          
          Ival2 = (du0(i,j)+du0(j,i))/2
          du0(i,j) = Ival2
          du0(j,i) = Ival2         
       end do
    end do
    
    
    
  end subroutine MatrixSquaringIntegrationWithBetaDeriv

  function FreeIntegrand(s,r,rp,xbar,tau,s_dim) result(value)
    use math
    implicit none
    integer :: s_dim,i
    real(kind=rk) :: r, rp, xbar,tau
    real(kind=rk) :: xs, ltau
    real(kind=rk), dimension(s_dim) :: s,value
    
    
    ltau = lambdaMS*tau
    
    do i=1,s_dim
       xs = xbar+s(i)
       value(i)=sqrt(2.0_rk/pi/ltau**3)*&
            xs**2*exp(-s(i)**2/(2.0_rk*ltau))*&
            m0(r*xs/(2.0_rk*ltau))*m0(xs*rp/(2.0_rk*ltau))&
            /m0(r*rp/(4.0_rk*ltau))
    end do
    

  end function FreeIntegrand

  function FreeIntegrandgk(s,r,rp,xbar,tau) result(value)
    use math
    implicit none
    real(kind=rk) :: r, rp, xbar,tau
    real(kind=rk) :: xs, ltau
    real(kind=rk) :: s,value
    
    
    ltau = lambdaMS*tau
    
    xs = xbar+s
    value=sqrt(2.0_rk/pi/ltau**3)*&
         xs**2*exp(-s**2/(2.0_rk*ltau))*&
         m0(r*xs/(2.0_rk*ltau))*m0(xs*rp/(2.0_rk*ltau))&
         /m0(r*rp/(4.0_rk*ltau))
    

  end function FreeIntegrandgk
  
  function TotalIntegrandSmallTau(s,r,rp,xbar,tau,x_grid,x_dim) result(value)
    use math
    implicit none
    integer :: x_dim
    real(kind=rk) :: r, rp, xbar,tau
    real(kind=rk) :: xs, ltau
    real(kind=rk) :: s,value
    real(kind=rk), dimension(x_dim) :: x_grid
       
    ltau = lambdaMS*tau
    
    xs = xbar + s
    value = calc_u0(x_grid,r,xs,tau,x_dim) &
         + calc_u0(x_grid,xs,rp,tau,x_dim)
    value = expm1(-value)
    value = value*sqrt(2.0_rk/pi/ltau**3)*&
         xs**2*exp(-s**2/(2.0_rk*ltau))*&
         m0(r*xs/(2.0_rk*ltau))*m0(xs*rp/(2.0_rk*ltau))&
         /m0(r*rp/(4.0_rk*ltau))
    
  end function TotalIntegrandSmallTau

  function TotalIntegrandLargerTau(s,r,rp,xbar,tau,x_grid,x_dim) result(value)
    use math
    implicit none
    integer :: x_dim
    real(kind=rk) :: r, rp, xbar,tau
    real(kind=rk) :: xs, ltau
    real(kind=rk) :: s,value
    real(kind=rk), dimension(x_dim) :: x_grid

    ltau = lambdaMS*tau
    
    xs = xbar + s
    value = calc_u0(x_grid,r,xs,tau,x_dim) &
         + calc_u0(x_grid,xs,rp,tau,x_dim)
    value = exp(-value)
    value= value*sqrt(2.0_rk/pi/ltau**3)*&
         xs**2*exp(-s**2/(2.0_rk*ltau))*&
         m0(r*xs/(2.0_rk*ltau))*m0(xs*rp/(2.0_rk*ltau))&
         /m0(r*rp/(4.0_rk*ltau))

  end function TotalIntegrandLargerTau

  function expUIntegrand(s,r,rp,xbar,tau,s_dim,x_grid,x_dim) result(value)
    use math
    implicit none
    integer :: s_dim,x_dim,i
    real(kind=rk) :: r, rp, xbar,tau
    real(kind=rk) :: xs
    real(kind=rk), dimension(s_dim) :: s,value
    real(kind=rk), dimension(x_dim) :: x_grid
    
    do i=1,s_dim       
       xs = xbar + s(i)
       value(i) = calc_u0(x_grid,r,xs,tau,x_dim) &
            + calc_u0(x_grid,xs,rp,tau,x_dim)
    end do
    value = exp(-value)    

  end function expUIntegrand

  function expm1UIntegrand(s,r,rp,xbar,tau,s_dim,x_grid,x_dim) result(value)
    use math
    implicit none
    integer :: s_dim,x_dim,i
    real(kind=rk) :: r, rp, xbar,tau
    real(kind=rk) :: xs
    real(kind=rk), dimension(s_dim) :: s,value
    real(kind=rk), dimension(x_dim) :: x_grid
    
    
    do i=1,s_dim       
       xs = xbar + s(i)
       value(i) = calc_u0(x_grid,r,xs,tau,x_dim) &
            + calc_u0(x_grid,xs,rp,tau,x_dim)
       value(i) = expm1(-value(i))
    end do

  end function expm1UIntegrand

  

  function SmallTau(s,r,rp,xbar,tau,s_dim,x_grid,x_dim) result(value)
    use math
    implicit none
    integer :: s_dim,x_dim,i
    real(kind=rk) :: r, rp, xbar,tau
    real(kind=rk) :: xs,ltau
    real(kind=rk), dimension(s_dim) :: s,value
    real(kind=rk), dimension(x_dim) :: x_grid

    ltau = lambdaMS*tau
    
    do i=1,s_dim       
       xs = xbar + s(i)
       value(i) = calc_u0(x_grid,r,xs,tau,x_dim) &
            + calc_u0(x_grid,xs,rp,tau,x_dim)
       value(i) = expm1(-value(i))
       value(i)=value(i)*sqrt(2.0_rk/pi/ltau**3)*&
            xs**2*exp(-s(i)**2/(2.0_rk*ltau))*&
            m0(r*xs/(2.0_rk*ltau))*m0(xs*rp/(2.0_rk*ltau))&
            /m0(r*rp/(4.0_rk*ltau))
    end do    

  end function SmallTau



  function LargeTau(s,r,rp,xbar,tau,s_dim,x_grid,x_dim) result(value)
    use math
    implicit none
    integer :: s_dim,x_dim,i
    real(kind=rk) :: r, rp, xbar,tau
    real(kind=rk) :: xs,ltau
    real(kind=rk), dimension(s_dim) :: s,value
    real(kind=rk), dimension(x_dim) :: x_grid

    ltau = lambdaMS*tau
    
    do i=1,s_dim       
       xs = xbar + s(i)
       value(i) = calc_u0(x_grid,r,xs,tau,x_dim) &
            + calc_u0(x_grid,xs,rp,tau,x_dim)
       value(i) = exp(-value(i))
       value(i)=value(i)*sqrt(2.0_rk/pi/ltau**3)*&
            xs**2*exp(-s(i)**2/(2.0_rk*ltau))*&
            m0(r*xs/(2.0_rk*ltau))*m0(xs*rp/(2.0_rk*ltau))&
            /m0(r*rp/(4.0_rk*ltau))
    end do

  end function LargeTau





  
  function TotalIntegrandSmallTau2(s,r,rp,xbar,tau,x_grid,x_dim) result(value)
    use math
    implicit none
    integer :: x_dim
    real(kind=rk) :: r, rp, xbar,tau
    real(kind=rk) :: xs, ltau
    real(kind=rk) :: s,value
    real(kind=rk), dimension(x_dim) :: x_grid
       
    ltau = lambdaMS*tau
    
    xs = xbar + s
    value = calc_u0(x_grid,r,xs,tau,x_dim) &
         + calc_u0(x_grid,xs,rp,tau,x_dim)
    value = expm1(-value)
    value = value*sqrt(2.0_rk/pi/ltau**3)*&
         xs**2*exp(-s**2/(2.0_rk*ltau))*&
         m0(r*xs/(2.0_rk*ltau))*m0(xs*rp/(2.0_rk*ltau))&
         /m0(r*rp/(4.0_rk*ltau))
    value = value*(C0coef(r,xs,tau)+C0coef(xs,rp,tau) &
         -2.0_rk*C0coef(r,rp,2.0_rk*tau)-calc_du0(x_grid,r,xs,tau,x_dim)&
         -calc_du0(x_grid,xs,rp,tau,x_dim))
    
  end function TotalIntegrandSmallTau2

  function TotalIntegrandLargerTau2(s,r,rp,xbar,tau,x_grid,x_dim) result(value)
    use math
    implicit none
    integer :: x_dim
    real(kind=rk) :: r, rp, xbar,tau
    real(kind=rk) :: xs, ltau
    real(kind=rk) :: s,value,t1,t2,t3
    real(kind=rk), dimension(x_dim) :: x_grid

    ltau = lambdaMS*tau
    
    xs = xbar + s
    value = calc_u0(x_grid,r,xs,tau,x_dim) &
         + calc_u0(x_grid,xs,rp,tau,x_dim)
    value = exp(-value)
    value= value*sqrt(2.0_rk/pi/ltau**3)*&
         xs**2*exp(-s**2/(2.0_rk*ltau))*&
         m0(r*xs/(2.0_rk*ltau))*m0(xs*rp/(2.0_rk*ltau))&
         /m0(r*rp/(4.0_rk*ltau))
    t1=r*xs/ltau/2
    t2=rp*xs/ltau/2
    t3=r*rp/ltau/4
    if (t1<10.0_rk .or. t2<10.0_rk .or. t3<10.0_rk) then
       value = value*(C0coef(r,xs,tau)+C0coef(xs,rp,tau) &
            -2.0_rk*C0coef(r,rp,2.0_rk*tau)-calc_du0(x_grid,r,xs,tau,x_dim)&
            -calc_du0(x_grid,xs,rp,tau,x_dim))
    else
       value = value*((s**2/ltau/2-0.5_rk)/tau &
            -calc_du0(x_grid,r,xs,tau,x_dim)&
            -calc_du0(x_grid,xs,rp,tau,x_dim))
    end if
    
    !value = value*(C0coef(r,xs,tau)+C0coef(xs,rp,tau) &
    !     -2.0_rk*C0coef(r,rp,2.0_rk*tau)-calc_du0(x_grid,r,xs,tau,x_dim)&
    !     -calc_du0(x_grid,xs,rp,tau,x_dim))
    !value = value*((s**2/ltau/2-0.5_rk)/tau &
    !     -calc_du0(x_grid,r,xs,tau,x_dim)&
    !     -calc_du0(x_grid,xs,rp,tau,x_dim))

  end function TotalIntegrandLargerTau2
  
  
  function SmallTau2(s,r,rp,xbar,tau,s_dim,x_grid,x_dim) result(value)
    use math
    implicit none
    integer :: s_dim,x_dim,i
    real(kind=rk) :: r, rp, xbar,tau
    real(kind=rk) :: xs,ltau
    real(kind=rk), dimension(s_dim) :: s,value
    real(kind=rk), dimension(x_dim) :: x_grid

    ltau = lambdaMS*tau
    
    do i=1,s_dim       
       xs = xbar + s(i)
       value(i) = calc_u0(x_grid,r,xs,tau,x_dim) &
            + calc_u0(x_grid,xs,rp,tau,x_dim)
       value(i) = expm1(-value(i))
       value(i)=value(i)*sqrt(2.0_rk/pi/ltau**3)*&
            xs**2*exp(-s(i)**2/(2.0_rk*ltau))*&
            m0(r*xs/(2.0_rk*ltau))*m0(xs*rp/(2.0_rk*ltau))&
            /m0(r*rp/(4.0_rk*ltau))
       value(i) = value(i)*(C0coef(r,xs,tau)+C0coef(xs,rp,tau) &
            -2.0_rk*C0coef(r,rp,2.0_rk*tau)-calc_du0(x_grid,r,xs,tau,x_dim)&
            -calc_du0(x_grid,xs,rp,tau,x_dim))
    end do


  end function SmallTau2



  function LargeTau2(s,r,rp,xbar,tau,s_dim,x_grid,x_dim) result(value)
    use math
    implicit none
    integer :: s_dim,x_dim,i
    real(kind=rk) :: r, rp, xbar,tau
    real(kind=rk) :: xs,ltau,t1,t2,t3
    real(kind=rk), dimension(s_dim) :: s,value
    real(kind=rk), dimension(x_dim) :: x_grid

    ltau = lambdaMS*tau
    
    do i=1,s_dim       
       xs = xbar + s(i)
       value(i) = calc_u0(x_grid,r,xs,tau,x_dim) &
            + calc_u0(x_grid,xs,rp,tau,x_dim)
       value(i) = exp(-value(i))
       value(i)=value(i)*sqrt(2.0_rk/pi/ltau**3)*&
            xs**2*exp(-s(i)**2/(2.0_rk*ltau))*&
            m0(r*xs/(2.0_rk*ltau))*m0(xs*rp/(2.0_rk*ltau))&
            /m0(r*rp/(4.0_rk*ltau))
       t1=r*xs/ltau/2
       t2=rp*xs/ltau/2
       t3=r*rp/ltau/4
       if (t1<10.0_rk .or. t2<10.0_rk .or. t3<10.0_rk) then
          value(i) = value(i)*(C0coef(r,xs,tau)+C0coef(xs,rp,tau) &
               -2.0_rk*C0coef(r,rp,2.0_rk*tau)-calc_du0(x_grid,r,xs,tau,x_dim)&
               -calc_du0(x_grid,xs,rp,tau,x_dim))
       else
          value(i) = value(i)*((s(i)**2/ltau/2-0.5_rk)/tau &
               -calc_du0(x_grid,r,xs,tau,x_dim)&
               -calc_du0(x_grid,xs,rp,tau,x_dim))
       end if

    end do


  end function LargeTau2





  function calc_u0_2(x_grid,r,rp,tau,x_dim) result(value)
    use math
    implicit none
    integer :: x_dim
    real(kind=rk), dimension(x_dim) :: x_grid
    real(kind=rk) :: r,rp,tau, value,x,y

    if (r<=x_grid(x_dim) .and. rp<=x_grid(x_dim)) then

       if (r<x_grid(1)) then
          x = x_grid(1)
       else
          x = r
       end if
       
       if (rp<x_grid(1)) then
          y = x_grid(1)
       else
          y = rp
       end if
       
       value = interpolate2dc(x,y,1)
       
    else
       if (r>x_grid(x_dim) .and. rp>x_grid(x_dim)) then
          value =  Uij(x_dim,x_dim,1)*&
               u0start(r,rp,tau)/&
               u0start(x_grid(x_dim),x_grid(x_dim),tau)
       elseif (r>x_grid(x_dim) .and. rp<=x_grid(x_dim)) then
          if (rp<x_grid(1)) then
             y = x_grid(1)
          else
             y = rp
          end if
          value = interpolate2dc(x_grid(x_dim),y,1)*&
               u0start(r,y,tau)/&
               u0start(x_grid(x_dim),y,tau)
       else
          if (r<x_grid(1)) then
             x = x_grid(1)
          else
             x = r
          end if
          value = interpolate2dc(x,x_grid(x_dim),1)*&
               u0start(x,rp,tau)/&
               u0start(x,x_grid(x_dim),tau)
       end if
       
    end if
    
    
  end function calc_u0_2

  function calc_u0(x_grid,r,rp,tau,x_dim) result(value)
    use math
    implicit none
    integer :: x_dim
    real(kind=rk), dimension(x_dim) :: x_grid
    real(kind=rk) :: r,rp,tau, value,x,y

    if (r<=x_grid(x_dim) .and. rp<=x_grid(x_dim)) then

       if (r<x_grid(1)) then
          x = x_grid(1)
       else
          x = r
       end if
       
       if (rp<x_grid(1)) then
          y = x_grid(1)
       else
          y = rp
       end if
       
       value = interpolate2dc(x,y,1)
       
    else
       if (r>x_grid(x_dim) .and. rp>x_grid(x_dim)) then
          value =  Uij(x_dim,x_dim,1)*&
               u0ep(r,rp,tau,x_grid,x_dim)/&
               u0ep(x_grid(x_dim),x_grid(x_dim),tau,x_grid,x_dim)
       elseif (r>x_grid(x_dim) .and. rp<=x_grid(x_dim)) then
          if (rp<x_grid(1)) then
             y = x_grid(1)
          else
             y = rp
          end if
          value = interpolate2dc(x_grid(x_dim),y,1)*&
               u0ep(r,y,tau,x_grid,x_dim)/&
               u0ep(x_grid(x_dim),y,tau,x_grid,x_dim)
       else
          if (r<x_grid(1)) then
             x = x_grid(1)
          else
             x = r
          end if
          value = interpolate2dc(x,x_grid(x_dim),1)*&
               u0ep(x,rp,tau,x_grid,x_dim)/&
               u0ep(x,x_grid(x_dim),tau,x_grid,x_dim)
       end if
       
    end if
    
    
  end function calc_u0


  function calc_du0(x_grid,r,rp,tau,x_dim) result(value)
    use math
    implicit none
    integer :: x_dim
    real(kind=rk), dimension(x_dim) :: x_grid
    real(kind=rk) :: r,rp,tau, value,x,y

    if (r<=x_grid(x_dim) .and. rp<=x_grid(x_dim)) then

       if (r<x_grid(1)) then
          x = x_grid(1)
       else
          x = r
       end if
       
       if (rp<x_grid(1)) then
          y = x_grid(1)
       else
          y = rp
       end if
       
       value = interpolate2dc(x,y,2)
       
    else

       if (r>x_grid(x_dim) .and. rp>x_grid(x_dim)) then
          value =  Uij(x_dim,x_dim,1)*&
               u0ep(r,rp,tau,x_grid,x_dim)/tau/&
               u0ep(x_grid(x_dim),x_grid(x_dim),tau,x_grid,x_dim) &
               + u0ep(r,rp,tau,x_grid,x_dim)/&
               u0ep(x_grid(x_dim),x_grid(x_dim),tau,x_grid,x_dim)* &
               Uij(x_dim,x_dim,2) &
               - u0ep(r,rp,tau,x_grid,x_dim)*Uij(x_dim,x_dim,1)/ &
               u0ep(x_grid(x_dim),x_grid(x_dim),tau,x_grid,x_dim)**2* &
               Uij(x_dim,x_dim,2) 
          
       elseif (r>x_grid(x_dim) .and. rp<=x_grid(x_dim)) then
          if (rp<x_grid(1)) then
             y = x_grid(1)
          else
             y = rp
          end if
          !value = interpolate2dc(x_grid(x_dim),y,2)*&
          !     u0start(r,y,tau)/&
          !     u0start(x_grid(x_dim),y,tau)

          value = interpolate2dc(x_grid(x_dim),y,1)*&
               u0ep(r,y,tau,x_grid,x_dim)/tau/&
               u0ep(x_grid(x_dim),y,tau,x_grid,x_dim) &
               + u0ep(r,y,tau,x_grid,x_dim)/&
               u0ep(x_grid(x_dim),y,tau,x_grid,x_dim)* &
               interpolate2dc(x_grid(x_dim),y,2) &
               - u0ep(r,y,tau,x_grid,x_dim)*interpolate2dc(x_grid(x_dim),y,1)/ &
               u0ep(x_grid(x_dim),y,tau,x_grid,x_dim)**2* &
               interpolate2dc(x_grid(x_dim),y,2)
          
       else
          if (r<x_grid(1)) then
             x = x_grid(1)
          else
             x = r
          end if
          !value = interpolate2dc(x,x_grid(x_dim),2)*&
          !     u0start(x,rp,tau)/&
          !     u0start(x,x_grid(x_dim),tau)

          value = interpolate2dc(x,x_grid(x_dim),1)*&
               u0ep(x,rp,tau,x_grid,x_dim)/tau/&
               u0ep(x,x_grid(x_dim),tau,x_grid,x_dim) &
               + u0ep(x,rp,tau,x_grid,x_dim)/&
               u0ep(x,x_grid(x_dim),tau,x_grid,x_dim)* &
               interpolate2dc(x,x_grid(x_dim),2) &
               - u0ep(x,rp,tau,x_grid,x_dim)*interpolate2dc(x,x_grid(x_dim),1)/ &
               u0ep(x,x_grid(x_dim),tau,x_grid,x_dim)**2* &
               interpolate2dc(x,x_grid(x_dim),2)
       end if
       
    end if
    
    
  end function calc_du0
  
  function calc_du0_2(x_grid,r,rp,tau,x_dim) result(value)
    use math
    implicit none
    integer :: x_dim
    real(kind=rk), dimension(x_dim) :: x_grid
    real(kind=rk) :: r,rp,tau, value,x,y

    if (r<=x_grid(x_dim) .and. rp<=x_grid(x_dim)) then

       if (r<x_grid(1)) then
          x = x_grid(1)
       else
          x = r
       end if
       
       if (rp<x_grid(1)) then
          y = x_grid(1)
       else
          y = rp
       end if
       
       value = interpolate2dc(x,y,2)
       
    else

       if (r>x_grid(x_dim) .and. rp>x_grid(x_dim)) then
          value =  Uij(x_dim,x_dim,1)*&
               u0start(r,rp,tau)/tau/&
               u0start(x_grid(x_dim),x_grid(x_dim),tau) &
               + u0start(r,rp,tau)/&
               u0start(x_grid(x_dim),x_grid(x_dim),tau)* &
               Uij(x_dim,x_dim,2) &
               - u0start(r,rp,tau)*Uij(x_dim,x_dim,1)/ &
               u0start(x_grid(x_dim),x_grid(x_dim),tau)**2* &
               Uij(x_dim,x_dim,2) 
          
       elseif (r>x_grid(x_dim) .and. rp<=x_grid(x_dim)) then
          if (rp<x_grid(1)) then
             y = x_grid(1)
          else
             y = rp
          end if
          !value = interpolate2dc(x_grid(x_dim),y,2)*&
          !     u0start(r,y,tau)/&
          !     u0start(x_grid(x_dim),y,tau)

          value = interpolate2dc(x_grid(x_dim),y,1)*&
               u0start(r,y,tau)/tau/&
               u0start(x_grid(x_dim),y,tau) &
               + u0start(r,y,tau)/&
               u0start(x_grid(x_dim),y,tau)* &
               interpolate2dc(x_grid(x_dim),y,2) &
               - u0start(r,y,tau)*interpolate2dc(x_grid(x_dim),y,1)/ &
               u0start(x_grid(x_dim),y,tau)**2* &
               interpolate2dc(x_grid(x_dim),y,2)
          
       else
          if (r<x_grid(1)) then
             x = x_grid(1)
          else
             x = r
          end if
          !value = interpolate2dc(x,x_grid(x_dim),2)*&
          !     u0start(x,rp,tau)/&
          !     u0start(x,x_grid(x_dim),tau)

          value = interpolate2dc(x,x_grid(x_dim),1)*&
               u0start(x,rp,tau)/tau/&
               u0start(x,x_grid(x_dim),tau) &
               + u0start(x,rp,tau)/&
               u0start(x,x_grid(x_dim),tau)* &
               interpolate2dc(x,x_grid(x_dim),2) &
               - u0start(x,rp,tau)*interpolate2dc(x,x_grid(x_dim),1)/ &
               u0start(x,x_grid(x_dim),tau)**2* &
               interpolate2dc(x,x_grid(x_dim),2)
       end if
       
    end if
    
    
  end function calc_du0_2




  subroutine u0_HighTempApprox2(r,u0,tau,r_dim)
    use math
    implicit none
    integer, intent(in) :: r_dim
    real(kind = rk), intent(in) :: tau
    real(kind = rk), dimension(r_dim), intent(in) :: r
    real(kind = rk), dimension(r_dim,r_dim), intent(inout) :: u0
    integer :: i,j
    
    do i = 1,r_dim
       do j = i,r_dim

          u0(i,j) = u0start(r(i),r(j),tau)
          u0(j,i) = u0(i,j)
          
       end do
    end do

  end subroutine u0_HighTempApprox2

  subroutine u0_HighTempApprox(r,u0,tau,r_dim)
    use math
    implicit none
    integer, intent(in) :: r_dim
    real(kind = rk), intent(in) :: tau
    real(kind = rk), dimension(r_dim), intent(in) :: r
    real(kind = rk), dimension(r_dim,r_dim), intent(inout) :: u0
    real(kind = rk) :: x,y
    integer :: i,j

    ! r(1) = 0.0_rk
    do i = 1,r_dim
       do j = i,r_dim
          x = maxval((/r(i),0.99_rk*minval((/r(2),r(r_dim-1)/))/))
          y = maxval((/r(j),0.99_rk*minval((/r(2),r(r_dim-1)/))/))

          !x = r(i)
          !y = r(j)

          u0(i,j) = u0start(x,y,tau)
          u0(j,i) = u0(i,j)
          
       end do
    end do

  end subroutine u0_HighTempApprox

  function Uaction(r1,r2,tau) result(U)
    use math
    implicit none
    real(kind = rk) :: r1,r2,U,tau !,dr    

    ! primitive
    U = 0.5_rk*tau*ZZMS*(1.0_rk/maxval((/r1,1.0e-6_rk/)) &
         + 1.0_rk/maxval((/r2,1.0e-6_rk/)))
        
  end function Uaction

  function u0ep(r,rp,tau,x_grid,x_dim) result(value)
    use math
    implicit none
    integer :: x_dim
    real(kind = rk) :: r,rp,tau,value !,x,y    
    real(kind=rk), dimension(x_dim) :: x_grid


    !value = Uaction(r,rp,tau)
    !value = u0start2(r,rp,tau)
    !return

    if (r>x_grid(x_dim) .and. rp>x_grid(x_dim)) then
       value = Uaction(r,rp,tau)
    else if (r<=x_grid(x_dim) .and. rp>x_grid(x_dim))then
       value = (interpolate2dc(r,r,1)+&
            interpolate2dc(x_grid(x_dim),x_grid(x_dim),1)/&
            Uaction(x_grid(x_dim),x_grid(x_dim),tau)*&
            Uaction(rp,rp,tau))/2
    else
       value = (interpolate2dc(rp,rp,1)+&
            interpolate2dc(x_grid(x_dim),x_grid(x_dim),1)/&
            Uaction(x_grid(x_dim),x_grid(x_dim),tau)*&
            Uaction(r,r,tau))/2
    end if
       
  end function u0ep

  
  function u0start(r,rp,tau) result(value)
    use math
    implicit none
    real(kind = rk) :: r,rp,tau,value !,x,y    

    !value = Uaction(r,rp,tau)
    !value = u0start2(r,rp,tau)
    !return

    !if (r<1.0e-10_rk .and. rp<1.0e-10_rk) then
    !   value =  NearOrigin(r,rp,tau)
       !elseif (r<1.0e-7 .and. rp>=1.0e-7) then
       !   value =  0.5_rk*NearOrigin(r,r,tau)+&
       !        0.5_rk*Uaction(rp,rp,tau)
       !elseif (rp<1.0e-7 .and. r>=1.0e-7) then
       !   value =  0.5_rk*NearOrigin(rp,rp,tau)+&
       !        0.5_rk*Uaction(r,r,tau)
    !else
    value = Uaction(r,rp,tau)
    !end if

  end function u0start

  function Uaction2(r1,r2,tau) result(U)
    use math
    implicit none
    real(kind = rk) :: r1,r2,U,tau,dr    
    
    ! semi-classical
    U = 0.0_rk
    dr = r2-r1
    if (abs(dr)>1.0e-10_rk) then
       U = tau*ZZMS*log(r2/r1)/dr
    else
       U = tau*ZZMS/r1
    end if
    
  end function Uaction2

  
  function u0start2(r,rp,tau) result(value)
    use math
    implicit none
    real(kind = rk) :: r,rp,tau,value !,x,y    
    
    if (r<1.0e-7 .and. rp<1.0e-7) then
       value =  NearOrigin(r,rp,tau)
    elseif (r<1.0e-7 .and. rp>=1.0e-7) then
       value =  0.5_rk*NearOrigin(r,r,tau)+&
            0.5_rk*Uaction2(rp,rp,tau)
    elseif (rp<1.0e-7 .and. r>=1.0e-7) then
       value =  0.5_rk*NearOrigin(rp,rp,tau)+&
            0.5_rk*Uaction2(r,r,tau)
    else
       value = Uaction2(r,rp,tau)
    end if

  end function u0start2



  function NearOrigin(r,rp,tau) result(V)
    use math
    implicit none
    real(kind = rk) :: r,rp,tau,V
    real(kind = rk), dimension(7) :: PP
    real(kind = rk) :: ii, vakio, vakio2
    integer :: i
    
    ! Let's also compute the known value at the origin
    PP = (/1.772453851_rk, -0.074137740_rk, 0.005834805_rk, -0.000382686_rk, &
         0.000008738_rk, 0.000002138_rk, -0.000000356_rk/)
    
    V = 0.0_rk
    vakio = (ZZMS)**2/lambdaMS*tau
    vakio2 = ZZMS/abs(ZZMS)
    do i = 1,7
       ii = real(i,kind=rk)
       V = V + vakio2**i*PP(i)*vakio**(ii/2.0_rk)
    end do
    V = V - ZZMS/2/lambdaMS*(r+rp)
    
    !write(*,*) 'value at origin', V

  end function NearOrigin


  function PairPotential(r,u0,tau,r_dim) result(u)
    !Compute u = u0 - log(1+rho0/rhoK*f), 
    !where f = 1/(4*pi*(x-y))*(d_x-d_y)u0(x,y)
    use math
    implicit none
    integer :: r_dim,i,j
    real(kind = rk), dimension(r_dim) :: r
    real(kind = rk), dimension(r_dim,r_dim) :: u0,u
    real(kind = rk) :: tau, twoltau, fourpi,value
    real(kind=rk) :: m1_8pi,m1_4pi,du0
    
    m1_8pi = 1.0_rk/pi/8
    m1_4pi = 1.0_rk/pi/4
    
    twoltau = 2.0_rk*lambdaMS*tau
    fourpi = 4.0_rk*pi

    do i=1,r_dim
       do j = i,r_dim
          
          if (i==j) then             
             du0 = m1_8pi*(d2Ux(r(i),r(i),1) + d2Uy(r(i),r(i),1) &
                  - 2.0_rk*d2Uxy(r(i),r(i),1))              
          else
             du0 = m1_4pi*(dUx(r(i),r(j),1)-dUy(r(i),r(j),1)) &
                  /(r(i)-r(j))
          end if
          
          value = fourpi*r(i)*r(j)*m0(r(i)*r(j)/twoltau)*du0          
          !u(i,j) = u0(i,j) - log(1.0_rk+value)
          u(i,j) = u0(i,j) - log1p(value)
          
          u(j,i) = u(i,j)
          
       end do
    end do
    
    
  end function PairPotential



  subroutine PairPotentials(r,u,du,tau,r_dim)
    !Compute u = u0 - log(1+rho0/rhoK*f), 
    !where f = 1/(4*pi*(x-y))*(d_x-d_y)u0(x,y)
    !and du/dbeta
    use math
    implicit none
    integer, intent(in) :: r_dim
    integer :: i,j
    real(kind = rk), dimension(r_dim), intent(in) :: r
    real(kind = rk), dimension(r_dim,r_dim), intent(out) :: u,du
    real(kind = rk), intent(in) :: tau
    real(kind = rk) :: twoltau, fourpi,value, CrelK, value1
    real(kind=rk) :: m1_8pi,m1_4pi,u_t,du_t
    
    m1_8pi = 1.0_rk/pi/8
    m1_4pi = 1.0_rk/pi/4
    
    twoltau = 2.0_rk*lambdaMS*tau
    fourpi = 4.0_rk*pi

    do i=1,r_dim
       do j = i,r_dim
          
          if (i==j) then
             u_t = m1_8pi*(d2Ux(r(i),r(i),1) + d2Uy(r(i),r(i),1) &
                  - 2.0_rk*d2Uxy(r(i),r(i),1)) 

             du_t = m1_8pi*(d2Ux(r(i),r(i),2) + d2Uy(r(i),r(i),2) &
                  - 2.0_rk*d2Uxy(r(i),r(i),2)) 
          else
             u_t = m1_4pi*(dUx(r(i),r(j),1)-dUy(r(i),r(j),1)) &
                  /(r(i)-r(j))

             du_t = m1_4pi*(dUx(r(i),r(j),2)-dUy(r(i),r(j),2)) &
                  /(r(i)-r(j))
          end if
          
          
          value1 = fourpi*r(i)*r(j)*m0(r(i)*r(j)/twoltau)
          value = value1*u_t          
          !u(i,j) = u0(i,j) - log(1.0_rk+value)
          u(i,j) = Uij(i,j,1) - log1p(value)
          u(j,i) = u(i,j)
          
          CrelK = CrelKcoef(r(i),r(j),tau)
          
          du_t = 1.0_rk*&
               ( (1.0_rk + value)*Uij(i,j,2)-CrelK-&
               value*C0coef(r(i),r(j),tau)-&
               value1*du_t )
          
          du(i,j) = CrelK + exp(u(i,j)-Uij(i,j,1))*du_t
          du(j,i) = du(i,j)                    
          

       end do
    end do
    
    
  end subroutine PairPotentials


  subroutine PairPotentials2(r,u,du,tau,r_dim)
    !Compute u = u0 - log(1+rho0/rhoK*f), 
    !where f = 1/(4*pi*(x-y))*(d_x-d_y)u0(x,y)
    !and du/dbeta
    use math
    implicit none
    integer, intent(in) :: r_dim
    integer :: i,j
    real(kind = rk), dimension(r_dim), intent(in) :: r
    real(kind = rk), dimension(r_dim,r_dim), intent(out) :: u,du
    real(kind = rk), intent(in) :: tau
    real(kind = rk) :: twoltau, fourpi,value,z, value1
    real(kind=rk) :: m1_8pi,m1_4pi,u_t,du_t
    
    m1_8pi = 1.0_rk/pi/8
    m1_4pi = 1.0_rk/pi/4
    
    twoltau = 2.0_rk*lambdaMS*tau
    fourpi = 4.0_rk*pi

    do i=1,r_dim
       do j = i,r_dim
          
          if (i==j) then
             u_t = m1_8pi*(d2Ux(r(i),r(i),1) + d2Uy(r(i),r(i),1) &
                  - 2.0_rk*d2Uxy(r(i),r(i),1)) 

             du_t = m1_8pi*(d2Ux(r(i),r(i),2) + d2Uy(r(i),r(i),2) &
                  - 2.0_rk*d2Uxy(r(i),r(i),2)) 
          else
             u_t = m1_4pi*(dUx(r(i),r(j),1)-dUy(r(i),r(j),1)) &
                  /(r(i)-r(j))

             du_t = m1_4pi*(dUx(r(i),r(j),2)-dUy(r(i),r(j),2)) &
                  /(r(i)-r(j))
          end if
          
          
          value1 = fourpi*r(i)*r(j)*m0(r(i)*r(j)/twoltau)
          value = value1*u_t          
          !u(i,j) = u0(i,j) - log(1.0_rk+value)
          u(i,j) = Uij(i,j,1) - log1p(value)
          u(j,i) = u(i,j)                    
          
          !du(i,j) = Uij(i,j,2)-&
          !     (value*(C0coef(r(i),r(j),tau)-CrelKcoef(r(i),r(j),tau)) &
          !     +value1*du_t)/(1.0_rk + value)
          
          z = r(i)*r(j)/twoltau

          du(i,j) = Uij(i,j,2)-&
               (value1*du_t - &
               (fourpi*r(i)*r(j)*dm0(z)*u_t - &
               value) &
               *z/tau) &
               /(1.0_rk + value)
          du(j,i) = du(i,j)                    
          
       end do
    end do
    
    
  end subroutine PairPotentials2

  
  subroutine PairPotentials2SRC(r,u,du,tau,r_dim,pot)
    !Compute u = u0 - log(1+rho0/rhoK*f),                                                   
    !where f = 1/(4*pi*(x-y))*(d_x-d_y)u0(x,y)                                              
    !and du/dbeta
    !
    !with short range correction in case of Ewald sums
    use math
    implicit none
    integer, intent(in) :: r_dim,pot
    integer :: i,j
    real(kind = rk), dimension(r_dim), intent(in) :: r
    !real(kind = rk), dimension(ulr_dim), intent(in) :: ulr,dulr,d2ulr
    real(kind = rk), dimension(r_dim,r_dim), intent(out) :: u,du
    real(kind = rk), intent(in) :: tau
    real(kind = rk) :: twoltau, fourpi,value,z, value1
    real(kind=rk) :: m1_8pi,m1_4pi,u_t,du_t,ulrt,dulrt
    
    m1_8pi = 1.0_rk/pi/8
    m1_4pi = 1.0_rk/pi/4
    
    twoltau = 2.0_rk*lambdaMS*tau
    fourpi = 4.0_rk*pi
    
    do i=1,r_dim
       do j = i,r_dim
          if (i==j) then
             u_t = m1_8pi*(d2Ux(r(i),r(i),1) + d2Uy(r(i),r(i),1) &
                  - 2.0_rk*d2Uxy(r(i),r(i),1))

             ulrt = m1_8pi*(evald2LR(r(i),pot,1) + evald2LR(r(i),pot,1))/2
                 ! - 2.0_rk*d2Uxy(r(i),r(i),1)) = 0                                        
             du_t = m1_8pi*(d2Ux(r(i),r(i),2) + d2Uy(r(i),r(i),2) &
                  - 2.0_rk*d2Uxy(r(i),r(i),2))
             
             dulrt = m1_8pi*(evald2LR(r(i),pot,2) + evald2LR(r(i),pot,2))/2


          else
             u_t = m1_4pi*(dUx(r(i),r(j),1)-dUy(r(i),r(j),1)) &
                  /(r(i)-r(j))

             ulrt = m1_4pi*(evaldLR(r(i),pot,1)-evaldLR(r(j),pot,1))/2 &
                  /(r(i)-r(j))

             du_t = m1_4pi*(dUx(r(i),r(j),2)-dUy(r(i),r(j),2)) &
                  /(r(i)-r(j))

             dulrt = m1_4pi*(evaldLR(r(i),pot,2)-evaldLR(r(j),pot,2))/2 &
                  /(r(i)-r(j))


          end if

          u_t = u_t - ulrt
          du_t = du_t - dulrt
          value1 = fourpi*r(i)*r(j)*m0(r(i)*r(j)/twoltau)
          value = value1*u_t
          !u(i,j) = u0(i,j) - log(1.0_rk+value)                                             
          !u(i,j) = Uij(i,j,1)-(ulr(i)+ulr(j))/2*tau*ZZMS &
          !     - log1p(value)
          u(i,j) = Uij(i,j,1)-log1p(value)
          u(j,i) = u(i,j)

          !du(i,j) = Uij(i,j,2)-&                                                           
          !     (value*(C0coef(r(i),r(j),tau)-CrelKcoef(r(i),r(j),tau)) &                   
          !     +value1*du_t)/(1.0_rk + value)                                              

          z = r(i)*r(j)/twoltau

          !du(i,j) = Uij(i,j,2) - (ulr(r(i))+ulr(r(j)))/2*ZZMS &
          du(i,j) = Uij(i,j,2) &
               - (value1*du_t - &
               (fourpi*r(i)*r(j)*dm0(z)*u_t - &
               value) &
               *z/tau) &
               /(1.0_rk + value)
          du(j,i) = du(i,j)
          
       end do
    end do


  end subroutine PairPotentials2SRC






  !---- for the calculation of the beta-derivative ----
  
  function zdm0_m0(z) result(value)
    use math
    implicit none
    real(kind = rk) :: z,value

    if (z<10.0_rk) then
       value = z*dm0(z)/m0(z)
    else
       value = z-1.0_rk
    end if


  end function zdm0_m0



  function C0coef(r,rp,tau) result(value)
    use math
    implicit none
    real(kind=rk) :: r, rp,tau,z
    real(kind=rk) :: ltau,value
    
    ltau = lambdaMS*tau
    z = r*rp/ltau/2
    
    if (z<10.0_rk) then
       value = ( (r**2+rp**2)/ltau/4 - 1.5_rk - z*dm0(z)/m0(z) )/tau
    else
       value = ( (r**2+rp**2)/ltau/4 - 1.5_rk - (z-1.0_rk) )/tau
    end if

  end function C0coef
  
  function CrelKcoef(r,rp,tau) result(value)
    use math
    implicit none
    real(kind=rk) :: r, rp,tau
    real(kind=rk) :: value
    
    value = ( (r-rp)**2/lambdaMS/tau/4 - 1.5_rk)/tau

  end function CrelKcoef


end module matrixsquaring
