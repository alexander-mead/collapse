PROGRAM collapse

  IMPLICIT NONE
  REAL :: om_m, om_w, om_nu, om_c, om_v
  REAL :: am, dm, w0, wm, mua, H0rc, wa, astar, nblip
  REAL :: dinit, ainit, amin, amax, vinit, ac, dc, dmin, dmax
  REAL :: av, a_rmax, d_rmax, Dv, rmax, rv, g, f, h
  REAL, ALLOCATABLE :: d(:), a(:), v(:)
  REAL, ALLOCATABLE :: dnl(:), vnl(:)
  REAL, ALLOCATABLE :: rl(:), rnl(:)
  REAL, ALLOCATABLE :: a_coll(:), r_coll(:)
  REAL, ALLOCATABLE :: growth(:), a_growth(:), rate(:), a_rate(:), bigG(:)
  INTEGER :: i, j, n, m, k, k2, na
  INTEGER :: iw, img
  INTEGER :: icol

  !This is a very basic and crude spherical-model code written by Alexander Mead
  !There must be much cleverer ways of doing this calculation,
  !but since speed and accuracy were not pressing issues for me I didn't worry about making it faster or cleverer
  !Comparing the EdS results to the exact solution I guess that the result is accurate to
  !4 s.f. for the absolute value of \delta_c and 3 s.f. for \Delta_v

  !Note that the user does not select the value of 'a' they require for \delta_c and \Delta_v
  !instead the code produces a table of values and the collapse time cannot be predetermined
  !Interpolation should be used for intermediate values 

  WRITE(*,*)
  WRITE(*,*) 'Spherical collapse integrator by Alexander Mead'
  WRITE(*,*) '==============================================='
  WRITE(*,*)

  !Print out the collapse data or not
  !0 - No
  !1 - Yes
  icol=0

  !Assign the cosmological model
  CALL assign_cosmology

  !Write out the cosmological model
  !a, Omega_i(a) and w(a) for DE
  WRITE(*,*) 'Writing cosmology.dat'
  OPEN(7,file='cosmology.dat')
  amin=0.001
  amax=10
  na=200
  CALL fill_table(amin,amax,a,na,1)
  DO i=1,na
     WRITE(7,*) a(i), omega_m(a(i)), omega_v(a(i)), omega_w(a(i)), w(a(i))
  END DO
  CLOSE(7)
  DEALLOCATE(a)
  WRITE(*,*) 'Done'
  WRITE(*,*)
  
  !Number of collapse 'a' to calculate
  !Actually you get fewer than this because some perturbations do not collapse due to dark energy
  !Set the minimum and maximum values. d_min=1e-7 and d_max=1e-3 are good, both very much within EdS regime
  dmin=1e-7
  dmax=1e-3
  m=200

  !BCs for integration. Note ainit=dinit means that collapse should occur around a=1 for dmin
  !amax should be slightly greater than 1 to ensure at least a few points for a>0.9 (i.e not to miss out a=1)
  !vinit=1 is EdS growing mode solution if dinit=ainit
  ainit=dmin
  amax=3.
  vinit=1.*(dmin/ainit)

  !Do the linear growth factor calculation
  WRITE(*,*) 'Solving growth factor ODE'
  CALL ode_adaptive(growth,v,a_growth,ainit,amax,dmin,vinit,fd,fv,0.0001,3,1)
  DEALLOCATE(v)
  !Do the linear logarithmic-growth rate calculation
  WRITE(*,*) 'Solving growth rate ODE'
  CALL ode_adaptive(rate,v,a_rate,dmin,amax,1.,1.,ff,one,0.0001,3,1)
  DEALLOCATE(v)

  !Table integration to calculate G(a)=int_0^a g(a')/a' da'
  ALLOCATE(bigG(SIZE(growth)))
  bigG=0.
  DO i=1,SIZE(growth)
     !Do the integral up to table position i
     bigG(i)=inttab_range(a_growth,growth/a_growth,SIZE(growth),1,i,1)
     !Add on the section that is missing from the beginning (NB. g(a=0)/0. = 1, so you just add on a rectangle of height g*a/a=g)
     bigG(i)=bigG(i)+growth(1)
  END DO

  !Number of points for ODE calculations (needs to be large to capture final stages of collapse nicely)
  !I tried doing a few clever things than this, but none worked well enough to beat the ignornant method
  n=100000

  !Now do the spherical-collapse calculation
  !WRITE(*,fmt='(I5,7F10.4)') j, ac, omega_m(ac), dc, Dv, g, f, h
  WRITE(*,*) '   i         a      Om_m        dc        Dv         g         f         G'
  WRITE(*,*) '=========================================================================='
  OPEN(8,file='dcDv.dat')
  DO j=1,m

     !Loop over a log range for delta_init
     dinit=exp(log(dmin)+log(dmax/dmin)*float(j-1)/float(m-1))

     !Solve the ODEs for both the linear and non-linear growth
     !Do both integrations with the same 'a' range (ainit->amax) and using the same number of time steps
     !This means that arrays 'a', and 'anl' will be identical, which simplifies the later calculation
     CALL ode_crass(dnl,vnl,a,ainit,amax,dinit,vinit,fd,fvnl,n,3,1)
     DEALLOCATE(a)
     CALL ode_crass(d,v,a,ainit,amax,dinit,vinit,fd,fv,n,3,1)

     !If this condtion is met then collapse occured some time a<amax
     IF(dnl(n)==0.) THEN

        !Find position in array when dnl(k)=0
        DO i=1,n
           IF(dnl(i)==0.) THEN
              !k is the new maxium size of the arrays
              k=i-1
              EXIT
           END IF
        END DO

        !Cut away parts of the arrays for a>ac
        CALL amputate(a,k)
        CALL amputate(d,k)
        CALL amputate(dnl,k)

        !Allocate arrays for the perturbation radius
        ALLOCATE(rl(k),rnl(k))

        !Relate the radius to density (missing some prefactors)
        rl=(1.+d)**(-1./3.)
        rnl=(1.+dnl)**(-1./3.)

        !Convert to physical radius from comoving
        rl=a*rl
        rnl=a*rnl

        IF(icol==1) THEN
           !Print out full calculation first time round if icol=1 above
           OPEN(7,file='collapse.dat')
           DO i=1,k
              WRITE(7,*) a(i), d(i), dnl(i), rl(i), rnl(i)
           END DO
           CLOSE(7)
           !So it only gets printed once
           icol=0 
        END IF

        !Now find the collapse point (very crude)
        !More accurate calculations seem to be worse
        !I think this is due to the fact that delta spikes very quickly
        !I tried changing to radius and calculating r->0 instead, but worked less well
        !Also tried using different variables for the ODEs, which also worked less well
                
        !Collapse has occured so use previous a as ac and d as dc
        !This is slightly crude, could probably interpolate.
        !I make up for laziness by using a large number of points for the ODEs
        ac=a(k)
        dc=d(k)

        !Now to Delta_v calculation

        !Find the 'a' value when the perturbation is maximum size
        !I could imagine models where the perturbation has mulitple maxima, so this might be correct for these models
        a_rmax=maximum(a,rnl,k)

        !Find the over-density at this point
        d_rmax=exp(find(log(a_rmax),log(a),log(dnl),k,1,3))

        !Find the maximum radius
        rmax=find(log(a_rmax),log(a),rnl,k,1,3)

        !The radius of the perturbation when it is virialised is half maximum is the EdS condition
        rv=rmax/2.

        !Need to assign new arrays for the collapse branch of r such that it is monotonic
        !There could be trouble here in DE models where the perurbation may stop growing and then start growing again
        k2=int_split(d_rmax,dnl,SIZE(dnl))

        !Allocate collapse branch arrays
        ALLOCATE(a_coll(k-k2+1),r_coll(k-k2+1))

        !Fill collapse branch arrays
        DO i=k2,k
           a_coll(i-k2+1)=a(i)
           r_coll(i-k2+1)=rnl(i)
        END DO

        !Find the scale factor when the perturbation has reached virial radius
        av=exp(find(rv,r_coll,log(a_coll),SIZE(r_coll),3,3))

        !Deallocate collapse branch arrays
        DEALLOCATE(a_coll,r_coll)

        !Spherical model approximation is that perturbation is at virial radius when
        !'collapse' is considered to have occured, the time of which has already been calculated
        !So Dv can now be found
        Dv=exp(find(log(av),log(a),log(dnl),k,1,3))*(ac/av)**3.
        Dv=Dv+1.

        !Get the values of the growth function, log-growth rate and integrated growth
        g=find(ac,a_growth,growth,SIZE(a_growth),3,3)
        f=find(ac,a_rate,rate,SIZE(a_growth),3,3)
        h=find(ac,a_growth,bigG,SIZE(a_growth),3,3)

        !Write out results
        WRITE(*,fmt='(I5,7F10.4)') j, ac, omega_m(ac), dc, Dv, g, f, h
        WRITE(8,*) ac, omega_m(ac), dc, Dv, g, f, h

        !Deallocate the radius arrays ready for the next calculation
        DEALLOCATE(rl,rnl)

     END IF

     !Deallocate arrays ready for next calculation
     DEALLOCATE(d,v,a)
     DEALLOCATE(dnl,vnl)

  END DO
  CLOSE(8)

CONTAINS

  FUNCTION maximum(x,y,n)

    IMPLICIT NONE
    REAL :: maximum
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: x(n), y(n)
    REAL :: x1, x2, x3, y1, y2, y3, a, b, c
    INTEGER :: i

    !From an array y(x) finds the first maxima

    DO i=1,n-1
       IF(y(i+1)<y(i)) THEN

          x1=x(i-1)
          x2=x(i)
          x3=x(i+1)

          y1=y(i-1)
          y2=y(i)
          y3=y(i+1)

          !Fit a quadratic around the maximum
          CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)

          !Read off x_max
          maximum=-b/(2.*a)

          EXIT
       !ELSE IF(i==n-1) THEN
       !   STOP 'MAXIMUM: Error, maximum not found'
       END IF
    END DO

  END FUNCTION maximum

  SUBROUTINE amputate(arr,n_new)

    IMPLICIT NONE
    REAL, ALLOCATABLE, INTENT(INOUT) :: arr(:)
    REAL, ALLOCATABLE :: hold(:)
    INTEGER, INTENT(IN) :: n_new
    INTEGER :: n_old, i

    !Chop an array down to a smaller size

    n_old=SIZE(arr)

    IF(n_old<n_new) STOP 'AMPUTATE: Error, new array should be smaller than the old one'

    ALLOCATE(hold(n_old))
    hold=arr
    DEALLOCATE(arr)
    ALLOCATE(arr(n_new))
    
    DO i=1,n_new
       arr(i)=hold(i)
    END DO
    
    DEALLOCATE(hold)

  END SUBROUTINE amputate

  SUBROUTINE fill8(x,n,x1,x2,ilog)

    IMPLICIT NONE
    REAL*8, INTENT(INOUT) :: x(:)
    REAL, INTENT(IN) :: x1, x2
    INTEGER :: i
    INTEGER, INTENT(IN) :: n, ilog

    !Fill an array either linearly or linear-log (ilog=0 or ilog=1)
    !Array of size 'n' should have already been allocated
    
    DO i=1,n
       IF(ilog==0) THEN
          x(i)=x1+(x2-x1)*float(i-1)/float(n-1)
       ELSE IF(ilog==1) THEN
          x(i)=exp(log(x1)+log(x2/x1)*float(i-1)/float(n-1))
       ELSE
          STOP 'FILL8: Error, ilog not properly specified'
       END IF
    END DO

  END SUBROUTINE fill8

  SUBROUTINE ode_crass(x,v,t,ti,tf,xi,vi,fx,fv,n,imeth,ilog)

    IMPLICIT NONE
    REAL, INTENT(IN) :: xi, vi, ti, tf
    REAL :: dt, x4, v4, t4
    REAL :: kx1, kx2, kx3, kx4, kv1, kv2, kv3, kv4
    REAL, ALLOCATABLE :: x(:), v(:), t(:)
    REAL*8, ALLOCATABLE :: x8(:), v8(:), t8(:)
    INTEGER :: i
    REAL, EXTERNAL :: fx, fv
    INTEGER, INTENT(IN) :: imeth, ilog, n

    !Solves 2nd order ODE x''(t) from ti to tf and writes out array of x, v, t values
    !xi and vi are the initial values of x and v (i.e. x(ti), v(ti))
    !fx is what x' is equal to
    !fv is what v' is equal to
    !n is the number of time-steps used in the calculation (user defined, thus crass)
    !imeth selects integration method

    ALLOCATE(x8(n),v8(n),t8(n))
    x8=0.d0
    v8=0.d0
    t8=0.d0
    
    x8(1)=xi
    v8(1)=vi
    CALL fill8(t8,n,ti,tf,ilog)

    DO i=1,n-1

       x4=x8(i)
       v4=v8(i)
       t4=t8(i)

       !Time step change
       dt=t8(i+1)-t8(i)

       IF(imeth==1) THEN

          !Crude method!
          x8(i+1)=x8(i)+fx(x4,v4,t4)*dt
          v8(i+1)=v8(i)+fv(x4,v4,t4)*dt

       ELSE IF(imeth==2) THEN

          !Mid-point method!
          x8(i+1)=x8(i)+fx(x4+fx(x4,v4,t4)*dt/2.,v4+fv(x4,v4,t4)*dt/2.,t4+dt/2.)*dt
          v8(i+1)=v8(i)+fv(x4+fx(x4,v4,t4)*dt/2.,v4+fv(x4,v4,t4)*dt/2.,t4+dt/2.)*dt

       ELSE IF(imeth==3) THEN

          !RK4 (Holy smoke, this is so fast compared to above methods)!
          kx1=dt*fx(x4,v4,t4)
          kv1=dt*fv(x4,v4,t4)
          kx2=dt*fx(x4+kx1/2.,v4+kv1/2.,t4+dt/2.)
          kv2=dt*fv(x4+kx1/2.,v4+kv1/2.,t4+dt/2.)
          kx3=dt*fx(x4+kx2/2.,v4+kv2/2.,t4+dt/2.)
          kv3=dt*fv(x4+kx2/2.,v4+kv2/2.,t4+dt/2.)
          kx4=dt*fx(x4+kx3,v4+kv3,t4+dt)
          kv4=dt*fv(x4+kx3,v4+kv3,t4+dt)

          x8(i+1)=x8(i)+(kx1+(2.*kx2)+(2.*kx3)+kx4)/6.
          v8(i+1)=v8(i)+(kv1+(2.*kv2)+(2.*kv3)+kv4)/6.

       END IF

       !Leave the calculation if the density gets too large (defines 'collapse time' crudely)
       IF(x8(i+1)>1e8) EXIT

    END DO

    ALLOCATE(x(n),v(n),t(n))
    x=x8
    v=v8
    t=t8

  END SUBROUTINE ode_crass

  FUNCTION one(f,v,a)

    !Function used for ODE solver
    IMPLICIT NONE
    REAL :: one
    REAL, INTENT(IN) :: f, v, a
    
    one=1.

  END FUNCTION one

  FUNCTION ff(f,v,a)

    !Function used for ODE solver
    IMPLICIT NONE
    REAL :: ff
    REAL, INTENT(IN) :: f, v, a
    REAL :: f1, f2, f3

    f1=3.*omega_m(a)*(1.+mu(a))/2.
    f2=f**2.
    f3=(1.+AH(a)/H2(a))*f

    ff=(f1-f2-f3)/a

  END FUNCTION ff

  FUNCTION fv(d,v,a)

    !Function used for ODE solver
    IMPLICIT NONE
    REAL :: fv
    REAL, INTENT(IN) :: d, v, a
    REAL :: f1, f2

    f1=3.*omega_c(a)*(1.+mu(a))*d/(2.*(a**2.))
    f2=-(2.+AH(a)/H2(a))*(v/a)

    fv=f1+f2

  END FUNCTION fv

  FUNCTION fvnl(d,v,a)

    !Function used for ODE solver
    IMPLICIT NONE
    REAL :: fvnl
    REAL, INTENT(IN) :: d, v, a
    REAL :: f1, f2, f3

    f1=3.*omega_c(a)*(1.+mu(a))*d*(1.+d)/(2.*(a**2.))
    f2=-(2.+AH(a)/H2(a))*(v/a)
    f3=4.*(v**2.)/(3.*(1.+d))

    fvnl=f1+f2+f3

  END FUNCTION fvnl

  FUNCTION fd(d,v,a)

    !Function used for ODE solver
    IMPLICIT NONE
    REAL :: fd
    REAL, INTENT(IN) :: d, v, a

    fd=v

  END FUNCTION fd

  FUNCTION mu(a)

    !Modified gravity G->(1+mu)*G
    IMPLICIT NONE
    REAL :: mu, a

    IF(img==0) THEN
       mu=0.
    ELSE IF(img==1) THEN
       mu=mua
    ELSE IF(img==2) THEN
       mu=mua*omega_w(a)/om_w
    ELSE IF(img==3) THEN
       mu=1./(3.*beta(a))
    ELSE
       STOP 'Modified gravity model not defined properly'
    END IF

  END FUNCTION mu

  FUNCTION beta(a)

    !DGP beta function
    IMPLICIT NONE
    REAL :: beta, a
    
    !H0rc is a dimensionless parameter

    beta=1.+(2./3.)*sqrt(H2(a))*H0rc*(2.+AH(a)/H2(a))

  END FUNCTION beta

  SUBROUTINE assign_cosmology

    !Routine to assing a cosmological model
    IMPLICIT NONE
    INTEGER :: imod

    om_m=1.0
    om_v=0.0
    om_w=0.0
    om_nu=0.
    w0=-1.
    wm=-1.
    am=0.
    dm=1.
    mua=0.
    wa=0.
    astar=1.
    nblip=1
    iw=0
    img=0

    WRITE(*,*) 'Choose cosmological model'
    WRITE(*,*) ' 1 - EdS'
    WRITE(*,*) ' 2 - Flat LCDM'
    WRITE(*,*) ' 3 - Flat w(a)CDM'
    WRITE(*,*) ' 4 - QUICC - INV1'
    WRITE(*,*) ' 5 - QUICC - INV2'
    WRITE(*,*) ' 6 - QUICC - SUGRA'
    WRITE(*,*) ' 7 - QUICC - 2EXP'
    WRITE(*,*) ' 8 - QUICC - AS'
    WRITE(*,*) ' 9 - QUICC - CNR'
    WRITE(*,*) '10 - Open or closed (Ov=0.)'
    WRITE(*,*) '11 - Non-flat w(a)CDM'
    WRITE(*,*) '12 - MG, constant mu'
    WRITE(*,*) '13 - Simpson mu parametrisation'
    WRITE(*,*) '14 - Blip dark energy'
    WRITE(*,*) '15 - Non-flat wCDM'
    WRITE(*,*) '16 - Flat DGP'
    WRITE(*,*) '17 - EdS DGP (Om_m=1.)'
    WRITE(*,*) '18 - Flat wCDM (constant w)'
    WRITE(*,*) '19 - nuLCDM'
    READ(*,*) imod
    WRITE(*,*)

    IF(imod==1) THEN
       !EdS
       om_m=1.
       om_v=0.
    ELSE IF(imod==2) THEN
       !Flat LCDM
       WRITE(*,*) 'Om_m:'
       READ(*,*) om_m
       WRITE(*,*)
       om_v=1.-om_m
    ELSE IF(imod==3) THEN
       !Flat w(a)CDM
       iw=2
       WRITE(*,*) 'Flat w(a)CDM'
       WRITE(*,*) 'w(a)=w0+wa*(1.-a)'
       WRITE(*,*) 'w0:'
       READ(*,*) w0
       WRITE(*,*) 'wa:'
       READ(*,*) wa
       WRITE(*,*) 'Om_m:'
       READ(*,*) om_m
       WRITE(*,*)
       om_w=1.-om_m
    ELSE IF(imod==4) THEN
       !QUICC - 
       om_m=0.26
       om_w=0.74
       iw=1
       w0=-0.4
       wm=-0.27
       am=0.18
       dm=0.5
    ELSE IF(imod==5) THEN
       !QUICC - 
       om_m=0.26
       om_w=0.74
       iw=1
       w0=-0.79
       wm=-0.67
       am=0.29
       dm=0.4
    ELSE IF(imod==6) THEN
       !QUICC - 
       om_m=0.26
       om_w=0.74
       iw=1
       w0=-0.82
       wm=-0.18
       am=0.1
       dm=0.7
    ELSE IF(imod==7) THEN
       !QUICC - 
       om_m=0.26
       om_w=0.74
       iw=1
       w0=-1.
       wm=0.01
       am=0.19
       dm=0.043
    ELSE IF(imod==8) THEN
       !QUICC - 
       om_m=0.26
       om_w=0.74
       iw=1
       w0=-0.96
       wm=-0.01
       am=0.53
       dm=0.13
    ELSE IF(imod==9) THEN
       !QUICC - 
       om_m=0.26
       om_w=0.74
       iw=1
       w0=-1.
       wm=0.1
       am=0.15
       dm=0.016
    ELSE IF(imod==10) THEN
       !Open or closed matter-dominated model
       WRITE(*,*) 'Om_m:'
       READ(*,*) om_m
       WRITE(*,*)
       om_w=0.
       om_v=0.
    ELSE IF(imod==11) THEN
       !Non-flat w(a)CDM
       WRITE(*,*) 'Non-flat w(a)CDM'
       WRITE(*,*) 'w(a)=w0+wa*(1.-a)'
       WRITE(*,*) 'w0:'
       READ(*,*) w0
       WRITE(*,*) 'wa:'
       READ(*,*) wa
       WRITE(*,*) 'Om_m:'
       READ(*,*) om_m
       WRITE(*,*) 'Om_w:'
       READ(*,*) om_w
       WRITE(*,*)
    ELSE IF(imod==12) THEN
       WRITE(*,*) 'MG: G -> (1.+mu)G, mu is constant'
       WRITE(*,*) 'mu:'
       READ(*,*) mua
       WRITE(*,*) 'Om_m:'
       READ(*,*) om_m
       WRITE(*,*)
       img=1
       om_v=1.-om_m
    ELSE IF(imod==13) THEN
       !Constant G change MG parameterisation
       img=2
       WRITE(*,*) 'MG: G -> (1.+mu)G, mu=mu0*om_v(a)/om_v'
       WRITE(*,*) 'mu0:'
       READ(*,*) mua
       WRITE(*,*)
    ELSE IF(imod==14) THEN
       !Blip dark energy
       iw=4
       WRITE(*,*) 'Om_m:'
       READ(*,*) om_m
       WRITE(*,*) 'Om_v:'
       READ(*,*) om_v
       WRITE(*,*) 'Om_w:'
       READ(*,*) om_w
       WRITE(*,*) 'a*:'
       READ(*,*) astar
       WRITE(*,*) 'n'
       READ(*,*) nblip
       WRITE(*,*)
    ELSE IF(imod==15) THEN
       !Non-flat wCDM model
       iw=3
       WRITE(*,*) 'Non-flat wCDM'
       WRITE(*,*) 'Om_m:'
       READ(*,*) om_m
       WRITE(*,*) 'Om_w:'
       READ(*,*) om_w
       WRITE(*,*) 'w:'
       READ(*,*) w0
       WRITE(*,*)       
    ELSE IF(imod==16) THEN
       !Flat DGP
       img=3
       WRITE(*,*) 'Flat DGP'
       WRITE(*,*) 'Om_m:'
       READ(*,*) om_m
       WRITE(*,*)
       om_v=1.-om_m
       WRITE(*,*) 'H0*rc [dimensionless]:'
       READ(*,*) H0rc
       WRITE(*,*)
    ELSE IF(imod==17) THEN
       img=3
       om_m=1.
       om_w=0.
       WRITE(*,*) 'EdS DGP (with Om_m=1.)'
       WRITE(*,*) 'H0*rc [dimensionless]:'
       READ(*,*) H0rc
       WRITE(*,*)
    ELSE IF(imod==18) THEN
       !Flat wCDM
       iw=3
       WRITE(*,*) 'Flat wCDM'
       WRITE(*,*) 'w:'
       READ(*,*) w0
       WRITE(*,*) 'Om_m:'
       READ(*,*) om_m
       WRITE(*,*)
       om_w=1.-om_m
    ELSE IF(imod==19) THEN
       !Massive neutrino model
       WRITE(*,*) 'nuCDM'
       WRITE(*,*) 'om_m:'
       READ(*,*) om_m
       WRITE(*,*) 'om_nu:'
       READ(*,*) om_nu
       WRITE(*,*)
       om_v=1.-om_m
    END IF

    !Omega cold is everything minus neutrinos
    om_c=om_m-om_nu

    !Write out cosmology
    WRITE(*,*) 'Parameters:'
    WRITE(*,*) 'Om_m:', om_m
    WRITE(*,*) 'Om_w:', om_w
    WRITE(*,*) 'Om_v:', om_v
    WRITE(*,*) 'Om_nu:', om_nu
    WRITE(*,*) 'Om_c:', om_c
    WRITE(*,*)
    IF(iw==0) THEN
       WRITE(*,*) 'Vacuum energy'
       WRITE(*,*) 'w:', -1
    ELSE IF(iw==1) THEN
       WRITE(*,*) 'QUICC prescription'
       WRITE(*,*) 'w0:', w0
       WRITE(*,*) 'wm:', wm
       WRITE(*,*) 'am:', am
       WRITE(*,*) 'dm:', dm
    ELSE IF(iw==2) THEN
       WRITE(*,*) 'w(a) = w0+wa(1.-a)'
       WRITE(*,*) 'w0:', w0
       WRITE(*,*) 'wa:', wa
    ELSE IF(iw==3) THEN
       WRITE(*,*) 'Constant w'
       WRITE(*,*) 'w0:', w0
    ELSE IF(iw==4) THEN
       WRITE(*,*) 'Blip dark energy'
       WRITE(*,*) 'a*:', astar
       WRITE(*,*) 'n:', nblip
    END IF
    WRITE(*,*)

  END SUBROUTINE assign_cosmology

  FUNCTION integrate(a,b,f,acc,ilog)

    !Integrates f(x) between x=a and x=b until desired accuracy is reached!

    IMPLICIT NONE
    INTEGER :: i, j, jmax, ilog
    REAL :: integrate, a, b, acc, dx, lima, limb
    INTEGER :: nint
    REAL :: x, fac
    REAL*8 :: sum1, sum2
    REAL, EXTERNAL :: f

    IF(a==b) THEN

       integrate=0.

    ELSE

       sum1=0.d0
       sum2=0.d0

       jmax=20

       IF(ilog==1) THEN
          lima=log(a)
          limb=log(b)
       ELSE
          lima=a
          limb=b
       END IF

       DO j=1,jmax

          nint=10*(2**j)

          DO i=1,nint

             x=lima+(limb-lima)*float(i-1)/float(nint-1)

             IF(ilog==1) THEN
                x=exp(x)
             END IF

             IF(i==1 .OR. i==nint) THEN
                !multiple of 1 for beginning and end and multiple of 2 for middle points!
                fac=1.d0
             ELSE
                fac=2.d0
             END IF

             IF(ilog==1) THEN
                sum2=sum2+fac*f(x)*x
             ELSE
                sum2=sum2+fac*f(x)
             END IF

          END DO

          dx=(limb-lima)/float(nint-1)
          sum2=sum2*dx/2.

          IF(j .NE. 1 .AND. ABS(-1.+sum2/sum1)<acc) THEN
             integrate=sum2
             EXIT
          ELSE IF(j==jmax) THEN
             WRITE(*,*) 'Integration timed out:', a, b
             STOP
          ELSE
             sum1=sum2
             sum2=0.d0
          END IF

       END DO

    END IF

  END FUNCTION integrate

  FUNCTION H2(a)

    !Calculates the dimensionless squared hubble parameter at 'a'
    IMPLICIT NONE
    REAL :: H2
    REAL, INTENT(IN) :: a
    
    !Ignores contributions from radiation (not accurate at high z)!
    !It is also dimensionless (no H_0 factors)
    
    H2=om_m*(a**(-3.))+om_w*X(a)+om_v+(1.-om_m-om_w-om_v)*(a**(-2.))

  END FUNCTION H2

  FUNCTION omega_m(a)

    !This calculates omega_m(a)
    IMPLICIT NONE
    REAL :: omega_m
    REAL, INTENT(IN) :: a

    omega_m=om_m*(a**(-3.))/H2(a)

  END FUNCTION omega_m

  FUNCTION omega_c(a)

    !This calculates omega_c(a)
    IMPLICIT NONE
    REAL :: omega_c
    REAL, INTENT(IN) :: a

    omega_c=om_c*(a**(-3.))/H2(a)

  END FUNCTION omega_c

  FUNCTION omega_nu(a)

    !This calculates omega_nu(a)
    IMPLICIT NONE
    REAL :: omega_nu
    REAL, INTENT(IN) :: a

    omega_nu=om_nu*(a**(-3.))/H2(a)

  END FUNCTION omega_nu

  FUNCTION omega_w(a)

    !This calculates omega_w(a)
    IMPLICIT NONE
    REAL :: omega_w
    REAL, INTENT(IN) :: a

    !X(a) is the redshift scaling of the dark energy
    omega_w=om_w*X(a)/H2(a)

  END FUNCTION omega_w

  FUNCTION omega_v(a)

    !This calculates omega_v(a)!
    IMPLICIT NONE
    REAL :: omega_v
    REAL, INTENT(IN) :: a

    omega_v=om_v/H2(a)

  END FUNCTION omega_v
  
  FUNCTION AH(a)

    !\ddot a/a
    IMPLICIT NONE
    REAL :: AH
    REAL, INTENT(IN) :: a

    !Gets no contribution from curvature

    AH=om_m*(a**(-3.))+om_w*(1.+3.*w(a))*x(a)-2.*om_v
    AH=-AH/2.

  END FUNCTION AH

  FUNCTION X(a)

    !Redshift scaling for dark energy (i.e., if w=0 X(a)=a^-3, if w=-1 X(a)=const etc.)
    IMPLICIT NONE
    REAL :: X
    REAL, INTENT(IN) :: a

    IF(iw==0) THEN
       !Vacuum
       X=1.   
    ELSE IF(iw==1) THEN
       !Generally true, doing this integration makes the spherical-collapse calculation very slow
       X=(a**(-3.))*exp(3.*integrate(a,1.,integrand,0.001,1))
    ELSE IF(iw==2) THEN
       !w(a)CDM
       X=(a**(-3.*(1.+w0+wa)))*exp(-3.*wa*(1.-a))
    ELSE IF(iw==3) THEN
       !wCDM
       X=a**(-3.*(1.+w0))
    ELSE IF(iw==4) THEN
       !Blip dark energy
       X=(((a/astar)**nblip+1.)/((1./astar)**nblip+1.))**(-6./nblip)
    ELSE
       STOP 'Dark energy model not defined'
    END IF

  END FUNCTION X

  FUNCTION integrand(a)

    !The integrand for one of my integrals
    IMPLICIT NONE
    REAL :: integrand
    REAL, INTENT(IN) :: a

    integrand=w(a)/a

  END FUNCTION integrand

  FUNCTION w(a)

    !Variations of the dark energy parameter w(a)
    IMPLICIT NONE
    REAL :: w
    REAL, INTENT(IN) :: a
    REAL :: p1, p2, p3, p4

    IF(iw==0) THEN
       !Vacuum
       w=-1.
    ELSE IF(iw==1) THEN
       !QUICC parameterisation
       p1=1.+exp(am/dm)
       p2=1.-exp(-(a-1.)/dm)
       p3=1.+exp(-(a-am)/dm)
       p4=1.-exp(1./dm)
       w=w0+(wm-w0)*p1*p2/(p3*p4)
    ELSE IF(iw==2) THEN
       !w(a)CDM
       w=w0+(1.-a)*wa
    ELSE IF(iw==3) THEN
       !wCDM
       w=w0
    ELSE IF(iw==4) THEN
       !Blip dark energy
       w=((a/astar)**nblip-1.)/((a/astar)**nblip+1.)
    ELSE
       STOP 'Dark energy model not defined'
    END IF

  END FUNCTION w

  FUNCTION find(x,xin,yin,n,iorder,imeth)

    IMPLICIT NONE
    REAL :: find
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: x, xin(n), yin(n)
    REAL, ALLOCATABLE ::  xtab(:), ytab(:)
    REAL :: a, b, c, d
    REAL :: x1, x2, x3, x4
    REAL :: y1, y2, y3, y4
    INTEGER :: i
    INTEGER, INTENT(IN) :: imeth, iorder

    !This version interpolates if the value is off either end of the array!
    !Care should be chosen to insert x, xtab, ytab as log if this might give better!
    !Results from the interpolation!

    !If the value required is off the table edge the interpolation is always linear

    !imeth = 1 => find x in xtab by crudely searching from x(1) to x(n)
    !imeth = 2 => find x in xtab quickly assuming the table is linearly spaced
    !imeth = 3 => find x in xtab using midpoint splitting (iterations=CEILING(log2(n)))

    !iorder = 1 => linear interpolation
    !iorder = 2 => quadratic interpolation
    !iorder = 3 => cubic interpolation

    ALLOCATE(xtab(n),ytab(n))

    xtab=xin
    ytab=yin

    IF(xtab(1)>xtab(n)) THEN
       !Reverse the arrays in this case
       CALL reverse(xtab,n)
       CALL reverse(ytab,n)
    END IF

    IF(x<xtab(1)) THEN

       !Do a linear interpolation beyond the table boundary

       x1=xtab(1)
       x2=xtab(2)

       y1=ytab(1)
       y2=ytab(2)

       CALL fit_line(a,b,x1,y1,x2,y2)
       find=a*x+b
       
    ELSE IF(x>xtab(n)) THEN

       !Do a linear interpolation beyond the table boundary
       
       x1=xtab(n-1)
       x2=xtab(n)

       y1=ytab(n-1)
       y2=ytab(n)

       CALL fit_line(a,b,x1,y1,x2,y2)
       find=a*x+b

    ELSE IF(iorder==1) THEN

       IF(n<2) STOP 'FIND: Not enough points in your table for linear interpolation'

       IF(x<=xtab(2)) THEN

          x1=xtab(1)
          x2=xtab(2)

          y1=ytab(1)
          y2=ytab(2)

       ELSE IF (x>=xtab(n-1)) THEN

          x1=xtab(n-1)
          x2=xtab(n)

          y1=ytab(n-1)
          y2=ytab(n)

       ELSE

          i=table_integer(x,xtab,n,imeth)
          
          x1=xtab(i)
          x2=xtab(i+1)

          y1=ytab(i)
          y2=ytab(i+1)

       END IF

       CALL fit_line(a,b,x1,y1,x2,y2)
       find=a*x+b

    ELSE IF(iorder==2) THEN

       IF(n<3) STOP 'FIND: Not enough points in your table'

       IF(x<=xtab(2) .OR. x>=xtab(n-1)) THEN

          IF(x<=xtab(2)) THEN

             x1=xtab(1)
             x2=xtab(2)
             x3=xtab(3)

             y1=ytab(1)
             y2=ytab(2)
             y3=ytab(3)

          ELSE IF (x>=xtab(n-1)) THEN

             x1=xtab(n-2)
             x2=xtab(n-1)
             x3=xtab(n)

             y1=ytab(n-2)
             y2=ytab(n-1)
             y3=ytab(n)

          END IF

          CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)

          find=a*(x**2.)+b*x+c

       ELSE

          i=table_integer(x,xtab,n,imeth)

          x1=xtab(i-1)
          x2=xtab(i)
          x3=xtab(i+1)
          x4=xtab(i+2)

          y1=ytab(i-1)
          y2=ytab(i)
          y3=ytab(i+1)
          y4=ytab(i+2)

          !In this case take the average of two separate quadratic spline values

          find=0.

          CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)
          find=find+(a*(x**2.)+b*x+c)/2.

          CALL fit_quadratic(a,b,c,x2,y2,x3,y3,x4,y4)
          find=find+(a*(x**2.)+b*x+c)/2.

       END IF

    ELSE IF(iorder==3) THEN

       IF(n<4) STOP 'FIND: Not enough points in your table'

       IF(x<=xtab(3)) THEN

          x1=xtab(1)
          x2=xtab(2)
          x3=xtab(3)
          x4=xtab(4)        

          y1=ytab(1)
          y2=ytab(2)
          y3=ytab(3)
          y4=ytab(4)

       ELSE IF (x>=xtab(n-2)) THEN

          x1=xtab(n-3)
          x2=xtab(n-2)
          x3=xtab(n-1)
          x4=xtab(n)

          y1=ytab(n-3)
          y2=ytab(n-2)
          y3=ytab(n-1)
          y4=ytab(n)

       ELSE

          i=table_integer(x,xtab,n,imeth)

          x1=xtab(i-1)
          x2=xtab(i)
          x3=xtab(i+1)
          x4=xtab(i+2)

          y1=ytab(i-1)
          y2=ytab(i)
          y3=ytab(i+1)
          y4=ytab(i+2)

       END IF

       CALL fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)
       find=a*x**3.+b*x**2.+c*x+d

    ELSE

       STOP 'FIND: Error, interpolation order specified incorrectly'

    END IF

  END FUNCTION find

  FUNCTION table_integer(x,xtab,n,imeth)

    IMPLICIT NONE
    INTEGER :: table_integer
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: x, xtab(n)
    INTEGER, INTENT(IN) :: imeth

    IF(imeth==1) THEN
       table_integer=linear_table_integer(x,xtab,n)
    ELSE IF(imeth==2) THEN
       table_integer=search_int(x,xtab,n)
    ELSE IF(imeth==3) THEN
       table_integer=int_split(x,xtab,n)
    ELSE
       STOP 'TABLE INTEGER: Method specified incorrectly'
    END IF

  END FUNCTION table_integer

  FUNCTION linear_table_integer(x,xtab,n)

    IMPLICIT NONE
    INTEGER :: linear_table_integer
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: x, xtab(n)
    REAL :: x1, x2, xn
    REAL :: acc

    !Returns the integer (table position) below the value of x
    !eg. if x(3)=6. and x(4)=7. and x=6.5 this will return 6
    !Assumes table is organised linearly (care for logs)

    x1=xtab(1)
    x2=xtab(2)
    xn=xtab(n)

    !Test for linear table
    acc=0.001

    IF(x1>xn) STOP 'LINEAR_TABLE_INTEGER :: table in the wrong order'
    IF(ABS(-1.+float(n-1)*(x2-x1)/(xn-x1))>acc) STOP 'LINEAR_TABLE_INTEGER :: table does not seem to be linear'

    linear_table_integer=1+FLOOR(float(n-1)*(x-x1)/(xn-x1))

  END FUNCTION linear_table_integer

  FUNCTION search_int(x,xtab,n)

    IMPLICIT NONE
    INTEGER :: search_int
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: x, xtab(n)
    INTEGER :: i

    IF(xtab(1)>xtab(n)) STOP 'SEARCH_INT: table in wrong order'

    DO i=1,n
       IF(x>=xtab(i) .AND. x<=xtab(i+1)) EXIT
    END DO

    search_int=i

  END FUNCTION search_int

  FUNCTION int_split(x,xtab,n)

    IMPLICIT NONE
    INTEGER :: int_split
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: x, xtab(n)
    INTEGER :: i1, i2, imid

    !Finds the position of the value in the table by continually splitting it in half

    IF(xtab(1)>xtab(n)) STOP 'INT_SPLIT: table in wrong order'

    i1=1
    i2=n

    DO
       
       imid=NINT((i1+i2)/2.)

       IF(x<xtab(imid)) THEN
          i2=imid
       ELSE
          i1=imid
       END IF

       IF(i2==i1+1) EXIT

    END DO
    
    int_split=i1

  END FUNCTION int_split

  SUBROUTINE fit_line(a1,a0,x1,y1,x2,y2)

    IMPLICIT NONE
    REAL, INTENT(OUT) :: a0, a1
    REAL, INTENT(IN) :: x1, y1, x2, y2

    !Given xi, yi i=1,2 fits a line between these points

    a1=(y2-y1)/(x2-x1)
    a0=y1-a1*x1

  END SUBROUTINE fit_line

  SUBROUTINE fit_quadratic(a2,a1,a0,x1,y1,x2,y2,x3,y3)

    IMPLICIT NONE
    REAL, INTENT(OUT) :: a0, a1, a2
    REAL, INTENT(IN) :: x1, y1, x2, y2, x3, y3

    !Given xi, yi i=1,2,3 fits a quadratic between these points

    a2=((y2-y1)/(x2-x1)-(y3-y1)/(x3-x1))/(x2-x3)
    a1=(y2-y1)/(x2-x1)-a2*(x2+x1)
    a0=y1-a2*(x1**2.)-a1*x1

  END SUBROUTINE fit_quadratic

  SUBROUTINE fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)

    IMPLICIT NONE
    REAL, INTENT(OUT) :: a, b, c, d
    REAL, INTENT(IN) :: x1, y1, x2, y2, x3, y3, x4, y4
    REAL :: f1, f2, f3

    !Given xi, yi i=1,2,3,4 fits a cubic between these points

    f1=(y4-y1)/((x4-x2)*(x4-x1)*(x4-x3))
    f2=(y3-y1)/((x3-x2)*(x3-x1)*(x4-x3))
    f3=(y2-y1)/((x2-x1)*(x4-x3))*(1./(x4-x2)-1./(x3-x2))

    a=f1-f2-f3

    f1=(y3-y1)/((x3-x2)*(x3-x1))
    f2=(y2-y1)/((x2-x1)*(x3-x2))
    f3=a*(x3+x2+x1)

    b=f1-f2-f3

    f1=(y4-y1)/(x4-x1)
    f2=a*(x4**2.+x4*x1+x1**2.)
    f3=b*(x4+x1)

    c=f1-f2-f3

    d=y1-a*x1**3.-b*x1**2.-c*x1

  END SUBROUTINE fit_cubic

  SUBROUTINE reverse(arry,n)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(INOUT) :: arry(n)
    INTEGER :: i
    REAL, ALLOCATABLE :: hold(:) 

    !This reverses the contents of arry!

    ALLOCATE(hold(n))

    hold=arry

    DO i=1,n
       arry(i)=hold(n-i+1)
    END DO

    DEALLOCATE(hold)

  END SUBROUTINE reverse

  SUBROUTINE ode_adaptive(x,v,t,ti,tf,xi,vi,fx,fv,acc,imeth,ilog)

    !Adaptive ODE solver
    IMPLICIT NONE
    REAL, INTENT(IN) :: xi, vi, ti, tf, acc
    REAL :: dt, x4, v4, t4
    REAL :: kx1, kx2, kx3, kx4, kv1, kv2, kv3, kv4
    REAL*8, ALLOCATABLE :: x8(:), t8(:), v8(:), xh(:), th(:), vh(:)
    REAL, ALLOCATABLE :: x(:), t(:), v(:)
    INTEGER :: i, j, n, k, np, ifail, kn, ninit
    INTEGER, INTENT(IN) :: imeth, ilog

    !Solves 2nd order ODE x''(t) from ti to tf and writes out array of x, v, t values
    !xi and vi are the initial values of x and v (i.e. x(ti), v(ti))
    !fx is what x' is equal to
    !fv is what v' is equal to
    !acc is the desired accuracy across the entire solution
    !imeth selects method

    INTERFACE

       REAL FUNCTION fx(xval,vval,tval)
         REAL, INTENT(IN) :: xval, tval, vval
       END FUNCTION fx
       
       REAL FUNCTION fv(xval,vval,tval)
         REAL, INTENT(IN) :: xval, tval, vval
       END FUNCTION fv
       
    END INTERFACE

    ninit=100

    DO j=1,30

       n=1+ninit*(2**(j-1))

       ALLOCATE(x8(n),v8(n),t8(n))

       x8=0.d0
       v8=0.d0
       t8=0.d0

       x8(1)=xi
       v8(1)=vi
       CALL fill8(t8,n,ti,tf,ilog)

       ifail=0

       DO i=1,n-1

          x4=real(x8(i))
          v4=real(v8(i))
          t4=real(t8(i))

          dt=t8(i+1)-t8(i)

          IF(imeth==1) THEN

             !Crude method!
             x8(i+1)=x8(i)+fx(x4,v4,t4)*dt
             v8(i+1)=v8(i)+fv(x4,v4,t4)*dt

          ELSE IF(imeth==2) THEN

             !Mid-point method!
             kx1=dt*fx(x4,v4,t4)
             kv1=dt*fv(x4,v4,t4)
             kx2=dt*fx(x4+kx1/2.,v4+kv1/2.,t4+dt/2.)
             kv2=dt*fv(x4+kx1/2.,v4+kv1/2.,t4+dt/2.)
             
             x8(i+1)=x8(i)+kx2
             v8(i+1)=v8(i)+kv2

          ELSE IF(imeth==3) THEN

             !RK4 (Holy smoke, this is so fast compared to above methods)!
             kx1=dt*fx(x4,v4,t4)
             kv1=dt*fv(x4,v4,t4)
             kx2=dt*fx(x4+kx1/2.,v4+kv1/2.,t4+dt/2.)
             kv2=dt*fv(x4+kx1/2.,v4+kv1/2.,t4+dt/2.)
             kx3=dt*fx(x4+kx2/2.,v4+kv2/2.,t4+dt/2.)
             kv3=dt*fv(x4+kx2/2.,v4+kv2/2.,t4+dt/2.)
             kx4=dt*fx(x4+kx3,v4+kv3,t4+dt)
             kv4=dt*fv(x4+kx3,v4+kv3,t4+dt)

             x8(i+1)=x8(i)+(kx1+(2.*kx2)+(2.*kx3)+kx4)/6.
             v8(i+1)=v8(i)+(kv1+(2.*kv2)+(2.*kv3)+kv4)/6.

          END IF
          
       END DO

       IF(j==1) ifail=1

       IF(j .NE. 1) THEN

          np=1+(n-1)/2

          DO k=1,1+(n-1)/2

             kn=2*k-1

             IF(ifail==0) THEN

                IF(xh(k)>acc .AND. x8(kn)>acc .AND. (ABS(xh(k)/x8(kn))-1.)>acc) ifail=1
                IF(vh(k)>acc .AND. v8(kn)>acc .AND. (ABS(vh(k)/v8(kn))-1.)>acc) ifail=1

                IF(ifail==1) THEN
                   DEALLOCATE(xh,th,vh)
                   EXIT
                END IF

             END IF
          END DO

       END IF

       IF(ifail==0) THEN
          WRITE(*,*) 'ODE: Integration complete in steps:', n-1
          WRITE(*,*)
          ALLOCATE(x(n),t(n),v(n))
          x=x8
          v=v8
          t=t8
          EXIT
       END IF

       WRITE(*,*) 'ODE: Integration at:', n-1
       ALLOCATE(xh(n),th(n),vh(n))
       xh=x8
       vh=v8
       th=t8
       DEALLOCATE(x8,t8,v8)

    END DO

  END SUBROUTINE ode_adaptive

   SUBROUTINE fill_table(min,max,arr,n,ilog)

    !Fills array 'arr' in equally spaced intervals
    IMPLICIT NONE
    INTEGER :: i
    REAL, INTENT(IN) :: min, max
    REAL :: a, b
    REAL, ALLOCATABLE :: arr(:)
    INTEGER, INTENT(IN) :: ilog, n

    !ilog=0 does linear spacing
    !ilog=1 does log spacing

    IF(ALLOCATED(arr)) DEALLOCATE(arr)

    ALLOCATE(arr(n))

    arr=0.

    IF(ilog==0) THEN
       a=min
       b=max
    ELSE IF(ilog==1) THEN
       a=log(min)
       b=log(max)
    END IF

    IF(n==1) THEN
       arr(1)=a
    ELSE IF(n>1) THEN
       DO i=1,n
          arr(i)=a+(b-a)*float(i-1)/float(n-1)
       END DO
    END IF

    IF(ilog==1) arr=exp(arr)

  END SUBROUTINE fill_table

  FUNCTION inttab_range(x,y,n,n1,n2,iorder)

    !Integrates tables y(x)dx
    IMPLICIT NONE
    REAL :: inttab_range
    INTEGER, INTENT(IN) :: n, n1, n2
    REAL, INTENT(IN) :: x(n), y(n)
    REAL :: a, b, c, d, h
    REAL :: q1, q2, q3, qi, qf
    REAL :: x1, x2, x3, x4, y1, y2, y3, y4, xi, xf
    REAL*8 :: sum
    INTEGER :: i, i1, i2, i3, i4
    INTEGER, INTENT(IN) :: iorder

    sum=0.d0

    IF(n2<n1) STOP 'INTTAB_RANGE: Error n2 must be greater than n1'

    IF(n1==n2) THEN

       sum=0.d0
    
    ELSE IF(iorder==1) THEN

       !Sums over all Trapezia (a+b)*h/2
       DO i=n1,n2-1
          a=y(i+1)
          b=y(i)
          h=x(i+1)-x(i)
          sum=sum+(a+b)*h/2.d0
       END DO

    ELSE IF(iorder==2) THEN

       DO i=n1,n2-2

          x1=x(i)
          x2=x(i+1)
          x3=x(i+2)

          y1=y(i)
          y2=y(i+1)
          y3=y(i+2)

          CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)

          q1=a*(x1**3.)/3.+b*(x1**2.)/2.+c*x1
          q2=a*(x2**3.)/3.+b*(x2**2.)/2.+c*x2
          q3=a*(x3**3.)/3.+b*(x3**2.)/2.+c*x3

          !Takes value for first and last sections but averages over sections where you
          !have two independent estimates of the area
          IF(n==3) THEN
             sum=sum+q3-q1
          ELSE IF(i==1) THEN
             sum=sum+(q2-q1)+(q3-q2)/2.d0
          ELSE IF(i==n-2) THEN
             sum=sum+(q2-q1)/2.d0+(q3-q2)
          ELSE
             sum=sum+(q3-q1)/2.
          END IF

       END DO

    ELSE IF(iorder==3) THEN

       DO i=n1,n2-1

          !First choose the integers used for defining cubics for each section
          !First and last are different because the section does not lie in the *middle* of a cubic

          IF(i==1) THEN

             i1=1
             i2=2
             i3=3
             i4=4

          ELSE IF(i==n-1) THEN

             i1=n-3
             i2=n-2
             i3=n-1
             i4=n

          ELSE

             i1=i-1
             i2=i
             i3=i+1
             i4=i+2

          END IF

          x1=x(i1)
          x2=x(i2)
          x3=x(i3)
          x4=x(i4)

          y1=y(i1)
          y2=y(i2)
          y3=y(i3)
          y4=y(i4)

          CALL fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)

          !These are the limits of the particular section of integral
          xi=x(i)
          xf=x(i+1)

          qi=a*(xi**4.)/4.+b*(xi**3.)/3.+c*(xi**2.)/2.+d*xi
          qf=a*(xf**4.)/4.+b*(xf**3.)/3.+c*(xf**2.)/2.+d*xf

          sum=sum+qf-qi

       END DO

    END IF

    inttab_range=sum

  END FUNCTION inttab_range

END PROGRAM collapse
