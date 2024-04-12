!##########################################################
!The following test provides a benchmark for 
!deuteration chemistry in starless cores. 
!It is based on the model of Walmsely et al. 2004.
!##########################################################
program evolve
    use krome_main
    use krome_user
    use krome_user_commons
    use omp_lib
    implicit none
  
    integer,parameter::nx=krome_nmols    !number of species in the krome array (i.e. the real number of species)
    integer::Npart,Nnew
    integer::i,j,irhol,irhoh,k,l,m
    integer::jmin,jmax
    real*8::x(nx),rnum,r
    real*8::dt,spy
    real*8::d2g,opr,Htotal
    real*8::Ms,pc,rhoi,Tgas_i
    real*8::Tdust_i, Av_i !new variables added to work
    real*4::tst,tend,frac
    real*4::t0,t1,rand
    real*8::xint(nx),rhol,rhoh
  
    real*4,allocatable::rho(:),xkrome(:,:),mass(:),xx(:),yy(:),zz(:),Tgas(:),Tdust(:),Av(:)
    real*4::rho_bin(100)
    integer,allocatable::list(:,:)
    integer::nbin(100)
    integer::kmin
    
    real*4::avg(nx),rmin
    integer::ncore
    character(len=10)::fstr
    real*8::pfrac
    real*8::bkgfrac
    call getarg(1,fstr)
    !read(fstr,*) pfrac
  
    open(UNIT=60, FILE='./BinaryFiles/rho.bin', STATUS='OLD', ACTION='READ',FORM='unformatted')
    open(UNIT=61, FILE='./BinaryFiles/xkrome.bin', STATUS='UNKNOWN',FORM='unformatted')
    open(UNIT=62, FILE='./BinaryFiles/mass.bin', STATUS='OLD', ACTION='READ',FORM='unformatted')
    open(UNIT=63, FILE='./BinaryFiles/pos.bin', STATUS='OLD', ACTION='READ',FORM='unformatted')
    open(UNIT=64, FILE='./BinaryFiles/tgas.bin',STATUS='OLD',ACTION='READ',FORM='unformatted')
    open(UNIT=65, FILE='./BinaryFiles/tdust.bin',STATUS='OLD',ACTION='READ',FORM='unformatted')
    open(UNIT=66, FILE='./BinaryFiles/av.bin',STATUS='OLD',ACTION='READ',FORM='unformatted')
  
    spy   = 3600. * 24. * 365. !seconds per year
    pc    = 3.086e18 !pc in cm
    Ms    = 1.989e33 !Solar mass in g
  
    read(60) Npart,t0
    read(62) 
    read(63) 
    read(64)
    read(65)
    read(66)
    print *,'Read initial values',Npart,t0
    
    allocate(rho(Npart))
    allocate(Tgas(Npart))
    allocate(Tdust(Npart))
    allocate(Av(Npart))
    allocate(mass(Npart))
    allocate(xkrome(Npart,nx))
    allocate(list(100,int(Npart*pfrac)))!/10)))
    allocate(xx(Npart))
    allocate(yy(Npart))
    allocate(zz(Npart))
    print*,'Done allocation'
  
    call srand(42)
    j=0
    read(60) rho
    print *,'rho read'
    read(62) mass
    print *,'mass read'
    read(63) xx
    read(63) yy
    read(63) zz
    print *,'pos read'
    read(64) Tgas
    print *,'Tgas read'
    read(65) Tdust
    print *,'Tdust read'
    read(66) Av
    print *,'Av read'

    close(62) 

    call krome_init()
  
    x(:)             = 1.0000000000e-20 !All
    x(KROME_idx_H2)  = 7.1136264572e-01 !H2
    x(KROME_idx_HE)  = 2.8445727019e-01 !He
    x(KROME_idx_H3j) = 3.8407250948e-04 !H3+
    x(KROME_idx_O)   = 1.5474475498e-03 !O
    x(KROME_idx_CO)  = 2.3894410696e-03 !CO
    x(KROME_idx_e)   = 1.2313256141e-12 !e-

  !$omp parallel do
    do i=1,Npart
      xkrome(i,:)=x(:)
    enddo
  !$omp end parallel do
   
    print *,'krome_nmols',krome_nmols
  
    write(61) krome_nmols
    write(61) Npart,t0
    !Start
    do 
     !print *,'xkrome',xkrome
     write(61) Npart,t0
     write(61) xkrome
     read(60,END=9999) Nnew,t1
     read(63,END=9999)
     read(64,END=9999)
     read(65,END=9999)
     read(66,END=9999)
    
     dt   = t1-t0
     print *,'begin from t ',t0,' for ',dt
     call cpu_time(tst)

     call krome_set_user_Tdust(1d1)
     call krome_set_user_crate(1.0d-17)

  !$omp parallel do private(x,rhoi,Tgas_i)
     do i = 1, Npart
       x          = xkrome(i,:)
       rhoi       = rho(i)
       Tgas_i     = Tgas(i)
       Tdust_i    = Tdust(i)
       Av_i       = Av(i)

       call krome_set_user_Av(Av_i)
       call krome_set_user_Tdust(Tdust_i)
       call krome_init_coef(Tgas_i)
       call krome(x(:), rhoi, Tgas_i, dt) !call KROME
  
       xkrome(i,:)=x
     enddo 
  !$omp end parallel do
     print*,'finish parallel'
     call cpu_time(tend)
  
     !Update density
     read(60) rho
     print *,'update rho'
     read(63) xx
     read(63) yy
     read(63) zz
     print *,'update pos'  
     read(64) Tgas
     print *,'update Tgas'  
     read(65) Tdust
     print *,'update Tdust'
     read(66) Av
     print *,'update Av'
  
     t0=t1
    enddo
  
  9999  close(60)
    close(61)
    close(63)
    close(64)
    close(65)
    close(66)
    deallocate(rho,xkrome,xx,yy,zz,Tgas,Tdust,list,Av)
    print*,'Done evolve file!'
  
  end program evolve
  