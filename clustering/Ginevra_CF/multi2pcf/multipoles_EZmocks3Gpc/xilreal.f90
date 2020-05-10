!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!2PCF real-space multipoles of the EZmocks
! for any issue please contact ginevra.favole@port.ac.uk
!-------------------------------------------------------
!compile as:
!  ifort -fopenmp xilreal.f90 -o xilreal
!execute as:
! ./xilrsd < params_real.inp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Program  main
  implicit none
  CHARACTER(150), PARAMETER   :: params_inp='params_real/params_real.inp'
  CHARACTER*250  CATALOG,DD,XI2D,XIL 
  !--------------------------------------------------------------------------------------------
  real*8 :: Lbox, redshift, omegam, omegal, scale_factor
  real*8, parameter  :: pi=3.14159265,hubble0=100.0, s_min=0.,s_max=200., mu_min=0., mu_max=1.
  integer, parameter :: Nth=32, Nheader=2, Ns_bins=200, Nmu_bins=120

  !---------------------------------------------------------------------------------------------
  integer   :: Ngalbin
  real      :: Ndd
  REAL*8, dimension(1:30e6) :: x, y, z
  real*8   :: s_bin_l(1:Ns_bins),s_bin_r(1:Ns_bins),s_bin_c(1:Ns_bins), &
              mu_bin_l(1:Nmu_bins),mu_bin_r(1:Nmu_bins),mu_bin_c(1:Nmu_bins),mu_bin_width(1:Nmu_bins), &
              dvol(1:Ns_bins,1:Nmu_bins)

  real*8   :: bucket_data_data(1:Ns_bins,1:Nmu_bins,1:Nth), bucket_random(1:Ns_bins,1:Nmu_bins), &
       bucket_data_data_new(1:Ns_bins,1:Nmu_bins,1:Nth),bucket_data_data_new_final(1:Ns_bins,1:Nmu_bins), &
       cf(1:Ns_bins,1:Nmu_bins),mono(1:Ns_bins),quad(1:Ns_bins),hexa(1:Ns_bins)
  real*8   :: x_temp,y_temp,z_temp,junk1

  real*8   :: step_mu,step_s,mu,s,dx,dy,dz,rz
  integer  :: i,j,k,n,io,l,beginning,rate,end

  integer  :: NTHREADS, TID, OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM,OMP_GET_NUM_PROCS
  integer, dimension(3) :: time

  write(*,'(/,/,/)')
  call system_clock(beginning, rate)        
      



  !----------------------- read positions and velocities from galaxy catalog:
   open(11,file=params_inp)
   read (11,'(A)') CATALOG
   read (11,'(A)') DD
   read (11,'(A)') XI2D
   read (11,'(A)') XIL
   read (11,*) Lbox
   read (11,*) redshift
   read (11,*) omegam
   read (11,*) omegal
   read (11,*) scale_factor
   close(11)
   
   write(*,*) 'read input'
   
hubbleEV=dble(hubble0*sqrt(omegam*(1.+redshift)**3.+omegal))

  OPEN(1,FILE=CATALOG,STATUS='OLD',action='read')
  i =1
  DO j=1,Nheader
  READ(1,*)
  ENDDO
  DO WHILE (io >= 0)
     READ(1,*,IOSTAT=io) x_temp,y_temp,z_temp
     x(i) = x_temp
     y(i) = y_temp
     z(i) = z_temp
     i = i+1
  ENDDO
  Ngalbin = i-2
  Ndd = dble(Ngalbin)
  write(*,*) 'Ndd=', Ngalbin
  close(1)


  !--------------------------------- set z-space separation bins log:
  step_s=dble((s_max-s_min)/Ns_bins)
  DO n =1, Ns_bins
     s_bin_r(n) = s_min+(dble(n)*step_s)
     s_bin_l(n) = s_min+(dble(n-1)*step_s)
     s_bin_c(n) = 0.5*(s_bin_l(n) + s_bin_r(n))
  ENDDO

  
  
  !--------------------------------- set mu separation bins linear (Hong+14):

  step_mu=dble((mu_max-mu_min)/Nmu_bins)
  DO k =1, Nmu_bins
        mu_bin_l(k)=mu_min+(dble(k-1)*step_mu)
        mu_bin_r(k)=mu_min+(dble(k)*step_mu)
        mu_bin_c(k) = 0.5*(mu_bin_l(k) + mu_bin_r(k))
        mu_bin_width(k) = dble(mu_bin_r(k) - mu_bin_l(k))
  ENDDO
  

  

  
  !--------------------------------- count DD pairs:

  bucket_data_data = 0
  !$omp parallel shared(x,y,z,dx,dy,dz, &
  !$omp&     bucket_data_data,s,mu,s_bin_l,s_bin_r,mu_bin_l,mu_bin_r) private(i,j,k,l,n,TID) num_threads(Nth)

  TID = OMP_GET_THREAD_NUM() !working thread 
  NTHREADS = OMP_GET_NUM_PROCS() !total available threads 
  write(*,*) 'TID', TID, 'NTHREADS',NTHREADS
  l=TID+1


  !$OMP DO
  DO i = 1, Ngalbin
     if (mod(i,100000)==0) write(*,*) 'step =', i
     DO j = i+1, Ngalbin
        dx = ABS(x(i) - x(j))
        dy = ABS(y(i) - y(j))
        dz = ABS(z(i) - z(j)) 
       
        !!periodic boundary conditions
        if (dx > (Lbox-dx)) then
           dx = Lbox - dx
        endif
        if (dy > (Lbox-dy)) then
           dy = Lbox - dy
        endif
        if (dz > (Lbox-dz)) then
           dz = Lbox - dz
        endif
       
        s=sqrt(dx**2+dy**2+dz**2)
        mu=dble(abs(dz)/s)
        if (mu < mu_bin_r(Nmu_bins)) then
           if (s >= s_bin_l(1) .and. s < s_bin_r(Ns_bins)) then
              n = 1 + INT((mu - mu_min) / step_mu)
              if (n >= 1 .and. n <= Nmu_bins) then
                 k = 1 + INT((s-s_min)/step_s)
                 if (k >= 1 .and. k <= Ns_bins) then
                    bucket_data_data(k,n,l) = bucket_data_data(k,n,l)+ 1.0
                 endif
              endif
           endif
        endif


    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL

write(*,*) "end parallel"


  open(3,file=DD)
  do k=1,Ns_bins
     do n=1,Nmu_bins
           do l=1,Nth
              bucket_data_data_new(k,n,l) = dble(2.0*bucket_data_data(k,n,l)/(Ndd)**2)
           enddo
           bucket_data_data_new_final(k,n) = sum(bucket_data_data_new(k,n,1:Nth))
           write(3,'(f12.5,x,f12.5,x,e)') s_bin_c(k),mu_bin_c(n),bucket_data_data_new_final(k,n)
     enddo
  enddo
  close(3)
  write(*,*) 'DD pairs done'


  !--------------------------------  random pairs using partial volumes shells:

  bucket_random =0
  do k=1,Ns_bins
     do n=1,Nmu_bins
         dvol(k,n) = (4./3.)*pi*(s_bin_r(k)**3 - s_bin_l(k)**3)*mu_bin_width(n) 
         bucket_random(k,n) =  (dvol(k,n)/Lbox**3)
     enddo
  enddo



  !------------------------ calculate 3D correlation function:
  open(71,file=XI2D)
  do k=1,Ns_bins
     do n=1,Nmu_bins
        cf(k,n) = dble(bucket_data_data_new_final(k,n)/bucket_random(k,n) - 1.0)
        write(71,'(f12.5,x,f12.5,x,f12.5)') s_bin_c(k),mu_bin_c(n),cf(k,n)                                                                              
     enddo
  enddo
close(71)

call system_clock(end)
print *, "elapsed time: ", real(end - beginning) / real(rate)

print *, 'computed xi(mu,s)'



!beginning=0
!rate=0
!end=0

!call system_clock(beginning, rate)

  open(72,file=XIL)
  mono=0.0
  quad=0.0
  hexa=0.0
  do k=1,Ns_bins
     do n=1,Nmu_bins
        mono(k) = mono(k)+cf(k,n)*mu_bin_width(n)
        quad(k) = quad(k)+5.*cf(k,n)*0.5*(3.*mu_bin_c(n)**2-1.)*mu_bin_width(n) 
        hexa(k) = hexa(k) +9.*cf(k,n)*(1./8.)*(35.*mu_bin_c(n)**4-30.*mu_bin_c(n)**2+3.)*mu_bin_width(n) 
     enddo
     write(72,'(f10.4,x,f20.8,x,f20.8,x,f20.8)') s_bin_c(k),mono(k),quad(k),hexa(k)
  enddo
 close(72)


!call system_clock(end)
!print *, "elapsed time: ", real(end - beginning) / real(rate)

end program main
