!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute mean 2PCF multipoles from Nres=1000 EZmocks
! for any issue please contact ginevra.favole@port.ac.uk
!-------------------------------------------------------
!compile as:
!  ifort -fopenmp mean_multipoles_real.f90 -o meanmreal
!execute as:
! ./meanmreal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main
implicit none
CHARACTER(len=150) :: FN
CHARACTER(len=150) :: GN
CHARACTER(150), PARAMETER   :: output='../multipoles_EZmocks3Gpc/multipoles_5Mpc/mean_multipoles_real.txt'

integer, parameter :: Nres=999,cov_num=40
integer :: i, j, k, io, N_halsum
real*8 :: mxi0(1:cov_num), mxi2(1:cov_num), mxi4(1:cov_num), s(1:cov_num), temp
real*8, dimension(:), allocatable  :: bin_mock, xi0_mock, xi2_mock, xi4_mock


!------------------------------------------------ format input read all
mxi0=0.0
mxi2=0.0
mxi4=0.0

DO I=0,Nres
WRITE(FN,10)I
WRITE(6,*)FN
j=1
OPEN(1,FILE=FN)
DO
READ(1,*,IOSTAT=io) temp
IF (io < 0) THEN
N_halsum = j-1
EXIT
ENDIF
j =j+1
ENDDO

allocate(bin_mock(1:N_halsum),xi0_mock(1:N_halsum),xi2_mock(1:N_halsum),xi4_mock(1:N_halsum))
rewind(1)
do k=1,N_halsum
read(1,*,IOSTAT=io) bin_mock(k),xi0_mock(k),xi2_mock(k),xi4_mock(k)
enddo

do k=1,N_halsum
s(k) = bin_mock(k)
mxi0(k) = mxi0(k)+xi0_mock(k)
mxi2(k) = mxi2(k)+xi2_mock(k)
mxi4(k) = mxi4(k)+xi4_mock(k)
enddo
CLOSE(1)
DEALLOCATE(bin_mock, xi0_mock, xi2_mock, xi4_mock)
ENDDO



open(unit=1020,file=output)
do k=1,N_halsum
mxi0(k) = mxi0(k)/dble(Nres+1.)
mxi2(k) = mxi2(k)/dble(Nres+1.)
mxi4(k) = mxi4(k)/dble(Nres+1.)
write(1020,'(f10.4,x,f20.8,x,f20.8,x,f20.8)') s(k), mxi0(k), mxi2(k), mxi4(k)
enddo
close(1020)


10 FORMAT('../multipoles_EZmocks3Gpc/multipoles_5Mpc/xilreal_',I3.3,'.txt')
end program
