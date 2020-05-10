!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Covariance matrices for redshift-space EZmocks multipoles
! for any issue please contact ginevra.favole@port.ac.uk
!-------------------------------------------------------
!compile as:
!  ifort -fopenmp covariances_rsd.f90 -o covrsd
!execute as:
! ./covrsd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main
  implicit none
  CHARACTER(150), PARAMETER   :: inputcat='multipoles/mean_multipoles_rsd.txt'
  CHARACTER(150), PARAMETER   :: output00='covariances/covariance_xi0xi0_rsd.txt'
  CHARACTER(150), PARAMETER   :: output22='covariances/covariance_xi2xi2_rsd.txt'
  CHARACTER(150), PARAMETER   :: output44='covariances/covariance_xi4xi4_rsd.txt'
  CHARACTER(150), PARAMETER   :: output02='covariances/covariance_xi0xi2_rsd.txt'
  CHARACTER(150), PARAMETER   :: output04='covariances/covariance_xi0xi4_rsd.txt'
  CHARACTER(150), PARAMETER   :: output24='covariances/covariance_xi2xi4_rsd.txt'
  CHARACTER(150), PARAMETER   :: outputsigma00='covariances/sigma_xi0_rsd.txt'
  CHARACTER(150), PARAMETER   :: outputsigma22='covariances/sigma_xi2_rsd.txt'
  CHARACTER(150), PARAMETER   :: outputsigma44='covariances/sigma_xi4_rsd.txt'

  CHARACTER(len=150) :: FN
  CHARACTER(len=150) :: GN
  integer, parameter :: Nres=999,cov_num=40
  real*8, dimension(:), allocatable  :: s, xi0m, xi2m, xi4m, s_mock, xi0_mock, xi2_mock, xi4_mock
  real*8, dimension(cov_num,cov_num) :: cov00,cov22,cov44,cov02,cov04,cov24,b1
  integer :: l, i, j, k, io, media_num, N_halsum
  real*8 :: temp


!----------------------------------------- read input mean catalogue:
  open(unit=1400,file=inputcat)
  media_num = 0
200 read(1400,*,end=201) b1
  media_num = media_num+1
  goto 200
201 close(1400)

  allocate(s(media_num),xi0m(media_num),xi2m(media_num),xi4m(media_num))

  open(unit=1400,file=inputcat)
  do l=1,media_num
     read(1400,*) s(l),xi0m(l),xi2m(l),xi4m(l)
  end do
  close(1400)


  !------------------------------------------------ compute covariances
cov00=0.
cov22=0.
cov44=0.
cov02=0.
cov04=0.
cov24=0.


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
allocate(s_mock(1:N_halsum),xi0_mock(1:N_halsum), xi2_mock(1:N_halsum), xi4_mock(1:N_halsum))
rewind(1)

do k=1,N_halsum
read(1,*,IOSTAT=io) s_mock(k), xi0_mock(k), xi2_mock(k), xi4_mock(k)
enddo


do l=1,N_halsum
do j=1,N_halsum
cov00(l,j)=cov00(l,j)+(xi0_mock(l)-xi0m(l))*(xi0_mock(j)-xi0m(j))
cov22(l,j)=cov22(l,j)+(xi2_mock(l)-xi2m(l))*(xi2_mock(j)-xi2m(j))
cov44(l,j)=cov44(l,j)+(xi4_mock(l)-xi4m(l))*(xi4_mock(j)-xi4m(j))
cov02(l,j)=cov02(l,j)+(xi0_mock(l)-xi0m(l))*(xi2_mock(j)-xi2m(j))
cov04(l,j)=cov04(l,j)+(xi0_mock(l)-xi0m(l))*(xi4_mock(j)-xi4m(j))
cov24(l,j)=cov24(l,j)+(xi2_mock(l)-xi2m(l))*(xi4_mock(j)-xi4m(j))
enddo
enddo

CLOSE(1)
DEALLOCATE(s_mock,xi0_mock, xi2_mock, xi4_mock)
ENDDO



!------------------------------------ outputs in matrix format (40 x 40)
open(unit=1570,file=output00)
open(unit=1571,file=output22)
open(unit=1572,file=output44)
open(unit=1573,file=output02)
open(unit=1574,file=output04)
open(unit=1575,file=output24)

open(unit=1576,file=outputsigma00)
open(unit=1577,file=outputsigma22)
open(unit=1578,file=outputsigma44)

do l=1,N_halsum
do j=1,N_halsum
!we have 1000 mocks with 0<=Nres<=999 and here we need to divide by 999:
cov00(l,j)=(1./dble(Nres))*cov00(l,j)
cov22(l,j)=(1./dble(Nres))*cov22(l,j)
cov44(l,j)=(1./dble(Nres))*cov44(l,j)
cov02(l,j)=(1./dble(Nres))*cov02(l,j)
cov04(l,j)=(1./dble(Nres))*cov04(l,j)
cov24(l,j)=(1./dble(Nres))*cov24(l,j)
enddo
write(1570,'(40(e,x))') (cov00(l,j), j=1,N_halsum)
write(1571,'(40(e,x))') (cov22(l,j), j=1,N_halsum)
write(1572,'(40(e,x))') (cov44(l,j), j=1,N_halsum)
write(1573,'(40(e,x))') (cov02(l,j), j=1,N_halsum)
write(1574,'(40(e,x))') (cov04(l,j), j=1,N_halsum)
write(1575,'(40(e,x))') (cov24(l,j), j=1,N_halsum)

write(1576,'(f20.8)') sqrt(cov00(l,l))
write(1577,'(f20.8)') sqrt(cov22(l,l))
write(1578,'(f20.8)') sqrt(cov44(l,l))
enddo
close(1570)
close(1571)
close(1572)
close(1573)
close(1574)
close(1575)
close(1576)
close(1577)
close(1578)


10 FORMAT('multipoles/xilrsd_',I3.3,'.txt')

deallocate(s, xi0m, xi2m, xi4m)
end program main

