program converter
use minitools
use initsolver
implicit none

integer, parameter                    :: numchains = 1
integer, parameter                    :: numcolSN  = 12
integer, parameter                    :: numcolSNBAO  = 16
integer                               :: nlines, c2ind, c3ind,c4ind,OMind,H0ind
real*8, dimension(:,:), allocatable   :: incolSN
real*8, dimension(:,:), allocatable   :: incolSNBAO
real*8                                :: fake, c2dh, c3dh, c4dh, OM, H0, betadh
character(100)                        :: chnum
    real*8            :: diff               !difference between H0 and computed H(z=0)
    real*8            :: time1, time2       !variables for time computation
    real*8            :: finalhubble        !H(z=0) after minizer
    !BRENT variables
    real*8            :: aa,amin, aastep
    real*8            :: bb,bmin
    real*8            :: value
    integer             :: status
    integer             :: stepmin
    real*8            :: arg

    real*8 :: condrealmat, condrealds, condreal   !beta>condreal gives complex initial conditions
    logical, parameter  :: minimizeme = .true.
    integer  :: iter
    integer, parameter :: maxiter = 1
    real*8, parameter :: minitol = 0.01

integer :: i,j,k

!Doing JLA
do i=1,numchains
   write(chnum,*) i
   
   nlines = 0 
   OPEN (666, file = '../chains/JLA_DHOST_'//trim(adjustl(chnum))//'.txt') 
   DO 
       READ (666,*, END=10) fake
       nlines = nlines + 1 
   END DO 
   10 CLOSE (666)

   allocate(incolSN(nlines,numcolSN+1))

   open(42,file='../chains/JLA_DHOST_'//trim(adjustl(chnum))//'.txt')
   do j=1,nlines
      read(42,*) (incolSN(j,k),k=1,numcolSN)
   end do
   close(42)

   !READ DHOST AND OMEGA AND H
   c2ind   = 5
   c3ind   = 6
   c4ind   = 7
   OMind   = 3
   H0ind   = 4
write(*,*) 'MINIMIZER STARTS!!'
do j=1,nlines
   !RUN MINIMIZER TO GET BETA FOR EACH CHAIN POINT
   c2dh = incolSN(j,c2ind)
   c3dh = incolSN(j,c3ind)
   c4dh = incolSN(j,c4ind)
   OM = incolSN(j,OMind)
   H0 = incolSN(j,H0ind)


   condrealmat = (c3dh**2./c2dh)-(48./9.)*c4dh !Real xi in matter era
   condrealds = (c3dh**2./(2.*c2dh))-(8./3.)*c4dh !Real xi in de Sitter

   condreal = min(condrealmat,condrealds)

   bb = condreal  !upper limit of the minimizing interval
   if (condreal.gt.0.d0) then
      aa = -5*bb
   else
      aa = 2.*bb
   end if
   stepmin = 0

   arg = aa!bb-iter*aastep
   call deinterface(OM,H0,c2dh,c3dh,c4dh,arg,diff)
   value =diff
   arg = bb
   call deinterface(OM,H0,c2dh,c3dh,c4dh,arg,diff)
   value =diff

   amin = aa!bb-iter*aastep
   bmin = bb
   status = 0

   do

      call local_min_rc ( amin, bmin, arg, status, value )

      if ( status < 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST_LOCAL_MIN_RC_ONE - Fatal error!'
        write ( *, '(a)' ) '  LOCAL_MIN_RC returned negative status.'
        exit
      end if

      call deinterface(OM,H0,c2dh,c3dh,c4dh,arg,diff)
      value =diff
      stepmin = stepmin + 1
      if ( status == 0 ) then
         exit
      end if

   end do
      
   betadh = arg
   call deinterface(OM,H0,c2dh,c3dh,c4dh,betadh,diff)
   call getH(1.d0,finalhubble)

   incolSN(j,numcolSN+1) = betadh
write(*,*) j, betadh
       !if (iter.eq.maxiter) then
   if (diff.gt.minitol) then
      write(*,*) 'THIS SHOULD NOT HAPPEN'
      stop
   end if
end do
   

   open(42,file='../converted_chains/JLA_DHOST_converted_'//trim(adjustl(chnum))//'.txt')
   do j=1,nlines
      write(42,"(200(1X,E15.8))") (incolSN(j,k),k=1,numcolSN+1)
   end do
   close(42)

   deallocate(incolSN)

end do

!Doing JLA+BAO
do i=1,numchains
   write(chnum,*) i
   
   nlines = 0 
   OPEN (666, file = '../chains/JLA+BAO_DHOST_'//trim(adjustl(chnum))//'.txt') 
   DO 
       READ (666,*, END=11) fake
       nlines = nlines + 1 
   END DO 
   11 CLOSE (666)

   allocate(incolSNBAO(nlines,numcolSNBAO+1))

   open(42,file='../chains/JLA+BAO_DHOST_'//trim(adjustl(chnum))//'.txt')
   do j=1,nlines
      read(42,*) (incolSNBAO(j,k),k=1,numcolSNBAO)
   end do
   close(42)

   !READ DHOST AND OMEGA AND H
   c2ind   = 5
   c3ind   = 6
   c4ind   = 7
   OMind   = 3
   H0ind   = 4
write(*,*) 'MINIMIZER STARTS!!'
do j=1,nlines
   !RUN MINIMIZER TO GET BETA FOR EACH CHAIN POINT
   c2dh = incolSNBAO(j,c2ind)
   c3dh = incolSNBAO(j,c3ind)
   c4dh = incolSNBAO(j,c4ind)
   OM = incolSNBAO(j,OMind)
   H0 = incolSNBAO(j,H0ind)


   condrealmat = (c3dh**2./c2dh)-(48./9.)*c4dh !Real xi in matter era
   condrealds = (c3dh**2./(2.*c2dh))-(8./3.)*c4dh !Real xi in de Sitter

   condreal = min(condrealmat,condrealds)

   bb = condreal  !upper limit of the minimizing interval
   if (condreal.gt.0.d0) then
      aa = -5*bb
   else
      aa = 2.*bb
   end if
   stepmin = 0

   arg = aa!bb-iter*aastep
   call deinterface(OM,H0,c2dh,c3dh,c4dh,arg,diff)
   value =diff
   arg = bb
   call deinterface(OM,H0,c2dh,c3dh,c4dh,arg,diff)
   value =diff

   amin = aa!bb-iter*aastep
   bmin = bb
   status = 0

   do

      call local_min_rc ( amin, bmin, arg, status, value )

      if ( status < 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST_LOCAL_MIN_RC_ONE - Fatal error!'
        write ( *, '(a)' ) '  LOCAL_MIN_RC returned negative status.'
        exit
      end if

      call deinterface(OM,H0,c2dh,c3dh,c4dh,arg,diff)
      value =diff
      stepmin = stepmin + 1
      if ( status == 0 ) then
         exit
      end if

   end do
      
   betadh = arg
   call deinterface(OM,H0,c2dh,c3dh,c4dh,betadh,diff)
   call getH(1.d0,finalhubble)

   incolSNBAO(j,numcolSNBAO+1) = betadh
write(*,*) j, betadh
       !if (iter.eq.maxiter) then
   if (diff.gt.minitol) then
      write(*,*) 'THIS SHOULD NOT HAPPEN'
      stop
   end if
end do
   

   open(42,file='../converted_chains/JLA_DHOST_converted_'//trim(adjustl(chnum))//'.txt')
   do j=1,nlines
      write(42,"(200(1X,E15.8))") (incolSNBAO(j,k),k=1,numcolSNBAO+1)
   end do
   close(42)

   deallocate(incolSNBAO)

end do



end program converter
