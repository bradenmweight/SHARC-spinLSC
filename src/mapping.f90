!...
!...This file contains all the subroutines relatated mapping variables
!...
!...initialize mapping variables
!...nstates=# of states
subroutine init_map(nstates, initState) ! This is true initial state, not focused foward backward. Only use for density matrix.
use definitions
use misc
implicit none
integer,intent(in)::nstates, initState
integer::i,j
real*8,allocatable::rands(:)
real*8,allocatable::r(:),phi(:)
REAL, PARAMETER :: PI_CONST = 3.1415927

if( initState==0 ) stop "spin-PLDM initial states are wrong. MCH is wrong! This is not Python"
write(*,*)"Zhou,nstates,initState",nstates,initState

allocate(Dmatrix(nstates,nstates))  !density matrix
Dmatrix=0.d0

! Allocate Mapping variables
!QD basis
allocate(z(nstates))
z=0.d0

!diabatic basis
allocate(zd(nstates))
zd=0.d0

allocate(z0(nstates))
z0=0.d0

! Allocate Radii
allocate(r(nstates))
r=0.d0

! Allocate Angles
allocate(phi(nstates))
phi=0.d0

!allocate random variables
allocate(rands(nstates)) ! Only the angles are random for spin-LSC

do i=1,nstates
   call random_number(rands(i))
   write(*,*) "Randoms:", rands(i)
enddo

! ZPE parameter, which is global
gw = (2/nstates) * ( sqrt( real(nstates) + 1.0) - 1 )

! Calculate radii for each oscillator
do i=1,nstates
   r(i) = sqrt( gw )
enddo

! Radii for initial forward and backward states are slightly bigger
r(initState) = sqrt( 2 + gw )

! Calculate initial mapping variables
do i=1,nstates
   phi(i) = rands(i) * 2 * PI_CONST
   z(i) = cmplx( r(i) * cos( phi(i) ) , r(i) * sin( phi(i) ) ) ! Complex number
enddo

! We will need to keep track of initial state for density matrix
! zF0,zB0 are global
do i=1,nstates
   z0(i) = z(i)
   write(*,*)"Z(i), i:",z(i)
enddo

! Check initial density matrix. We need true initial state now
do i=1,nstates
   do j=1,nstates
      Dmatrix(i,j) = 0.5d0 * ( conjg(z(i)) * z(j) - gw  )
   enddo
enddo

write(*,*) "Initial Density Matrix:", Dmatrix(:,:)


deallocate(rands)
return
end subroutine init_map
!
!...propagate mapping varialbes using Verlet method
!
subroutine propagate_mapping(traj,ctrl)
   use definitions
   use matrix
   implicit none
   type(trajectory_type) :: traj
   type(ctrl_type) :: ctrl
   real*8,allocatable::Hnow(:,:),Hold(:,:),Hqd(:,:),Sss(:,:),Hdiff(:,:)
   real*8,allocatable::zreal(:),zimag(:)
   real*8::dt,Hmax,dtmax
   integer::istate,jstate,i,n,stepmin
   
   write(*,*)"Verlet"
   !in sharc the Hamiltonian matrix is complex
   !we change it to real.
   
   !write(*,*)"Zhou,nstates=",ctrl%nstates
   allocate(Hnow(ctrl%nstates,ctrl%nstates))  !current Hamiltonian in adiabatic basis t2
   allocate(Hold(ctrl%nstates,ctrl%nstates))  !Old Hamiltonian in QD t1
   allocate(Hqd(ctrl%nstates,ctrl%nstates))   !QD Hamiltonian
   allocate(Sss(ctrl%nstates,ctrl%nstates))   !overlap matrix
   allocate(Hdiff(ctrl%nstates,ctrl%nstates)) !=HQD(t2)-HQD(t1)

   allocate( zreal(ctrl%nstates) )
   allocate( zimag(ctrl%nstates) )

  
   Hnow=real(traj%H_MCH_ss)
   Hold=real(traj%H_MCH_old_ss)
   Sss=real(traj%overlaps_ss)  !change complex to real
   
   !adjust the nsubsteps
   Hmax=-1.d0
   do istate=1,ctrl%nstates
     if(abs(Hnow(istate,istate))>Hmax)then
        Hmax=abs(Hnow(istate,istate))
     endif
   enddo 
   dtmax=2.d0*pi/Hmax/200.d0  !period/200
   stepmin=nint(ctrl%dtstep/dtmax)
   write(*,*)"stepmin=",stepmin

   call transform(ctrl%nstates,Hnow,Sss,'uaut')

   Hdiff=Hnow-Hold  !H(t2)-H(t1) at QD basis at t1

   dt=ctrl%dtstep/ctrl%nsubsteps !in a.u. 
   
   zreal = real(z)
   zimag = imag(z)

   do i=1,ctrl%nsubsteps
      
      !interpolate Hamiltonian at t in [t1,t2]
      Hnow=Hold+Hdiff*dble(i)/dble(ctrl%nsubsteps)     ! This removes interpolation

      zimag = zimag - 0.5 * matmul( Hnow, zreal ) * dt
      zreal = zreal + 1.0 * matmul( Hnow, zimag ) * dt
      zimag = zimag - 0.5 * matmul( Hnow, zreal ) * dt
   
   enddo
   
   z = cmplx( zreal, zimag )
   
   deallocate(Hnow,Hold,Hqd,Sss,Hdiff)
   return
   end subroutine propagate_mapping

!
!...transform mapping variables into t2 QD basis
!
subroutine transform_mapping(traj,ctrl)
use definitions
use matrix
implicit none
type(trajectory_type) :: traj
type(ctrl_type) :: ctrl
real*8,allocatable::Sss(:,:)
real*8,allocatable::zreal(:),zimag(:)
real*8,allocatable::zrealtmp(:),zimagtmp(:)
integer::istate,jstate

allocate(Sss(ctrl%nstates,ctrl%nstates))   !overlap matrix
Sss=real(traj%overlaps_ss)

allocate( zreal(ctrl%nstates) )
allocate( zimag(ctrl%nstates) )

allocate( zrealtmp(ctrl%nstates) )
allocate( zimagtmp(ctrl%nstates) )

zrealtmp = 0.d0
zimagtmp = 0.d0


zreal = real(z)
zimag = imag(z)

zreal = matmul( transpose(Sss), zreal )
zimag = matmul( transpose(Sss), zimag )

z = cmplx( zreal, zimag )

deallocate(Sss); deallocate(zreal,zimag)
return
end subroutine transform_mapping

