!
!...This code can do nonadiabatic dyanmics simulation of nulear and eletron on
!...equal footing using quasi-diabatic (QD) propagation method based on
!...SHARC. So we borrowed many subroutines from SHARC.
!...The algorithm and the details of the method are given in
!...        A. Mandal et al. JCTC 2018, 14, 1828-1840.
!...We need Hamiltonian, grdients and overlap Sik=<Phi_i(t0)|Phi_k(t1)>
!...which can be obatined using SHARC.
!...Wanghuai Zhou,11/1/2018
!
program main
use definitions
use electronic
use input
use matrix
use misc
use nuclear
use qm
use restart
use output
implicit none

!> \param traj Contains all data which would be private to each trajectory in an
!ensemble
type(trajectory_type) :: traj
!> \param ctrl Contains all data which would be shared in an ensemble
type(ctrl_type) :: ctrl
!> \param i_step Loop variable for the dynamics loop
integer :: i_step
!> \param time Define the integer function time()
integer :: time
integer::istate,jstate

traj%time_start=time()
traj%time_last=traj%time_start

!read input file
call read_input(traj,ctrl)
call allocate_lapack(ctrl%nstates)

!initialize mapping variables
call init_map(ctrl%nstates,traj%state_mch,ctrl%init_cond)
 
!call get_Dmatrix(ctrl)  ! BMW, we skip this one because we did it in init_map.

!Do initial QM calculation to get 
!Hamiltonian matrix elements: H_MCH_ss
!Gradients: grad_MCH_sad
!Nonadiabatic coupling vectors: NACdr_ssad
!Overlap matrix: overlaps_ss
!Construct Gmatrix: Gmatrix_ssqd--diagonal gradients + non-adiabatic coupling

call do_first_qm(traj,ctrl)  !quantum chemistry calculation


!Calculate the initial Force

call calc_force(traj,ctrl)

!Save old variables
call save_old(traj)

call Calculate_etot(traj,ctrl) !traj%Etot=traj%Ekin+traj%Epot
call set_time(traj)            !
call write_dat(u_dat, traj, ctrl) 
call write_list_line(u_lis,traj,ctrl)
call write_geom(u_geo, traj, ctrl)
call write_Dmatrix(traj,ctrl) !output some matrix properties

traj%overlaps_ss = traj%H_MCH_ss * 0.d0
do istate=1,ctrl%nstates
    traj%overlaps_ss(istate,istate) = cmplx(1.0,0.0)
enddo


!evolution loop
do i_step=traj%step+1,ctrl%nsteps
!do i_step=traj%step+1,1
   traj%step=i_step
   write(*,*)"Zhou,i_step=",i_step
   call write_logtimestep(u_log,i_step)

  !...Calculate Force of previous position with half-step mapping variables
   call calc_force(traj,ctrl)

   !...Verlet step for nuclear positions
   call VelocityVerlet_xstep(traj,ctrl) 

   !...Do QM calculation to get
   !...H_MCH_ss and grad_MCH_sad
   !...In sharc, grad_MCH_sad is transformed into Gmatrix_ss and then Gmatrix_ssad
   !call do_qm_calculations(traj,ctrl)
   !Construct Gmatrix: diagonal gradients + non-adiabatic coupling
   call do_qm_calc_qd(traj,ctrl)

   call calc_solve_polariton_Hamiltonian(traj,ctrl)

  !...Verlet for nuclear velocity
  !call VelocityVerlet_vstep(traj,ctrl)
   call VelocityVerlet_vstep_Half(traj,ctrl)

   !...begin substep for mapping variable
   !...Calcualate the HQD(t) throught interpelation
   !...Propagate the forward and backward mapping variables by solving Hamilton's
   !...equation with Hamiltonian hmF and hmB

   !...Calculate Force of new position with half-step mapping with new position
   call calc_force(traj,ctrl)

   !...Verlet for nuclear velocity
   !call VelocityVerlet_vstep(traj,ctrl)
   call VelocityVerlet_vstep_Half(traj,ctrl)

   if(ctrl%integrator==1)then
    call propagate_mapping(traj,ctrl)
  endif

  !...Transforming mapping variables into t2 QD basis
  call transform_mapping(traj,ctrl)

   !...Save variables
   call save_old(traj)
   call set_time(traj)
   call write_list_line(u_lis,traj,ctrl)
   call write_dat(u_dat, traj, ctrl)
   call write_geom(u_geo, traj, ctrl)
   call allflushqd()
  
   !...calculate the density matrix (adiabatic population)
   call get_Dmatrix(ctrl,traj)
   call write_Dmatrix(traj,ctrl)

enddo !i_step loop

call write_final(traj)

end program main
