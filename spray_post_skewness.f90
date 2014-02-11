subroutine dump_spray_parameters(io)
! NABEEL QAZI (NQ)- NOV 28, 2013
! routine post processes results for spray/solid particles
! Calculates:-
!	 1. average kinetic energy
! 	 2. mean velocities
! 	 3. maximum velocities
! 	 4. Relative velocity of particles
! 	 4. Particle dispersion

!testing with GitHub
use runtime_m, only: time,run_title
use reference_m
use topology_m
use spray_m!, only: ARRAY_SIZE,fill

implicit none
 integer, intent(in) :: io
 integer, parameter :: iopt = 826
 integer, parameter :: avgopt = 827
 integer, parameter :: maxopt = 828
 integer, parameter :: relopt = 832
 integer, parameter :: dispopt = 835

 integer n,fillsum,pz,tr_counter,num

! real, dimension(ARRAY_SIZE,3) :: vel	!spray velocity
 real spray_ke,u_sum,v_sum,w_sum,u_max,v_max,w_max		!for each processor
 real spray_ke_avg,u_avg,v_avg,w_avg,spray_ke_sum,u_max_t,v_max_t,w_max_t	!for p0. average ke at a time
 real u_rel_sum,v_rel_sum,w_rel_sum
 real u_rel_m_sq,v_rel_m_sq,w_rel_m_sq
 real x_displacement_sum,y_displacement_sum,z_displacement_sum
 real x_dispersion,y_dispersion,z_dispersion
 character time_ext*10,myid_ext*5

! call compute_physical_velocity(vel,1)	!get velocities
! call advance_spray(1)	!get velocities

! write(time_ext,'(1pe10.4)') time*time_ref
! write(myid_ext,'(i5.5)') xz_id

 u_sum = 0.0
 v_sum = 0.0
 u_rel_sum = 0.0
 v_rel_sum = 0.0
 w_rel_sum = 0.0
 x_displacement_sum = 0.0
 y_displacement_sum = 0.0
 z_displacement_sum = 0.0
 
 spray_ke = 0.0
 spray_ke_sum = 0.0 
 tr_counter = 0
 
 if (myid .eq. 0) then		!initialise stuff on master node
	spray_ke_avg = 0.0; num = 0
	u_avg = 0.0;
	v_avg = 0.0;
	u_rel_m_sq = 0.0;
	v_rel_m_sq = 0.0;
	w_rel_m_sq = 0.0;
	x_dispersion = 0.0
	y_dispersion = 0.0
	z_dispersion = 0.0
 endif

 do n = 1,fill
!  if (count(state(1:fill) .eq. HEALTHY) .eq. 0) cycle
   if (state(n) .ne. HEALTHY) cycle
  spray_ke = spray_ke+0.5*((vel_d(n,1))**2 + (vel_d(n,2))**2 + (vel_d(n,3))**2)	!0.5*mass*vel^2 for spray/solid particles
  u_sum = u_sum + vel_d(n,1)
  v_sum = v_sum + vel_d(n,2)
  w_sum = w_sum + vel_d(n,3)
  u_rel_sum = u_rel_sum + (vel_d(n,1) - vel_g(n,1))**2
  v_rel_sum = v_rel_sum + (vel_d(n,2) - vel_g(n,2))**2
  w_rel_sum = w_rel_sum + (vel_d(n,3) - vel_g(n,3))**2
  x_displacement_sum = x_displacement_sum + (loc(n,1) - loc_init(n,1))**2
  y_displacement_sum = y_displacement_sum + (loc(n,2) - loc_init(n,2))**2
  z_displacement_sum = z_displacement_sum + (loc(n,3) - loc_init(n,3))**2

	tr_counter = tr_counter+1
 end do
 
 u_max = maxval(vel_d(:,1))
 v_max = maxval(vel_d(:,2))
 w_max = maxval(vel_d(:,3))

!debug check NQ 09/02/2014
! u_max = maxval(vel_g(:,1))
! v_max = maxval(vel_g(:,2))
! w_max = maxval(vel_g(:,3))

 call mpi_barrier(gcomm,ierr)

 !send parameters to p0
 call MPI_allreduce(fill,fillsum,1,MPI_INTEGER,MPI_SUM,gcomm,ierr)
 call MPI_allreduce(tr_counter,num,1,MPI_INTEGER,MPI_SUM,gcomm,ierr)

 call MPI_allreduce(u_sum,u_avg,1,MPI_REAL8,MPI_SUM,gcomm,ierr)
 call MPI_allreduce(v_sum,v_avg,1,MPI_REAL8,MPI_SUM,gcomm,ierr)
 call MPI_allreduce(w_sum,w_avg,1,MPI_REAL8,MPI_SUM,gcomm,ierr)

 call MPI_allreduce(u_max,u_max_t,1,MPI_REAL8,MPI_MAX,gcomm,ierr)
 call MPI_allreduce(v_max,v_max_t,1,MPI_REAL8,MPI_MAX,gcomm,ierr)
 call MPI_allreduce(w_max,w_max_t,1,MPI_REAL8,MPI_MAX,gcomm,ierr)

 call MPI_allreduce(u_rel_sum,u_rel_m_sq,1,MPI_REAL8,MPI_SUM,gcomm,ierr)
 call MPI_allreduce(v_rel_sum,v_rel_m_sq,1,MPI_REAL8,MPI_SUM,gcomm,ierr)
 call MPI_allreduce(w_rel_sum,w_rel_m_sq,1,MPI_REAL8,MPI_SUM,gcomm,ierr)

 call MPI_allreduce(x_displacement_sum,x_dispersion,1,MPI_REAL8,MPI_SUM,gcomm,ierr)
 call MPI_allreduce(y_displacement_sum,y_dispersion,1,MPI_REAL8,MPI_SUM,gcomm,ierr)
 call MPI_allreduce(z_displacement_sum,z_dispersion,1,MPI_REAL8,MPI_SUM,gcomm,ierr)

 call MPI_allreduce(spray_ke,spray_ke_sum,1,MPI_REAL8,MPI_SUM,gcomm,ierr)
! spray_ke_avg = spray_ke_sum/fillsum

  spray_ke_avg = spray_ke_sum/num
  u_avg = u_avg/num
  v_avg = v_avg/num
  w_avg = w_avg/num
  u_rel_m_sq = u_rel_m_sq/num
  v_rel_m_sq = v_rel_m_sq/num
  w_rel_m_sq = w_rel_m_sq/num
  x_dispersion = x_dispersion/num
  y_dispersion = y_dispersion/num
  z_dispersion = z_dispersion/num

 if (myid .eq. 0) then
	open(unit=iopt, file = '../post/spray_ke.dat',access='append')
!        open(unit=avgopt, file = '../post/spray_vel_avg.dat',access='append')
	open(unit=maxopt, file = '../post/spray_vel_max.dat',access='append')
	open(unit=relopt, file = '../post/spray_rel_vel.dat',access='append')
	open(unit=dispopt, file = '../post/spray_dispersion.dat',access='append')
!	    write(iopt,*) time,' ',fillsum
!	    write(iopt,*) time,' ',spray_ke_avg,' ',fillsum
!	    write(iopt,*) time,' ',spray_ke_avg,' ',num
	    write(iopt,'(2(1pe12.5,1x),2(1i5.5))') time,spray_ke_avg,num,fillsum
!	    write(avgopt,'(4(1pe12.5,1x))') time,u_avg,v_avg,w_avg
	    write(maxopt,'(4(1pe12.5,1x))') time,u_max_t,v_max_t,w_max_t
	    write(relopt,'(4(1pe12.5,1x))') time,u_rel_m_sq,v_rel_m_sq,w_rel_m_sq
	    write(dispopt,'(4(1pe12.5,1x))') time,u_rel_m_sq,v_rel_m_sq,w_rel_m_sq
	close(iopt)
!	close(avgopt)
	close(maxopt)	 
	close(relopt)
	close(dispopt)
 end if	!p0 writes the file

! LOOP_PZ: do pz = 0, npes-1
!  call mpi_barrier(gcomm,ierr)
!  if (pz .ne. myid) cycle
!  if (count(state(1:fill) .eq. HEALTHY) .eq. 0) cycle
!   open(unit=iopt, file = '../post/tecplot/fill_spray-'//time_ext//'.dat',access='append')
!    write(iopt,*) pz,' - ',fill
!   close(iopt)
! end do LOOP_PZ

! call mpi_barrier(gcomm,ierr)
 !Find total ! number of sprays
! call MPI_allreduce(fill,fillsum,1,MPI_INTEGER,MPI_SUM,gcomm,ierr) 

 ! For debug
! if (myid .eq. 0) then
! 	open(unit=iopt, file = '../post/tecplot/fill_spray-'//time_ext//'.dat',access='append')
! 	 write(iopt,*) fillsum
! 	close(iopt)
! end if
return	
end subroutine dump_spray_parameters



