program main
  use, intrinsic :: iso_fortran_env
  use stdlib_kinds, only: dp
  use stdlib_strings, only: to_string
  use arnie, only: print_arnie
  use io, only: &
    init_random_seed,&
    read_indat, &
    set_potential, &
    print_indat, &
    setup_sim, &
    print_xyz, &
    count_lines, &
    init_target, &
    load_target => load_target_fv
  use force, only: &
    varystep,&
    varystep_vv,&
    update_ion,& 
    !update_ion => update_ion_old,&
    calc_ion_props
  use potentials, only: print_pot
  use mpi_f08

  implicit none
  !include 'mpif.h'
  real(dp), parameter :: tol=1.0e-15_dp,&
                         len_fact=1.0_dp/0.529_dp,&
                         mass_fact=1822.89_dp,&
                         e_fact= 27.211_dp
  integer(kind=int32), parameter :: expr = range(1.0_dp),&
                                    dpr  = precision(1.0_dp)

  real(dp), allocatable :: a_pos(:,:), a_vel(:,:), a_acel(:,:),&
                           a_zz(:), a_mass(:), ion_iv(:), ion_xy(:,:),&
                           ion_ke_arr(:), ke_tar_arr(:), tan_phi_arr(:), &
                           r_min_arr(:), tan_psi_arr(:), &
                           ion_xy_arr(:,:), a_pos_og(:,:), xx(:), yy(:), &
                           zz(:), ion_xx(:), ion_yy(:)
  real(dp) :: ion_zz, ion_mass, ion_ke, ff, fwhm_qout, sigma_therm, &
              frozen_par, alpha_max, ion_zi, ion_zf, dx_step, acc, &
              n_exp, cell(3), cell_scaled(3), factor, &
              n_cor, n_sta, n_cap, r_min, r0, dt, dt_max, t, ion_trj_len, ddr,&
              ion_trj_max, lam_a, lam_mu, r_cut, gam_p, gam_c, gam_s, gam_cut, &
              tan_phi, tan_psi, ion_ispeed, vp, ke_tar, ke_ion, tmp
  integer, allocatable :: ion_qout_arr(:) 
  integer :: i, j, k, n, nion, verbose, v_type, natom, logfile, &
             nprint, count, is_xyz, ion_qin, ion_qout, i_ion, iofile, &
             myid, ncpu, stat(MPI_STATUS_SIZE), err, nchunk, istart, istop, ii, icpu, &
             outputfile, ntargetatom, rem
  character(len=:), allocatable :: fname_input, fname_target, prename, ion_elem, v_typename,&
    logfilename, xyzfilename, outfilename
  character(len=*),parameter :: mkdir = 'mkdir -p '
  logical :: exists
  
  ! Standard MPI initialisation 
  type(MPI_Datatype) :: dptype
  call MPI_INIT(err)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,err)
  call MPI_Type_Create_F90_real(dpr, expr, dptype, err)
  
  !how many procs do we have? 
  if (myid == 0 .and. verbose > 0) print*, ncpu, ' processors for calculation'
  do icpu=0,ncpu-1
    if ((verbose > 0) .and. (myid == icpu)) print*, 'hello from proc', myid
    call MPI_BARRIER( MPI_COMM_WORLD, err)
  end do  

  ! hard code for now
  fname_input = 'indat.in'

  ! call with 1 for same seed 
  call init_random_seed(0, myid)

  do icpu =0, ncpu 
    if (myid == icpu) then
      call read_indat(fname_input, prename, ion_elem, ion_zz, ion_mass, ion_ke, &
        ion_qin, factor, gam_p, gam_c, gam_s, gam_cut, fwhm_qout, sigma_therm, &
        frozen_par, alpha_max, ion_zi, ion_zf, dx_step, acc, nion, verbose, &
        is_xyz, v_typename, fname_target)
    end if 
    call MPI_BARRIER( MPI_COMM_WORLD, err)
  end do 

  ! some constants
  lam_a = 1.0_dp
  lam_mu = 1.0_dp
  r_cut = 10.0_dp  !*len_fact
  
  ! set potential type switch
  v_type = set_potential(v_typename)
  if ((verbose >0 ) .and. (myid == 0)) then 
    print*, 'Hallo, dis ist Arnold... your instructor.'
    print*, 'VERBOSE MODE ACTIVATED NYGARAHRARHAARHR!'
    call system( mkdir // 'logs')
    !call print_arnie()
  end if 
  
  if (verbose == 2) then
    nprint = 1
  end if
  
  ! check if a file of initial ion coordinates exists  
  ! if yes, then we do as many ions as we have coordinates 
  if (myid == 0) then
    inquire(file='ion.xy', exist=exists)
    if (exists) then
      nion = 0
      if (verbose>0 .and. myid == 0) print*, "initial ion coordinates exists!"
      call count_lines('ion.xy', nion)
      print*,'proc ',myid,'has nion', nion
      allocate(ion_xx(nion))
      allocate(ion_yy(nion))
      allocate(ion_xy(nion,2))

      open(unit=iofile, file='ion.xy', status='old') 
      do i=1, nion 
        read(iofile,*) ion_xx(i), ion_yy(i) 
        ion_xy(i,1) = ion_xx(i)
        ion_xy(i,2) = ion_yy(i)
      end do 
      close(iofile)
    ! otherwise, use number of ions provided
    else 
      allocate(ion_xx(nion))
      allocate(ion_yy(nion))
      allocate(ion_xy(nion,2))
      ion_xx = -1.0_dp
      ion_yy = -1.0_dp
      ion_xy = -1.0_dp
    end if 
    ! okay proc 0 has the initial ion coordinates, tell the other procs
    do icpu =1, ncpu -1
      call MPI_SEND(nion, 1, MPI_INT, icpu, 1, MPI_COMM_WORLD, err)
      call MPI_SEND(ion_xx, size(ion_xx), dptype, icpu, 2, MPI_COMM_WORLD, err)
      call MPI_SEND(ion_yy, size(ion_yy), dptype, icpu, 3, MPI_COMM_WORLD, err)
      print*, 'sent nion',nion, ' to proc',icpu
    end do
  else 
    call MPI_RECV(nion, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
    allocate(ion_xx(nion))
    allocate(ion_yy(nion))
    allocate(ion_xy(nion,2))
    call MPI_RECV(ion_xx, size(ion_xy), dptype, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
    call MPI_RECV(ion_yy, size(ion_xy), dptype, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
  end if 

  do i = 1, nion
    ion_xy(i,1) = ion_xx(i)
    ion_xy(i,2) = ion_yy(i)
  end do 
    
  ! now we partition the work between all the procs
  nchunk = int(nion/ncpu)
  rem = int(mod(nion,ncpu))
  istart = 0
  istop = 0
  do i = 0, ncpu-1
    if (myid == i) then 
      istart = i*nchunk + 1
      istop = (i+1)*nchunk
    end if 
  end do 
  if ((myid == ncpu-1) .and. (rem>0)) then
    istop = istop +  rem
    nchunk = nchunk + rem
  end if 

  ! round call, who has how many ion coordinates
  if (verbose>0) then 
    do i =0, ncpu-1
      if (icpu == myid) then 
        print*, 'proc ',myid, ' istart:',istart,' istop', istop, ' nchunk', nchunk, ' rem',rem
      end if 
      call MPI_BARRIER( MPI_COMM_WORLD, err)
    end do 
  end if 

  ! each proc now needs to allocate for the ion runs
  allocate(ion_ke_arr(nchunk))
  allocate(ke_tar_arr(nchunk))
  allocate(tan_phi_arr(nchunk))
  allocate(ion_qout_arr(nchunk))
  allocate(r_min_arr(nchunk))
  allocate(tan_psi_arr(nchunk))
  allocate(ion_xy_arr(nchunk,2))

  ion_xy_arr = ion_xy(istart:istop,:)
  do icpu=0,ncpu-1
    if (myid == icpu) then
      print*, 'proc:', myid, ' has ion_xy:'
      do i=1,nchunk
        print*, ion_xy_arr(i,:)
      end do 
    end if 
    call MPI_BARRIER( MPI_COMM_WORLD, err)
  end do 
  
  ! setup the arrays for the target atoms
  if (myid == 0) then
    call init_target(fname_target, natom) 
    ! Now natom is the total target atoms + the ion 
    natom = natom + 1
    do icpu=1, ncpu-1
      call MPI_SEND(natom, 1, MPI_int, icpu, 1, MPI_COMM_WORLD, err)
    end do
  else 
    call MPI_RECV(natom, 1, MPI_int, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
  end if 

  allocate(xx(natom))
  allocate(yy(natom))
  allocate(zz(natom))
  allocate(a_pos(natom,3))
  allocate(a_mass(natom))
  allocate(a_zz(natom))
  allocate(a_vel(natom,3))
  allocate(a_acel(natom,3))

  ! read in atom coordinates and distribute them among the procs
  if (myid == 0) then
    call load_target(fname_target, xx, yy, zz, a_mass, a_zz, natom, cell, cell_scaled)
    do icpu= 1, ncpu-1
      call MPI_SEND(a_mass, natom, dptype, icpu, 2, MPI_COMM_WORLD, err)
      call MPI_SEND(a_zz, natom, dptype, icpu, 3, MPI_COMM_WORLD, err)
      call MPI_SEND(cell, 3, dptype, icpu, 4, MPI_COMM_WORLD, err)
      call MPI_SEND(cell_scaled, 3, dptype, icpu, 5, MPI_COMM_WORLD, err)

      call MPI_SEND(xx, natom, dptype, icpu, 6, MPI_COMM_WORLD, err)
      call MPI_SEND(yy, natom, dptype, icpu, 7, MPI_COMM_WORLD, err)
      call MPI_SEND(zz, natom, dptype, icpu, 8, MPI_COMM_WORLD, err)
    end do 
  else 
    call MPI_RECV(a_mass, natom, dptype, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
    call MPI_RECV(a_zz, natom, dptype, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
    call MPI_RECV(cell, 3, dptype, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
    call MPI_RECV(cell_scaled, 3, dptype, 0, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
    
    call MPI_RECV(xx, natom, dptype, 0, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
    call MPI_RECV(yy, natom, dptype, 0, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
    call MPI_RECV(zz, natom, dptype, 0, 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
  end if 
  
  ! create an output file for the results  
  outfilename = prename // '_out.txt'
  if (myid == 0) then
    open(newunit=outputfile, file=outfilename, action='write')
      write(outputfile, *) '#ion_id ion_x ion_y ion_KE_i ion_KE_f tar_KE qout r_min tan_phi tan_psi'
    close(outputfile)
  end if 

  !===============================!
  ! ion fly event starts here ! 
  !===============================!
  ii = 0
  do i_ion = istart, istop 
    ii = ii + 1
    
    a_pos(:,1) = xx
    a_pos(:,2) = yy
    a_pos(:,3) = zz
    a_vel = 0.0_dp
    a_acel = 0.0_dp

    ! logging stuff
    if (verbose > 0) print*, 'Shoot ion ', to_string(i_ion), ' of ', nion
    if (verbose > 1) then
      logfilename = 'logs/' // prename // '_' // to_string(i_ion) // '.log'
      open(newunit=logfile, file=logfilename, action='write')
    end if 
      
    xyzfilename = 'logs/ion_' // to_string(i_ion) // '.xyz'


    call setup_sim(ion_zi, ion_zf, ion_xy_arr(ii,1), ion_xy_arr(ii,2), ion_zz, ion_mass, ion_qin, ion_ke, a_pos, &
                   a_mass, a_zz, a_vel, a_acel, cell, cell_scaled, n_cor, n_sta, n_cap, &
                   factor, ff, r0, r_min, vp, sigma_therm)
    
    if (is_xyz==1) call print_xyz(a_pos, a_mass, a_zz, cell, 0.0_dp, xyzfilename, 'new')

    ! set some things up
    ddr = 0.01_dp
    t = a_pos(1,3)/vp
    ion_trj_max = abs(ion_zf - ion_zi) *len_fact ! stupid fucking factor
    ion_trj_len = abs(a_pos(1,3) - ion_zi * len_fact) ! fucking factor
    dt_max = abs(dx_step/vp)
    tan_phi = 0.0_dp
    tan_psi = 0.0_dP
    ion_qout = ion_qin

    if (verbose > 1) then
      write(logfile,*)  '#step time dt ion_x ion_y ion_z N_core N_stable N_cap ion_disp' 
    end if 

    count = 0
    do while (ion_trj_len < ion_trj_max)
      !RK method
      call varystep(t, a_pos, a_vel, a_acel, a_mass, a_zz, cell_scaled, vp, acc, &
        dt_max, v_type, ff, ddr, n_cap, n_sta, n_cor, r_cut, r0, dt, verbose)

      !VV integrator
      !call varystep_vv(t, a_pos, a_vel, a_acel, a_mass, a_zz, cell_scaled, vp, acc, &
      !  dt_max, v_type, ff, ddr, n_cap, n_sta, n_cor, r_cut, r0, dt, verbose)

      call update_ion(dt, t, a_pos, a_zz, n_sta, n_cap, n_cor, factor, ff, &
        r0, ion_ispeed, lam_a, frozen_par, lam_mu, alpha_max, r_min, gam_p, &
        gam_c, gam_s, gam_cut, cell_scaled, logfile, verbose)
      
      ! Current ion z-displacement w.r.t. initial position
      ion_trj_len = abs(a_pos(1,3) - ion_zi * len_fact)
      t = t + dt

      if ((mod(count,nprint) == 0) .and. (is_xyz == 1)) call print_xyz(a_pos, a_mass, a_zz, cell, t, xyzfilename, 'old')
      if (verbose >1) write(logfile,*)  count, t, dt, a_pos(1,1), a_pos(1,2), a_pos(1,3), n_cor, n_sta, n_cap, ion_trj_len 
      
      count = count + 1
    end do 

    ke_ion = 0.5_dp * a_mass(1) * sum(a_vel(1,:)**2)
    ke_tar = 0.0_dp
    do i=2,natom
      ke_tar = ke_tar + 0.5_dp * a_mass(i) * sum(a_vel(i,:)**2)
    end do

    ! calc ion properties at end of trj
    call calc_ion_props(fwhm_qout, a_vel, a_zz, ion_iv, ion_qin, n_cor, n_sta, n_cap, &
      tan_phi, tan_psi, ion_qout)

    ! Output good bits in a.u. for comparison with pascal code
    ion_ke_arr(ii) = (ion_ke*1000.0/e_fact-ke_ion)*e_fact
    ke_tar_arr(ii) = ke_tar*e_fact
    tan_phi_arr(ii) = tan_phi
    ion_qout_arr(ii) = ion_qout
    r_min_arr(ii) = r_min
    tan_psi_arr(ii) = tan_psi

    if (verbose > 0) then
      write(6, '(i5, 5(f15.5), i4, 2(f15.5))') i_ion, ion_xy_arr(ii,1)*len_fact, ion_xy_arr(ii,2)*len_fact, &
        (ion_ke*1000.0/e_fact-ke_ion)*e_fact, ke_tar*e_fact, tan_phi, ion_qout, r_min, tan_psi
    end if

    if (verbose > 1) close(logfile)
  end do

  do icpu=0, ncpu-1
    if (myid == icpu) then 
      open(newunit=outputfile, file=outfilename, position="append", status='old', action='write')
      do ii=1, nchunk
        i_ion = istart + ii-1
        write(outputfile, '(i6, 4(f15.5), i4, 3(f15.5))') i_ion, ion_xy_arr(ii,1), ion_xy_arr(ii,2), ion_ke_arr(ii), &
          ke_tar_arr(ii), ion_qout_arr(ii), r_min_arr(ii), tan_phi_arr(ii), tan_psi_arr(ii) 
      end do 
      close(outputfile)
    end if 
    call MPI_BARRIER( MPI_COMM_WORLD, err)
  end do 
  call MPI_FINALIZE(err)
end program 
