program main
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

  implicit none
  include 'mpif.h'
  real(dp), parameter :: tol=1.0e-15_dp,&
                         len_fact=1.0_dp/0.529_dp,&
                         mass_fact=1822.89_dp,&
                         e_fact= 27.211_dp

  real(dp), allocatable :: a_pos(:,:), a_vel(:,:), a_acel(:,:),&
                           a_zz(:), a_mass(:), ion_iv(:), ion_xy(:,:),&
                           ion_ke_arr(:), ke_tar_arr(:), tan_phi_arr(:), &
                           ion_qout_arr(:), r_min_arr(:), tan_psi_arr(:), &
                           ion_xy_arr(:,:)
  real(dp) :: ion_zz, ion_mass, ion_ke, ff, fwhm_qout, sigma_therm, &
              frozen_par, alpha_max, ion_zi, ion_zf, dx_step, acc, surface_cov, &
              n_exp, exp_arg_limit, dx_exp, cell(3), cell_scaled(3), factor, &
              n_cor, n_sta, n_cap, r_min, r0, dt, dt_max, t, ion_trj_len, ddr,&
              ion_trj_max, lam_a, lam_mu, r_cut, gam_p, gam_c, gam_s, gam_cut, &
              tan_phi, tan_psi, ion_ispeed, vp, ke_tar, ke_ion

  integer :: i, j, k, n, nion, verbose, v_type, natom, logfile, &
             nprint, count, is_xyz, ion_qin, ion_qout, i_ion, iofile, &
             myid, ncpu, stat(MPI_STATUS_SIZE), err, nchunk, istart, istop, ii, icpu, &
             outputfile, ntargetatom
  character(len=:), allocatable :: fname_input, fname_target, prename, ion_elem, v_typename,&
    logfilename, xyzfilename, outfilename
  character(len=*),parameter :: mkdir = 'mkdir -p '
  logical :: exists
  
  ! Standard MPI initialisation 
  call MPI_INIT(err)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,err)

  ! hard code for now
  fname_input = 'indat.in'

  ! call with 1 for same seed 
  call init_random_seed(0, myid)

  !========================!
  ! read in the input file !
  !========================!
  call read_indat(fname_input, prename, ion_elem, ion_zz, ion_mass, ion_ke, &
      ion_qin, factor, gam_p, gam_c, gam_s, gam_cut, fwhm_qout, sigma_therm, &
      frozen_par, alpha_max, ion_zi, ion_zf, dx_step, acc, nion, verbose, &
      surface_cov, v_typename, fname_target)
  ! all input values are unscaled !

  ! some constants
  lam_a = 1.0_dp
  lam_mu = 1.0_dp
  r_cut = 10.0_dp  !*len_fact
  
  
  ! set potential type switch
  v_type = set_potential(v_typename)
  if ((verbose >0 ) .and. (myid == 0)) then 
    print*, 'Hallo, dis ist Arnold... your instructor.'
    print*, 'VERBOSE MODE ACTIVATED NYGARAHRARHAARHR!'
    !call print_arnie()
  end if 

  call system( mkdir // 'logs')
  
  if (verbose == 2) then
    is_xyz = 1
    nprint = 100
  end if
  
  ! check if a file of initial ion coordinates exists  
  ! if yes, then we do as many ions as we have coordinates 
  inquire(file='ion.xy', exist=exists)
  if (exists) then
    nion = 0
    if (verbose>0) print*, "initial ion coordinates exists!"
    nion = count_lines('ion.xy')
    allocate(ion_xy(nion,2))

    open(unit=iofile, file='ion.xy', status='old') 
    do i=1, nion 
      read(iofile,*) ion_xy(i,1), ion_xy(i,2) 
    end do 
    close(iofile)

  ! otherwise, use number of ions provided
  else 
    allocate(ion_xy(nion,2))
    ion_xy = -1.0_dp
  end if 

  !how many procs do we have? 
  if (myid == 0 .and. verbose > 0) print*, ncpu, ' processors for calculation'
  if (verbose > 0) print*, 'hello from proc',myid
  
  nchunk = nion/ncpu
  istart = 0
  istop = 0
  do i = 0, ncpu-1
    if (myid == i) then 
      istart = i*nchunk + 1
      istop = (i+1)*nchunk
    end if 
    if ((myid == ncpu-1) .and. (mod(1.0*nion, 1.0*ncpu)>0)) then
      istop = istop +  int(mod(1.0*nion, 1.0*ncpu))
      nchunk = nchunk + int(mod(1.0*nion, 1.0*ncpu))
    end if 
    if (myid == i) print*, 'proc ',myid, ' istart:',istart,' istop', istop, ' nchunk', nchunk

  end do 

  allocate(ion_ke_arr(nchunk))
  allocate(ke_tar_arr(nchunk))
  allocate(tan_phi_arr(nchunk))
  allocate(ion_qout_arr(nchunk))
  allocate(r_min_arr(nchunk))
  allocate(tan_psi_arr(nchunk))
  allocate(ion_xy_arr(nchunk,2))

  !===============================!
  ! one ion fly event starts here ! 
  !===============================!
  call init_target(fname_target, a_pos, a_mass, a_zz, a_vel, a_acel, natom) 
  ! Now natom is the total target atoms + the ion 
  natom = natom + 1

  ii = 0
  do i_ion = istart, istop 
    ii = ii + 1
    ! logging stuff
    if (verbose > 0) print*, 'Shoot ion ', to_string(i_ion), ' of ', nion
    xyzfilename = 'logs/ion_' // to_string(i_ion) // '.xyz'

    logfilename = 'logs/' // prename // '_' // to_string(i_ion) // '.log'
    open(newunit=logfile, file=logfilename, action='write')

    ! Load target atoms and allocate arrays for target atoms, sets z==0 at the surface
    ! Note that all atom lengths and masses are now scaled, except the input values
    call load_target(fname_target, sigma_therm, a_pos, a_vel, a_acel, a_mass, &
      a_zz, natom, cell, cell_scaled)
    
    ! here, natoms is the number of target atoms + the ion
    if (verbose > 0) print*, natom-1, ' target atoms in sample'

    call setup_sim(ion_zi, ion_zf, ion_xy(i_ion,:), ion_zz, ion_mass, ion_qin, ion_ke, a_pos, &
                   a_mass, a_zz, a_vel, a_acel, cell, cell_scaled, n_cor, n_sta, n_cap, &
                   factor, ff, r0, r_min, vp)
    
    ion_xy_arr(ii,:) = a_pos(1,:2)/len_fact
  
    if (is_xyz==1) call print_xyz(a_pos, a_mass, a_zz, cell, xyzfilename, 'new')

    ! set some things up
    ddr = 0.01_dp
    t = a_pos(1,3)/vp
    ion_trj_max = abs(ion_zf - ion_zi) *len_fact ! stupid fucking factor
    ion_trj_len = abs(a_pos(1,3) - ion_zi * len_fact) ! fucking factor
    dt_max = abs(dx_step/vp)
    tan_phi = 0.0_dp
    tan_psi = 0.0_dP
    ion_qout = ion_qin

    if (verbose >= 2) then
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

      if ((mod(count,nprint) == 0) .and. (is_xyz == 1)) call print_xyz(a_pos, a_mass, a_zz, cell, xyzfilename, 'old')
      if (verbose >=2) write(logfile,*)  count, t, dt, a_pos(1,1), a_pos(1,2), a_pos(1,3), n_cor, n_sta, n_cap, ion_trj_len 
      
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
      print*, i_ion, ion_xy(i_ion,1)*len_fact, ion_xy(i_ion,2)*len_fact, (ion_ke*1000.0/e_fact-ke_ion)*e_fact, &
        ke_tar*e_fact, tan_phi, ion_qout, r_min, tan_psi
    end if 
    close(logfile)
  end do

  outfilename = prename // '_out.txt'
  close(outputfile)
  do icpu=0, ncpu-1
    if (myid == icpu) then 
      if (icpu == 0) then
        open(newunit=outputfile, file=outfilename, action='write')
        write(outputfile, *) '#ion_x ion_y ion_KE_i ion_KE_f tar_KE tan_phi qout r_min tan_psi'
      else 
        open(newunit=outputfile, file=outfilename, position="append", status='old', action='write')
      end if 

      do ii=1, nchunk
        write(outputfile, *) ion_xy_arr(ii,1), ion_xy_arr(ii,2), ion_ke_arr(ii), ke_tar_arr(ii), ion_qout_arr(ii),&
          r_min_arr(ii), tan_phi_arr(ii), tan_psi_arr(ii), ion_qout_arr(ii)
      end do 
      close(outputfile)
    end if 
  end do 
  call MPI_FINALIZE(err)
end program 
