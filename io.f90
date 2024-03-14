module io
  use stdlib_kinds, only: sp, dp
  use stdlib_random, only: random_seed
  use stdlib_stats_distribution_uniform, only: uni => rvs_uniform
  use stdlib_stats_distribution_normal, only: norm => rvs_normal

  implicit none 
  real(dp), parameter :: tol=1.0e-15_dp,&
                         len_fact = 1.0_dp/0.529_dp,&
                         mass_fact=1822.89_dp,&
                         e_fact= 27.211_dp

  private 
  public :: init_random_seed,&
            read_indat,&
            print_indat,&
            set_potential,&
            setup_sim,&
            print_xyz,&
            load_target_fv,&
            count_lines,&
            init_target

  contains

  function count_lines(filename) result(nlines)
  character :: filename
  integer :: nlines, iofile, io
    open(iofile,file=filename)
    nlines = 0
    do
      read(iofile,*,iostat=io)
      nlines = nlines + 1
      if (io/=0) exit
    end do
    close(10)
  end function

  subroutine init_random_seed(consist, myid)
    integer :: seed_put, consist, seed_get, myid
    if (consist == 1) then
      seed_put = 666
    else
      call system_clock(count=seed_put)
    end if
    call random_seed(seed_put+myid, seed_get)
  end subroutine

  subroutine read_indat(fname_input, prename, ion_elem, ion_zz, ion_mass, ion_ke, &
    ion_qin, ff, gam_p, gam_c, gam_s, gam_cut, fwhm_qout, sigma_therm, frozen_par,&
    alpha_max, ion_zi, ion_zf,dx_step, acc, nhist, log_mode, surface_cov, v_typename,&
    fname_target)
    
    character(len=:), allocatable, intent(in) :: fname_input
    character(len=:), allocatable :: prename, ion_elem, fname_target, v_typename
    character(len=512) :: tmpstr
    real(dp) :: ion_zz, ion_mass, ion_ke, ff, fwhm_qout, sigma_therm, &
      frozen_par, alpha_max, ion_zi, ion_zf, dx_step, acc, surface_cov, &
      gam_p, gam_c, gam_s, gam_cut
    integer :: log_mode, nhist, io, stat, ion_qin
    
    open(newunit=io, file=fname_input, status='old', action='read')
      read(io,*) tmpstr
      prename = trim(tmpstr)
      read(io,*) tmpstr
      v_typename = trim(tmpstr)
      read(io,*) tmpstr 
      ion_elem = trim(tmpstr)
      read(io,*) ion_zz, ion_mass, ion_ke, ion_qin
      read(io,*) ff
      read(io,*) fwhm_qout, sigma_therm, frozen_par, alpha_max
      read(io,*) ion_zi, ion_zf, dx_step, acc, nhist
      read(io,*) gam_p, gam_c, gam_s, gam_cut
      read(io,*) log_mode, surface_cov
      read(io,*,iostat=stat) tmpstr
      if (stat == 0) then
        fname_target = trim(tmpstr)
      else 
        fname_target= 'structure.dat'
      end if 
    close(io)
  end subroutine

  subroutine print_indat(fname_input, prename, ion_elem, ion_zz, ion_mass, ion_ke, &
    ion_qin, ff, gam_p, gam_c, gam_s, gam_cut, &
    fwhm_qout, sigma_therm, frozen_par, alpha_max, ion_zi, ion_zf,dx_step, &
    acc, nhist, log_mode, surface_cov, v_typename, fname_target)
    
    character(len=:), allocatable, intent(in) :: fname_input, prename, ion_elem, fname_target, v_typename
    real(dp), intent(in) :: ion_zz, ion_mass, ion_ke, ff, fwhm_qout, sigma_therm, &
      frozen_par, alpha_max, ion_zi, ion_zf, dx_step, acc, surface_cov, gam_p, gam_c, gam_s, gam_cut
    integer, intent(in) :: log_mode, nhist, ion_qin
    
    print*, fname_input
    print*, prename
    print*, v_typename
    print*, ion_elem
    print*, ion_zz, ion_mass, ion_ke, ion_qin
    print*, ff
    print*, fwhm_qout, sigma_therm, frozen_par, alpha_max
    print*, ion_zi, ion_zf, dx_step, acc, nhist
    print*, gam_p, gam_c, gam_s, gam_cut
    print*, log_mode, surface_cov 
    print*, fname_target
  end subroutine

  function set_potential(vname) result(val)
    character(len=*), intent(in) :: vname
    integer :: val

    select case (vname)
    case ('hollow-krc')
      val = 1
    case default 
      val = 1
    end select
  end function

  subroutine get_random_ion_xy(ion_xy, cell)
    real(dp) :: ion_xy(2), cell(3)
    
    ion_xy(1) = uni(cell(1))
    ion_xy(2) = uni(cell(2))
  end subroutine
  
  subroutine init_target(fname, x, m, zz, v, a, natom) 
    character(len=:), allocatable, intent(in) :: fname
    real(dp), allocatable :: x(:,:), v(:,:), a(:,:), m(:), zz(:)
    integer :: io, natom
    open(newunit=io, file=fname, status='old', action='read')
      ! natom here is the number of target atoms
      read(io,*) natom
    close(io)
      allocate(x(natom+1,3))
      allocate(m(natom+1))
      allocate(zz(natom+1))
      allocate(v(natom+1,3))
      allocate(a(natom+1,3))
  end subroutine 

  subroutine load_target_fv(fname, sigma_therm, x, v, a, m, zz, natom, cell, cell_scaled)
    integer :: natom, io, i, ind
    character(len=:), allocatable, intent(in) :: fname
    real(dp), intent(in) :: sigma_therm
    real(dp) :: z_max, cell_scaled(3), cell(3), x(natom,3), v(natom,3), a(natom,3), m(natom), zz(natom)
    
    open(newunit=io, file=fname, status='old', action='read')
    
      ! natom here is the number of target atoms
      read(io,*) 
      
      z_max = -1.0_dp * huge(1.0_dp)
      read(io,*) cell(1), cell(2), cell(3)
      do i=2, natom
        read(io,*) x(i,1), x(i,2), x(i,3), zz(i), m(i)
        do ind=1, 3
          x(i,ind) = x(i,ind) * len_fact
        end do 
        m(i) = m(i) * mass_fact
        z_max = max(x(i,3),z_max)
      end do
    close(io)

    !shift target so that z=0 and add thermal motion if specified
    do i=2,natom
      if (abs(sigma_therm) > tol) then
        do ind=1,3
          x(i,ind) = norm(x(i,ind), sigma_therm) 
        end do 
      end if 
      x(i,3) = x(i,3) - z_max
    end do 

    v = 0.0_dp
    a = 0.0_dp
    cell_scaled = cell*len_fact
  end subroutine

  subroutine setup_sim(ion_zi, ion_zf, ion_xy, ion_zz, ion_m, ion_qin, ion_ke, a_pos, &
                       a_mass, a_zz, a_v, a_a, cell, cell_scaled, n_cor, n_sta, n_cap, &
                       factor, ff, r0, r_min, vp)
    real(dp), intent(in) :: ion_zi, ion_zf, ion_m, ion_zz, factor, ion_ke, cell(3), cell_scaled(3)
    real(dp) :: a_pos(:,:), a_mass(:), a_zz(:), ion_xy(2), delta(2), dir, &
      n_cor, n_sta, n_cap, ff, r0, r_min, a_v(:,:), a_a(:,:), vp
    integer :: natom, ion_qin, i, ind

    natom = size(a_mass)
    ! If we provide initial xy ion coordinates of negative then pick a random location
    if ((ion_xy(1) < 0.0_dp) .and. (ion_xy(2) < 0.0_dp)) then
      print*, 'Get random ion XY coordinate'
      call get_random_ion_xy(ion_xy, cell)
    end if

    ! setup ion as first atom in array
    a_pos(1,:2) = ion_xy *len_fact
    a_pos(1,3) = ion_zi  *len_fact
    a_mass(1) = ion_m    *mass_fact
    a_zz(1) = ion_zz

    ! Initialize ion electron counts
    n_cor = ion_zz - 1.0_dp * ion_qin
    n_sta = 0.0_dp
    n_cap = 0.0_dp
    ff = factor
    r0 = (1.81_dp + 1.6_dp * sqrt(1.0_dp*ion_qin)) * len_fact ! TD-DFT calc of graphene
    r_min = huge(1.0_dp)

    a_v = 0.0_dp
    a_a = 0.0_dp

    ! Set initial ion velocity from supplied kinetic energy
    vp = sqrt((2.0_dp * ion_ke * 1000.0_dp / e_fact)/(a_mass(1)))
    print*, ion_ke, vp, ion_m, a_mass(1)
    print*, 'E_k (eV/u) = ', ion_ke/(a_mass(1)/mass_fact), ' vp= ', vp
    
    dir = ion_zf - ion_zi
    a_v(1,3) = sign(vp, dir)

    ! place ion in the center of the target
    do i=2, natom
      delta(:) = a_pos(i,:2) - a_pos(1,:2)
      do ind=1,2
        if (abs(delta(ind))>(cell_scaled(ind)*0.5_dp)) then
          a_pos(i,ind) = a_pos(i,ind) - cell_scaled(ind)* dnint(delta(ind)/cell_scaled(ind))
        end if  
      end do 
    end do
    print*, 'X coord extent: ', minval(a_pos(:,1)), maxval(a_pos(:,1))
    print*, 'Y coord extent: ', minval(a_pos(:,2)), maxval(a_pos(:,2))
  end subroutine

  subroutine print_xyz(a_pos, a_mass, a_zz, cell, fname, status)
    real(dp) :: a_pos(:,:), cell(:), a_zz(:), a_mass(:), zero
    integer :: io, i, ind, natom
    character(len=*):: fname, status

    natom = size(a_mass)
    zero = 0.0_dp

    if (status == 'old') then
      open(newunit=io, file=fname, status=status, position='append', action='write')
    else 
      open(newunit=io, file=fname, action='write')
    end if 

      write(io,*) natom
      write(io,*) 'Lattice="',cell(1), zero, zero, &
                              zero , cell(2), zero, &
                              zero, zero, cell(3), &
                              '" Properties=pos:R:3:type:I:1:mass:R:1'
      do i=1, natom 
        write(io,*) a_pos(i,1)/len_fact, a_pos(i,2)/len_fact, a_pos(i,3)/len_fact, a_zz(i), a_mass(i)/mass_fact
      end do 
    close(io)
  end subroutine



end module
