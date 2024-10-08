module io
  use mod_types, only: dp

  implicit none 
  real(dp), parameter :: tol = 1.0e-15_dp, &
                         len_fact = 1.0_dp/0.529_dp, &
                         mass_fact = 1822.89_dp, & ! 1 amu = 1822.89 a.u.
                         e_fact = 27.211_dp, &
                         pi = 4.d0*datan(1.d0)

  private
  public :: str,&
            uni,&
            normal,&
            read_indat, &
            print_indat, &
            set_potential, &
            setup_sim, &
            print_xyz, &
            load_target_fv, &
            count_lines, &
            init_target, &
            export_e_number

  contains
  
  !   "Convert an integer to string."
  character(len=20) function str(k)
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
    str = trim(str)
  end function str

  real(dp) function uni(a,b)
    real(dp), intent(in) :: a, b
    call random_number(uni)
    uni = uni * (b-a) + a
  end function uni 

  real(dp) function normal(a,b)
  real(dp), intent(in) :: a, b
  real(dp) :: r1, r2, x

  call random_number(r1)
  call random_number(r2)
  x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )
  normal = a + b * x

  end function normal

  subroutine count_lines(filename, nlines)
    character(len=*), intent(in) :: filename
    integer :: nlines, iofile, stat

    nlines = 0
    open(newunit= iofile, file=filename, action='read')
    stat = 0
    do while (stat == 0)
      read(iofile, *, iostat=stat)
      if (stat == 0) nlines = nlines + 1
    end do
    close(iofile)
  end subroutine

!  subroutine init_random_seed(consist, myid)
!    integer :: seed_put, consist, seed_get, myid
!  
!    if (consist == 1) then
!      seed_put = 666
!    else
!      call system_clock(count=seed_put)
!    end if
!    call random_seed(seed_put + myid, seed_get)
!  end subroutine

  subroutine read_indat(fname_input, prename, ion_elem, ion_zz, ion_mass, ion_ke, &
    ion_qin, ff, gam_p, gam_c, gam_s, gam_cut, fwhm_qout, sigma_therm, frozen_par, &
    alpha_max, ion_zi, ion_zf, dx_step, acc, nhist, log_mode, is_xyz, v_typename, &
    chi_min, chi_max, fname_target)

    character(len=:), allocatable, intent(in) :: fname_input
    character(len=:), allocatable :: prename, ion_elem, fname_target, v_typename
    character(len=512) :: tmpstr
    real(dp) :: ion_zz, ion_mass, ion_ke, ff, fwhm_qout, sigma_therm, &
      frozen_par, alpha_max, ion_zi, ion_zf, dx_step, acc, &
      gam_p, gam_c, gam_s, gam_cut, &
      chi_min, chi_max
    integer :: log_mode, nhist, io, stat, ion_qin, is_xyz

    open(newunit=io, file=fname_input, status='old', action='read')
      read(io, *) tmpstr
      prename = trim(tmpstr)
      read(io, *) tmpstr
      v_typename = trim(tmpstr)
      read(io, *) tmpstr
      ion_elem = trim(tmpstr)
      read(io, *) ion_zz, ion_mass, ion_ke, ion_qin
      read(io, *) ff
      read(io, *) fwhm_qout, sigma_therm, frozen_par, alpha_max
      read(io, *) ion_zi, ion_zf, dx_step, acc, nhist
      read(io, *) chi_min, chi_max
      read(io, *) gam_p, gam_c, gam_s, gam_cut
      read(io, *) log_mode, is_xyz
      read(io, *, iostat=stat) tmpstr
      if (stat == 0) then
        fname_target = trim(tmpstr)
      else
        fname_target= 'structure.dat'
      end if
    close(io)
  end subroutine

  subroutine print_indat(fname_input, prename, ion_elem, ion_zz, ion_mass, ion_ke, &
    ion_qin, ff, gam_p, gam_c, gam_s, gam_cut, &
    fwhm_qout, sigma_therm, frozen_par, alpha_max, ion_zi, ion_zf, dx_step, &
    acc, nhist, log_mode, is_xyz, v_typename, fname_target)

    character(len=:), allocatable, intent(in) :: fname_input, prename, ion_elem, fname_target, v_typename
    real(dp), intent(in) :: ion_zz, ion_mass, ion_ke, ff, fwhm_qout, sigma_therm, &
      frozen_par, alpha_max, ion_zi, ion_zf, dx_step, acc, gam_p, gam_c, gam_s, gam_cut
    integer, intent(in) :: log_mode, nhist, ion_qin, is_xyz

    print *, fname_input
    print *, prename
    print *, v_typename
    print *, ion_elem
    print *, ion_zz, ion_mass, ion_ke, ion_qin
    print *, ff
    print *, fwhm_qout, sigma_therm, frozen_par, alpha_max
    print *, ion_zi, ion_zf, dx_step, acc, nhist
    print *, gam_p, gam_c, gam_s, gam_cut
    print *, log_mode, is_xyz
    print *, fname_target
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

  subroutine get_random_ion_xy(ion_x, ion_y, cell)
    real(dp) :: ion_x, ion_y, cell(3)

    ion_x = uni(0.0_dp, cell(1))
    ion_y = uni(0.0_dp, cell(2))
  end subroutine

  subroutine init_target(fname, natom) 
    character(len=:), allocatable, intent(in) :: fname
    integer :: io, natom, myid, ncpu

    open(newunit=io, file=fname, status='old', action='read')
      ! natom here is the number of target atoms
      read(io, *) natom
    close(io)
  end subroutine 

  subroutine load_target_fv(fname, x, y, z, m, zz, natom, cell, cell_scaled)
    integer :: natom, io, i, ind
    character(len=:), allocatable, intent(in) :: fname
    real(dp) :: z_max, cell_scaled(3), cell(3), x(natom), y(natom), z(natom), m(natom), zz(natom)

    open(newunit=io, file=fname, status='old', action='read')
      ! natom here is the number of target atoms
      read(io, *)

      z_max = -1.0_dp * huge(1.0_dp)
      read(io, *) cell(1), cell(2), cell(3)
      do i = 2, natom
        read(io, *) x(i), y(i), z(i), zz(i), m(i)
        x(i) = x(i) * len_fact
        y(i) = y(i) * len_fact
        z(i) = z(i) * len_fact
        m(i) = m(i) * mass_fact
        z_max = max(z(i), z_max)
      end do
    close(io)

    ! shift target so that z=0 and add thermal motion if specified
    do i = 2, natom
      z(i) = z(i) - z_max
    end do
    cell_scaled = cell * len_fact
  end subroutine

  subroutine setup_sim(ion_zi, ion_zf, ion_x, ion_y, ion_zz, ion_m, ion_qin, ion_ke, a_pos, &
                       a_mass, a_zz, a_v, a_a, cell, cell_scaled, n_cor, n_sta, n_cap, &
                       factor, ff, r0, r_min, vp, sigma_therm, &
                       chi_min, chi_max, chi)
    real(dp), intent(in) :: ion_zi, ion_zf, ion_m, ion_zz, factor, ion_ke, cell(3), cell_scaled(3), sigma_therm
    real(dp) :: a_pos(:,:), a_mass(:), a_zz(:), ion_x, ion_y, delta(2), dir, &
      n_cor, n_sta, n_cap, ff, r0, r_min, a_v(:,:), a_a(:,:), vp, &
      chi_min, chi_max, chi
    integer :: natom, ion_qin, i, ind

    natom = size(a_mass)
    ! if we provide initial xy ion coordinates of negative then pick a random location
    if ((ion_x < 0.0_dp) .and. (ion_y < 0.0_dp)) then
      call get_random_ion_xy(ion_x, ion_y, cell)
    end if

    ! setup ion as first atom in array
    a_pos(1,1) = ion_x * len_fact
    a_pos(1,2) = ion_y * len_fact
    a_pos(1,3) = ion_zi * len_fact
    a_mass(1) = ion_m * mass_fact
    a_zz(1) = ion_zz

    ! initialize ion electron counts
    n_cor = ion_zz - 1.0_dp * ion_qin
    n_sta = 0.0_dp
    n_cap = 0.0_dp
    ff = factor
    r0 = (1.81_dp + 1.6_dp * sqrt(1.0_dp * ion_qin)) * len_fact ! TD-DFT calc of graphene
    r_min = huge(1.0_dp)

    a_v = 0.0_dp
    a_a = 0.0_dp

    ! set initial ion velocity from supplied kinetic energy
    vp = sqrt((2.0_dp * ion_ke * 1000.0_dp / e_fact) / (a_mass(1)))

    dir = ion_zf - ion_zi

    if (chi_min == chi_max) then
      chi = chi_min
    else
      chi = uni(chi_min, chi_max - chi_min)
    end if

    if (chi > 45.0_dp) then
      chi = 45.0_dp
    end if

    ! from degree to rad
    chi = chi*pi/180.0_dp

    ! calculate velocity vector
    a_v(1,1) = -1.0_dp*sign(vp, dir)*sin(chi)
    a_v(1,3) = sign(vp, dir)*cos(chi)

    ! place ion in the center of the target
    do i = 2, natom
      delta(:) = a_pos(i,:2) - a_pos(1,:2)
      do ind = 1, 2
        if (abs(delta(ind)) > (cell_scaled(ind) * 0.5_dp)) then
          a_pos(i,ind) = a_pos(i,ind) - cell_scaled(ind) * dnint(delta(ind)/cell_scaled(ind))
        end if
      end do
      ! rotate target atoms
      !a_pos(i,1) = a_pos(i,1)*cos(omega) - a_pos(i,2)*sin(omega)
      !a_pos(i,2) = a_pos(i,1)*sin(omega) + a_pos(i,2)*cos(omega)
    end do

    ! translate ion to correct coordinates
    a_pos(1,1) = a_pos(1,1) - a_pos(1,3)*tan(chi)

    ! rotate ion
    !a_pos(1,1) = a_pos(1,1)*cos(omega) - a_pos(1,2)*sin(omega)
    !a_pos(1,2) = a_pos(1,1)*sin(omega) + a_pos(1,2)*cos(omega)

    ! from rad to degree
    chi = chi/pi*180.0_dp

    if (sigma_therm > 0.0) then
      do i = 2, natom
        do ind = 1, 3
          a_pos(i,ind) = a_pos(i,ind) + normal(a_pos(i,ind), sigma_therm)
        end do
      end do
    end if
  end subroutine

  subroutine print_xyz(a_pos, a_mass, a_zz, cell, time, fname, status)
    real(dp) :: a_pos(:,:), cell(:), a_zz(:), a_mass(:), zero, time
    integer :: io, i, ind, natom
    character(len=*) :: fname, status

    natom = size(a_mass)
    zero = 0.0_dp

    if (status == 'old') then
      open(newunit=io, file=fname, status=status, position='append', action='write')
    else
      open(newunit=io, file=fname, action='write')
    end if

      write(io, *) natom
      write(io, *) 'Lattice="',cell(1), zero, zero, &
                              zero , cell(2), zero, &
                              zero, zero, cell(3), '"',&
                              ' Time=', time, &
                              ' Properties=pos:R:3:type:I:1:mass:R:1'
      do i=1, natom 
        write(io,*) a_pos(i,1)/len_fact, a_pos(i,2)/len_fact, a_pos(i,3)/len_fact, a_zz(i), a_mass(i)/mass_fact
      end do 
    close(io)
  end subroutine

  subroutine export_e_number(n_cor, n_sta, n_cap, count, a_pos_z, method, i_ion)
    real(dp) :: n_cor, n_sta, n_cap, a_pos_z
    integer :: count, io, method, i_ion
    logical :: exists

    if (count == 0) then
      inquire(file="logs/e-number_"//trim(str(i_ion))//".txt", exist=exists)
      if (exists) then
        open(newunit=io, file="logs/e-number_"//trim(str(i_ion))//".txt", status="replace", action="write")
          write(io, *) a_pos_z, n_cor, n_sta, n_cap
        close(io)
      else
        open(newunit=io, file="logs/e-number_"//trim(str(i_ion))//".txt", status="new", action="write")
          write(io, *) a_pos_z, n_cor, n_sta, n_cap
        close(io)
      end if
    else
      open(newunit=io, file="logs/e-number_"//trim(str(i_ion))//".txt", position="append", action="write")
        write(io, *) a_pos_z, n_cor, n_sta, n_cap
      close(io)
    end if
  end subroutine

end module
