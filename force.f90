module force
  use stdlib_kinds, only: sp, dp
  use stdlib_strings, only: to_string
  use potentials, only: calc_vprime
  use stdlib_stats_distribution_normal, only: norm_samp => rvs_normal

  implicit none
  real(dp), parameter :: tol = 1.0e-15_dp, &
                         len_fact = 1.0_dp/0.529_dp, &
                         mass_fact = 1822.89_dp, &
                         e_fact = 27.211_dp, &
                         force_fact = 8.2387235038_dp * 10.0e-8_dp ! E_h/a_0, from a.u. to SI
  
  private
  public :: varystep, update_ion, calc_ion_props, update_ion_old

  contains

  function dist(x0, x1, pbc) result(val)
    real(dp), intent(in) :: x0(:), x1(:), pbc(:)
    real(dp) :: val, delta(3), half_pbc(3), a(3), b(3)
    integer :: ind

    val = 0.0_dp
    !half_pbc = pbc * 0.5_dp
    !a = x0 - half_pbc
    !b = x1 - half_pbc

    !do ind=1, 3
      !delta(ind) = x0(ind) - x1(ind)
      !if (ind < 3) then
      !  delta(ind) = delta(ind) - pbc(ind) * dnint(delta(ind) / pbc(ind))
      !end if

      !val = val + delta(ind)**2
    !end do
    val = sqrt(sum((x1 - x0)**2))
  end function 

  function lambda_f(z, r0, froz_p, a, mu) result(val)
    real(dp), intent(in) :: z, r0, froz_p, a, mu
    real(dp) :: val
    val = (froz_p * a) / 4.0_dp * (erf(mu*(z+r0)) + 1.0_dp) * (erf(mu*(r0-z)) + 1.0_dp)
  end function

  function gamma_f(r, p, c, s) result(val)
    real(dp), intent(in) :: r, p, c, s
    real(dp) :: val
    val = abs(-1.0_dp * 0.5_dp * p * (erf(s*(r-c)) - 1.0_dp))
  end function

  ! x = a_pos, v = a_vel, a = a_acel, mass = a_mass, z = a_zz
  ! cell = cell_scaled, ion_iv = vp, vtype = v_type
  subroutine varystep(t, x, v, a, mass, z, cell, ion_iv, acc, dt_max, &
      vtype, ff, ddr, n_cap, n_sta, n_cor, r_cut, r0, dt, verbose, count, x_init, method, i_ion, springs, export_files)
    real(dp), intent(in) :: acc, dt_max, ion_iv, ddr, r0
    real(dp) :: x(:,:), v(:,:), a(:,:), &
                mass(:), z(:), cell(:), dt, dt2, acel, tmp, t, ff, &
                n_cap, n_sta, n_cor, r_cut, r, e_pot, e_kin_ion, k, x_init(:,:)
    real(dp), allocatable :: xold(:,:), vold(:,:), aold(:,:), vprime(:), potential_energies(:), kinetic_energies(:), &
                             e_kin_target(:)
    integer :: run, run_end, i, j, ind, natom, vtype, verbose, method, io, count, i_ion
    logical :: exists, springs, export_files
    natom = size(mass)

    if (.not. allocated(xold)) allocate(xold(natom,3))
    if (.not. allocated(vold)) allocate(vold(natom,3))
    if (.not. allocated(vprime)) allocate(vprime(natom))
    if (.not. allocated(potential_energies)) allocate(potential_energies(natom))
    if (.not. allocated(kinetic_energies)) allocate(kinetic_energies(natom))
    if (.not. allocated(e_kin_target)) allocate(e_kin_target(natom))

    xold = x
    vold = v
    aold = a

    ! k = 7 mdyn/Angstrom = 7 * 10**-8 N/Angstrom
    k = 7.0_dp * 10.0e-8_dp / force_fact * len_fact

    acel = sqrt(sum(a(1,:)**2)) + tiny(1.0_dp) ! + tiny(1.0_dp) to prevent division by zero
    dt = min(dt_max, acc * ion_iv / acel) ! dt_max = abs(dx_step/vp)
    dt2 = dt * 0.5_dp

    if (method == 1) then
      run_end = 2
    else
      run_end = 1
    end if

    do run = 1, run_end
      ! i = 1 -> ion, i > 1 -> target atoms
      do i = 1, natom
        a(i,:) = 0.0_dp
        ! ion
        if (i == 1) then
          do j = 2, natom
            r = dist(x(i,:), x(j,:), cell) ! distance between ion i and target atom j
            ! tmp: F = W/d, [F] = kg*m/s**2
            tmp = calc_vprime(vtype, r, r0, r_cut, n_sta, n_cor, n_cap, ff, z(1), z(j), ddr, e_pot)
            potential_energies(j) = e_pot
            ! vprime: F/l... force per length, [F/l] = kg/s**2
            vprime(j) = -1.0_dp * tmp / r
            do ind = 1, 3
              ! [dx * vprime / mass] = m * kg/s**2 / kg = m/s**2
              a(1,ind) = a(1,ind) + (x(1,ind) - x(j,ind)) * vprime(j) / mass(1)
            end do
          end do
        ! target atoms
        else
          do ind = 1, 3
            a(i,ind) = (x(i,ind) - x(1,ind)) * vprime(i) / mass(i)
            if (springs .eqv. .true.) then
              a(i,ind) = a(i,ind) - k * (x(i,ind) - x_init(i,ind)) / mass(i)
            end if
          end do
        end if

        if (run == 1) then
          if (method == 1) then
            ! half-step algorithm
            do ind = 1, 3
              v(i,ind) = vold(i,ind) + dt2 * a(i,ind)
              x(i,ind) = xold(i,ind) + 0.5_dp * (v(i,ind) + vold(i,ind)) * dt2
            end do
          else
            ! velocity verlet
            do ind = 1, 3
              v(i,ind) = vold(i,ind) + 0.5 * (a(i,ind) + aold(i,ind)) * dt
              x(i,ind) = xold(i,ind) + vold(i,ind) * dt + 0.5 * aold(i,ind) * dt**2
            end do
          end if
        else
          ! midpoint algorithm
          do ind = 1, 3
            v(i,ind) = vold(i,ind) + dt * a(i,ind)
            x(i,ind) = xold(i,ind) + 0.5 * (v(i,ind) + vold(i,ind)) * dt
          end do
        end if

        if (i == 1) then
          e_kin_ion = 0.5_dp * mass(i) * sum(v(i,:)**2) * e_fact
        else
          e_kin_target(i) = 0.5_dp * mass(i) * sum(v(i,:)**2) * e_fact
        end if

        kinetic_energies(i) = 0.5_dp * mass(i) * sum(v(i,:)**2) * e_fact

      end do
    end do

    if ((verbose > 0) .and. (export_files .eqv. .true.)) then
      if (count == 0) then
        inquire(file="log_energies_"//to_string(method)//"_"//to_string(i_ion)//".txt", exist=exists)
        if (exists) then
          open(newunit=io, file="log_energies_"//to_string(method)//"_"//to_string(i_ion)//".txt", status="replace", action="write")
            write(io, *) x(1,3), sum(potential_energies), sum(kinetic_energies), &
            sum(potential_energies) + sum(kinetic_energies), e_kin_ion, sum(e_kin_target)
          close(io)
        else
          open(newunit=io, file="log_energies_"//to_string(method)//"_"//to_string(i_ion)//".txt", status="new", action="write")
            write(io, *) x(1,3), sum(potential_energies), sum(kinetic_energies), &
            sum(potential_energies) + sum(kinetic_energies), e_kin_ion, sum(e_kin_target)
          close(io)
        end if
      else
        open(newunit=io, file="log_energies_"//to_string(method)//"_"//to_string(i_ion)//".txt", position="append", action="write")
          write(io, *) x(1,3), sum(potential_energies), sum(kinetic_energies), &
          sum(potential_energies) + sum(kinetic_energies), e_kin_ion, sum(e_kin_target)
        close(io)
      end if
    end if

  end subroutine

  subroutine update_ion(dt, t, x, z, n_sta, n_cap, n_cor, factor, ff, r0, ion_iv, &
      lam_a, froz_p, lam_mu, alpha_max, r_min, gam_p, gam_c, gam_s, gam_cut, cell, &
      logfile, verbose)
    real(dp), intent(in) :: x(:,:), cell(:), r0, lam_a, froz_p, lam_mu, dt, t, &
      ion_iv, factor, alpha_max, gam_p, gam_c, gam_s, gam_cut
    real(dp) :: z(:), n_sta, n_cap, n_cor, lam, tmp, gam, ff, r_min, r, &
      tmp_n_cap, tmp_n_sta, blah, tmp_gamma
    integer :: logfile, j, natom, verbose, nneigh, nelectron

    natom = size(x(:,1))
    nelectron = dint(n_cor + n_sta + n_cap)

    r = huge(1.0_dp)

    tmp_n_sta = 0.0_dp
    tmp_n_cap = 0.0_dp
    gam = 0.0_dp

    ! calculate distances between ion and target atoms
    do j = 2, natom
      tmp = dist(x(1,:), x(j,:), cell)
      r = min(tmp, r)
      ! if the gamma cutoff is +ve, calculate gamma between all targets within gam_cut
      if (gam_cut > 0.0_dp) then
        if (tmp < gam_cut) then
          gam = gam + gamma_f(tmp, gam_p, gam_c, gam_s)
        end if
      end if
    end do

    ! else if gam_cut is -ve, then just between ion and nearest target atom
    if (gam_cut <= 0.0) then
      gam = gamma_f(r, gam_p, gam_c, gam_s)
    end if

    lam = lambda_f(r, r0, froz_p, lam_a, lam_mu)
    n_cap = n_cap + lam * dt * (z(1) - n_cor - n_cap - n_sta) - gam * dt * n_cap
    n_sta = n_sta + gam * dt * n_cap

    ! record closest approach distance
    if (r < r_min) r_min = r

    if (verbose == 2) then
      print*, 'tmp_sta: ', tmp_n_sta, '  tmp_cap: ',tmp_n_cap
      print*, 'Core: ', n_cor, '  Captured: ', n_cap, '  Stabilised: ', n_sta, ' blah:', blah
    end if 

    tmp = exp(alpha_max * ion_iv * t)
    ff = tmp / (1.0_dp + tmp) + factor / (1.0_dp + tmp)
  end subroutine

  subroutine update_ion_old(dt, t, x, z, n_sta, n_cap, n_cor, factor, ff, r0, ion_iv, &
      lam_a, froz_p, lam_mu, alpha_max, r_min, gam_a, gam_b, gam_c, cutoff, cell, &
      logfile, verbose)
    real(dp), intent(in) :: x(:,:), cell(:), r0, lam_a, froz_p, lam_mu, dt, t, &
      ion_iv,factor, alpha_max, gam_a, gam_b, gam_c, cutoff
    real(dp) :: z(:), n_sta, n_cap, n_cor, lam, tmp, gam, ff, r_min, r
    integer :: logfile, j, natom, verbose, nneigh

    natom = size(x(:,1))

    r = huge(1.0_dp)
    nneigh = 0
    do j = 2, natom
      tmp = dist(x(1,:), x(j,:), cell)
      if (tmp / len_fact < 1.5) nneigh = nneigh + 1
      r = min(tmp, r)
    end do
    lam = lambda_f(r, r0, froz_p, lam_a, lam_mu)
    gam = gamma_f(r, gam_a, gam_b, gam_c)

    if (r < r_min) r_min = r

    tmp = n_cap + lam * dt * (z(1) - n_cor - n_cap - n_sta) - gam * dt * n_cap
    n_sta = n_sta + gam * dt * n_cap
    n_cap = tmp

    if (verbose == 2) then
      print *, 'Core: ', n_cor, '  Captured: ', n_cap, '  Stabilised: ', n_sta
    end if

    tmp = exp(alpha_max * ion_iv * t)
    ff = tmp / (1.0_dp + tmp) + factor / (1.0_dp + tmp)
  end subroutine

  subroutine calc_ion_props(fwhm_qout, a_v, a_zz, ion_iv, ion_qin, n_cor, n_sta, n_cap, tan_phi, tan_psi, ion_qout)
    real(dp), intent(in) :: a_v(:,:), a_zz(:), ion_iv(:)
    integer, intent(in) :: ion_qin
    real(dp) :: fwhm_qout, n_cor, n_sta, n_cap, tan_phi, tan_psi, blah, rand
    integer :: ion_qout

    blah = a_zz(1) - n_cor - n_sta
    rand = norm_samp(blah, fwhm_qout / 2.355_dp)

    ! if fwhm_qout == 0
    if (abs(fwhm_qout) < tol) then
      ion_qout = max(0, nint( blah ))
    else if (max(0,nint(rand))< ion_qin) then
      ion_qout = max(0, nint(rand))
    else
      ion_qout = ion_qin
    end if

    tan_phi = huge(1.0_dp)
    tan_psi = sqrt(a_v(1,1)**2 + a_v(1,2)**2) / a_v(1,3)
    if (abs(a_v(1,2)) > tol) tan_phi = a_v(1,1) / a_v(1,2)
  end subroutine

end module
