module force
  use stdlib_kinds, only: sp, dp
  use potentials, only: calc_vprime
  use stdlib_stats_distribution_normal, only: norm_samp => rvs_normal

  implicit none
  real(dp), parameter :: tol=1.0e-15_dp,&
                         len_fact=1.0_dp/0.529_dp,&
                         mass_fact=1822.89_dp,&
                         e_fact= 27.211_dp
  
  private
  public :: varystep, update_ion, calc_ion_props, update_ion_old

  contains 

  function dist(x0,x1,pbc) result(val)
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
      !  delta(ind) = delta(ind) - pbc(ind) * dnint(delta(ind)/pbc(ind))
      !end if 

    !  val = val + delta(ind)**2
    !end do
    val = sqrt(sum((x1-x0)**2))
  end function 

  function lambda_f(z, r0, froz_p, a, mu) result(val)
    real(dp), intent(in) :: z, r0, froz_p, a, mu
    real(dp) :: val
    val = (froz_p * a)/4.0_dp * (erf(mu*(z+r0))+1.0_dp) * (erf(mu*(r0-z))+1.0_dp)
  end function

  function gamma_f(r, p, c, s) result(val)
    real(dp), intent(in) :: r, p, c, s
    real(dp) :: val
    val = abs(-1.0_dp * 0.5_dp * p * (erf(s*(r-c))-1.0_dp))
  end function

  subroutine varystep(t, x, v, a, mass, z, cell, ion_iv, acc, dt_max, &
      vtype, ff, ddr, n_cap, n_sta, n_cor, r_cut, r0, dt, verbose)
    real(dp), intent(in) :: acc, dt_max, ion_iv, ddr, r0
    real(dp) :: x(:,:), v(:,:), a(:,:), &
                mass(:), z(:), cell(:), dt, dt2, acel, tmp, t, ff, &
                n_cap, n_sta, n_cor, r_cut, r
    real(dp), allocatable :: xold(:,:), vold(:,:), vprime(:)
    integer :: i, j, ind, natom, vtype, verbose
    natom = size(mass)

    if (.not. allocated(xold)) allocate(xold(natom,3))
    if (.not. allocated(vold)) allocate(vold(natom,3))
    if (.not. allocated(vprime)) allocate(vprime(natom))
    
    xold = x
    vold = v
    
    acel = sqrt(sum(a(1,:)**2))+tiny(1.0_dp)
    dt =min(dt_max, acc * ion_iv/acel)
    dt2 = dt * 0.5_dp

    do i=1, natom
      ! test with pascal code
      !a(i, :) = 0.0_dp
      if (i == 1) then 
        ! for similarity 
        a(1, :) = 0.0_dp
        do j=2, natom
          r = dist(x(i,:), x(j,:), cell) 
          tmp = calc_vprime(vtype, r, r0, r_cut, n_sta, n_cor, n_cap, ff, z(1), z(j), ddr)
          vprime(j) = -1.0_dp * tmp/r
          do ind=1, 3
            a(1,ind) = a(1,ind) + (x(1,ind) - x(j,ind)) * vprime(j)/mass(1)
          end do 
        end do
      else
        do ind=1, 3
          a(i, ind) = (x(i,ind)-x(1,ind))*vprime(i)/mass(i)
        end do
      end if 

      do ind=1, 3
        v(i,ind) = vold(i,ind) + dt2 * a(i,ind)
        x(i,ind) = xold(i,ind) + 0.5_dp * (v(i,ind) + vold(i,ind))*dt2
      end do 
    end do 

    do i=1, natom 
      ! test with pascal code
      !a(i, :) = 0.0_dp
      if (i == 1) then 
        a(1, :) = 0.0_dp
        do j=2, natom
          r = dist(x(1,:), x(j,:), cell) 
          tmp = calc_vprime(vtype, r, r0, r_cut, n_sta, n_cor, n_cap, &
                            ff, z(1), z(j), ddr)
          vprime(j) = -1.0_dp * tmp/r
          do ind=1, 3
            a(1,ind) = a(1,ind) + (x(1,ind) - x(j,ind)) * vprime(j)/mass(1)
          end do 
        end do
      else
        do ind=1, 3
          a(i, ind) = (x(i,ind)-x(1,ind))* vprime(i)/mass(i)
        end do
      end if 

      do ind=1, 3
        v(i,ind) = vold(i,ind) + dt * a(i,ind)
        x(i,ind) = xold(i,ind) + 0.5 * (v(i,ind) + vold(i,ind))*dt
      end do 
    end do 
  end subroutine

  subroutine update_ion(dt, t, x, z, n_sta, n_cap, n_cor, factor, ff, r0, ion_iv, &
      lam_a, froz_p, lam_mu, alpha_max, r_min, gam_p, gam_c, gam_s, gam_cut, cell, &
      logfile, verbose)
    real(dp), intent(in) :: x(:,:), cell(:), r0, lam_a, froz_p, lam_mu, dt, t, &
      ion_iv, factor, alpha_max, gam_p, gam_c, gam_s, gam_cut
    real(dp) :: z(:), n_sta, n_cap, n_cor, lam, tmp, gam, ff, r_min, r, & 
      tmp_n_cap, tmp_n_sta, blah
    integer :: logfile, j, natom, verbose, nneigh, nelectron

    natom = size(x(:,1))
    nelectron = dint(n_cor+n_sta+n_cap)
    
    lam = lambda_f(x(1,3), r0, froz_p, lam_a, lam_mu) 
    r = huge(1.0_dp)

    tmp_n_sta = 0.0
    tmp_n_cap = 0.0

    ! calculate distances between ion and target atoms
    do j=2, natom
      tmp = dist(x(1,:), x(j,:), cell) 
      r = min(tmp, r)
      ! if the gamma cutoff is +ve, calculate gamma and electron transfer between all targets within gam_cut
      if ((tmp < gam_cut) .and. (gam_cut>0.0) .and. nelectron < z(1)) then
        gam = gamma_f(tmp, gam_p, gam_c, gam_s)
        tmp_n_cap = tmp_n_cap + lam*dt*(z(1) - n_cor - n_cap - n_sta) - gam*dt*n_cap
        tmp_n_sta = tmp_n_sta + gam *dt*n_cap
      end if 
    end do

    ! else if gam_cut is -ve, then just calculate electron transfer between ion and nearest target atom
    if (gam_cut <= 0.0 .and. nelectron< z(1)) then
      gam = gamma_f(r, gam_p, gam_c, gam_s)
      tmp_n_cap = lam*dt*(z(1) - n_cor - n_cap - n_sta) - gam*dt*n_cap
      tmp_n_sta = gam *dt*n_cap
    end if 

    ! update the ion electron counts 
    if (nelectron< z(1)) then
      blah = z(1)-n_cor-n_sta-n_cap
      ! if ion wants more electrons than Z, then scale this so that electrons <= Z
      if (tmp_n_sta + tmp_n_cap > blah) then 
        tmp_n_sta = tmp_n_sta/blah * (tmp_n_sta + tmp_n_cap)
        tmp_n_cap = tmp_n_cap/blah * (tmp_n_sta + tmp_n_cap) 
      end if 
    end if 
    n_sta = n_sta + tmp_n_sta
    n_cap = n_cap + tmp_n_cap


    ! record closesest approach distance 
    if (r<r_min) r_min = r

    print*, 'tmp_sta: ', tmp_n_sta, '  tmp_cap: ',tmp_n_cap
    print*, 'Core: ', n_cor, '  Captured: ', n_cap, '  Stabilised: ', n_sta, ' blah:', blah

    ! old way
    !do j=2, natom
    !  tmp = dist(x(1,:), x(j,:), cell)
    !  r = min(tmp, r)
    !end do 
    !if (r<r_min) r_min = r

    !gam = gamma_f(r, gam_p, gam_c, gam_s)
    !tmp = n_cap + lam*dt*(z(1) - n_cor-n_cap-n_sta) - gam*dt*n_cap
    !n_sta = n_sta + gam*dt*n_cap
    !n_cap = tmp

    tmp = exp(alpha_max * ion_iv * t)
    ff = tmp/(1.0_dp +tmp) + factor/(1.0_dp+tmp)
  end subroutine


  subroutine update_ion_old(dt, t, x, z, n_sta, n_cap, n_cor, factor, ff, r0, ion_iv, &
      lam_a, froz_p, lam_mu, alpha_max, r_min, gam_a, gam_b, gam_c, cell, &
      logfile, verbose)
    real(dp), intent(in) :: x(:,:), cell(:), r0, lam_a, froz_p, lam_mu, dt, t, &
      ion_iv,factor, alpha_max, gam_a, gam_b, gam_c
    real(dp) :: z(:), n_sta, n_cap, n_cor, lam, tmp, gam, ff, r_min, r
    integer :: logfile, j, natom, verbose, nneigh

    natom = size(x(:,1))

    lam = lambda_f(x(1,3), r0, froz_p, lam_a, lam_mu)

    r = huge(1.0_dp)
    nneigh = 0
    do j=2, natom
      tmp = dist(x(1,:), x(j,:), cell)
      if (tmp/len_fact < 1.5) nneigh = nneigh + 1
      r = min(tmp, r)
    end do
    !print*, 'Ion has ', nneigh, ' neighbours'

    if (r<r_min) r_min = r

    gam = gamma_f(r, gam_a, gam_b, gam_c)

    tmp = n_cap + lam*dt*(z(1) - n_cor - n_cap - n_sta) - gam*dt*n_cap
    n_sta = n_sta + gam *dt*n_cap
    n_cap = tmp

    tmp = exp(alpha_max * ion_iv * t)
    ff = tmp/(1.0_dp +tmp) + factor/(1.0_dp+tmp)
  end subroutine
  
  subroutine calc_ion_props(fwhm_qout, a_v, a_zz, ion_iv, ion_qin, n_cor, n_sta, n_cap, tan_phi, tan_psi, ion_qout)
    real(dp), intent(in) :: a_v(:,:), a_zz(:), ion_iv(:)
    integer, intent(in) :: ion_qin
    real(dp) :: fwhm_qout, n_cor, n_sta, n_cap, tan_phi, tan_psi, blah, rand
    integer :: ion_qout

    blah = a_zz(1)-n_cor-n_sta 
    rand = norm_samp(blah, fwhm_qout/2.355_dp)

    ! if fwhm_qout == 0
    if (abs(fwhm_qout) < tol) then
      ion_qout = max(0, nint( blah ))
    else if (max(0,nint(rand))< ion_qin) then 
      ion_qout = max(0, nint(rand))
    else
      ion_qout = ion_qin
    end if

    tan_phi = huge(1.0_dp)
    tan_psi = sqrt(a_v(1,1)**2 + a_v(1,2)**2)/a_v(1,3)
    if (abs(a_v(1,3)) > tol) tan_phi = a_v(1,1)/a_v(1,2)
    print*, tan_phi, a_v(1,2)/a_v(1,1)
  end subroutine 

end module
