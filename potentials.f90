module potentials
  use stdlib_kinds, only: sp, dp
  
  implicit none

  private 
  public calc_vprime

  contains 

  function v_hollow_krc(r, r0, r_cut, n_sta, n_cor, n_cap, ff, zi, zj) result(val)
    real(dp), intent(in) :: r, r0, r_cut, n_sta, n_cor, n_cap, ff, zi, zj
    real(dp) :: a1, a2, phi1, phi2, phi3, tmp, x, x0, x1, x2, val

     a1   = 0.8854_dp/( (n_cor+n_sta)**0.23_dp + (zj/ff)**0.23_dp)
     x    = r/a1
     phi1 =   0.190945_dp*exp(-0.278544_dp*x)
     phi2 =   0.473674_dp*exp(-0.637174_dp*x)
     phi3 =   0.335381_dp*exp(-1.919249_dp*x)
     tmp = (phi1 + phi2 + phi3)*(n_cor + n_sta)* zj/r

     a2 = 0.8854_dp/((zj/ff)**(1/3))

     x = r/a2
     phi1 =   0.190945_dp*exp(-0.278544_dp*x)
     phi2 =   0.473674_dp*exp(-0.637174_dp*x)
     phi3 =   0.335381_dp*exp(-1.919249_dp*x)
     tmp = tmp + (phi1 + phi2 + phi3)*(zi - n_sta - n_cor - n_cap)* zj/r

     x0 = r0/a2
     x1 = min(x,x0)
     x2 = max(x,x0)
     phi1= phi1-0.190945_dp*exp(-0.278544_dp*x2)*sinh(0.278544_dp*x1)/(0.278544_dp*x1)*x/x2
     phi2= phi2-0.473674_dp*exp(-0.637174_dp*x2)*sinh(0.637174_dp*x1)/(0.637174_dp*x1)*x/x2
     phi3= phi3-0.335381_dp*exp(-1.919249_dp*x2)*sinh(1.919249_dp*x1)/(1.919249_dp*x1)*x/x2

     val = tmp + (phi1 + phi2 + phi3) * n_cap * zj/r
  end function

  function calc_vprime(vtype, r, r0, r_cut, n_sta, n_cor, n_cap, ff, zi, zj, ddr) result(val)
    real(dp), intent(in) :: r, r0, n_sta, n_cor, n_cap, ff, zi, zj, r_cut, ddr
    integer, intent(in) :: vtype
    real(dp) :: val, v1, v2, dr

    if (r > r_cut) then
      val = 0.0_dp
    else 
      dr = ddr*r
      select case (vtype)
        case (1)
          v1 = v_hollow_krc(r+dr, r0, r_cut, n_sta, n_cor, n_cap, ff, zi, zj)
          v2 = v_hollow_krc(r, r0, r_cut, n_sta, n_cor, n_cap, ff, zi, zj)

        case default 
          v1 = v_hollow_krc(r+dr, r0, r_cut, n_sta, n_cor, n_cap, ff, zi, zj)
          v2 = v_hollow_krc(r, r0, r_cut, n_sta, n_cor, n_cap, ff, zi, zj)
      end select

      val = (v1-v2)/dr
    end if
  end function
end module
