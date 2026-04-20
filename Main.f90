program main_project
  !use omp_lib
  !use omp_lib_kinds
  use constants
  use tight_binding
  implicit none

  complex(kind=dp), dimension(:,:), allocatable:: htn
  type(unit_cell)                              :: cell
  integer                                      :: sys_size,istat=0,i=1,&!j=1,&
                                                 &filling,size
  real(kind=dp)                                :: a_val, k_val,phase,&
                                                 &theta_val, phi=0.0_dp!,k_pos
  real(kind=dp), dimension(:),allocatable      :: egn,loc_temp
  real(kind=dp), dimension(:,:),allocatable    :: dat_array,grad_array,grad_sum,loc_array
  complex(kind=dp)                             :: kexp, aexp, pexp!,ctn
  real(kind=dp)                                :: d_phi,fermi
  integer                                      :: phi_max,k_max
  ! Variables which are used for the many body calculation are defined under here
  logical                                      :: mb=.false. ! Defines if the system is doing a many body calculation
  integer(kind=i64)                            :: d_tot
  integer(kind=i64),dimension(:,:),allocatable :: hbt_tot
  real(kind=dp)                                :: U_val=1.0_dp
  integer,dimension(2)                         :: many_body
  complex(kind=dp), dimension(:,:), allocatable:: hmb

  call init_sys(cell,sys_size,filling,phase,a_val,phi_max,&
               &many_body,U_val,k_max)

  if(sum(many_body)>0) mb=.true.

  if(mb)then

    print*, 'Many body calculation'

    size=sys_size*cell%n_sites

    ! NB this is an arbitrary safety feature
    ! I have not done the calculation to determine this
    if(size>=10) stop 'System is too large'

    call get_basis(2,2,size,hbt_tot,d_tot)

    k_val=phase*real_pi

        size=sys_size*cell%n_sites

    if(filling>size.or.filling==0) then 
      print*, 'invalid filling value, 1/2 filled current calculated'
      filling=int(real(size,kind=dp)/2.0_dp)
    end if
    
    k_val=phase*real_pi

    allocate(hmb(d_tot, d_tot), stat=istat)
    if(istat/=0) stop 'error allocating hmb matrix'

    allocate(loc_array(1:phi_max+1,d_tot+1), stat=istat)
    if(istat/=0) stop 'error allocating loc_array array'

    allocate(loc_temp(d_tot),stat=istat)
    if(istat/=0) stop 'error allocating loc_temp array'

    allocate(dat_array(1:phi_max+1,size+1),stat=istat)
    if(istat/=0) stop 'error allocating dat_array'

    allocate(grad_array(1:phi_max+1,size+1),stat=istat)
    if(istat/=0) stop 'error allocating grad_array'

    grad_array=0.0_dp

    aexp=exp(cmplx_i*a_val)

    do i=1,phi_max+1

      print*, 'iteration', i

      ! System now iterates -pi -> pi
      phi=2.0_dp*real(i-1,kind=dp)/real(phi_max,kind=dp)-1.0_dp

      !defining phi as $\frac{\phi}{\phi_0}$
      theta_val=((real_pi*2)*phi)
      pexp=exp(cmplx_i*theta_val)

        kexp=aexp**(k_val)

        call make_hmb(sys_size, cell, kexp, pexp, U_val, d_tot, hbt_tot, hmb)

        call zheev_evals(hmb,egn)

        dat_array(i,1) =phi
        grad_array(i,1)=phi
        loc_array(i,1)=phi
        dat_array(i,2:filling+1)=egn(1:filling) ! Only want the lower evals

        call get_local(hmb,d_tot,loc_temp)
        
        loc_array(i,2:)=loc_temp(:)

        deallocate(hmb, stat=istat)
        if(istat/=0) stop 'error deallocating hmb array'

        deallocate(egn, stat=istat)
        if(istat/=0) stop 'error deallocating egn array'
  
    end do

    d_phi = (2.0_dp * real_pi) / real(phi_max, kind=dp)

    call get_grad(dat_array(:, 2:filling+1), grad_array(:, 2:filling+1),&
                &d_phi, (phi_max+1), filling)

    allocate(grad_sum(1:phi_max+1,2),stat=istat)
    if(istat/=0) stop 'error allocating grad_sum'

    fermi=maxval(dat_array(:,filling))

    grad_sum(:,1)=grad_array(:,1)
    do i=1,phi_max+1
      grad_sum(i,2) =sum(grad_array(i,2:filling+1))
    !  dat_array(i,2:)=dat_array(i,2:)-fermi
    end do

    grad_sum(:,2:)=grad_sum(:,2:)/maxval(grad_sum(:,2:))

    call dat_write('tbtest.dat',dat_array,13)

    call dat_write('currenttest.dat',grad_sum,12)

    call dat_write('loctest.dat',loc_array,11)

    deallocate(dat_array, stat=istat)
    if(istat/=0) stop 'error deallocating dat_array array'

    deallocate(grad_array,stat=istat)
    if(istat/=0) stop 'error deallocating grad_array array'

    deallocate(grad_sum, stat=istat)
    if(istat/=0) stop 'error deallocating grad_sum array'

    deallocate(loc_array,stat=istat)
    if(istat/=0) stop 'error deallocating loc_array array'

    deallocate(loc_temp, stat=istat)
    if(istat/=0) stop 'error deallocating loc_temp array'

  else ! Default is to retain previous algorithm for single particle

    size=sys_size*cell%n_sites

    if(filling>size.or.filling==0) then 
      print*, 'invalid filling value, 1/2 filled current calculated'
      filling=int(real(size,kind=dp)/2.0_dp)
    end if
    
    k_val=phase*real_pi

    allocate(dat_array(1:phi_max+1,size+1),stat=istat)
    if(istat/=0) stop 'error allocating dat_array'

    allocate(loc_array(1:phi_max+1,size+1),stat=istat)
    if(istat/=0) stop 'error allocating loc_array array'

    allocate(loc_temp(size),stat=istat)
    if(istat/=0) stop 'error allocating loc_temp array'

    allocate(grad_array(1:phi_max+1,size+1),stat=istat)
    if(istat/=0) stop 'error allocating grad_array'

    grad_array=0.0_dp

    aexp=exp(cmplx_i*a_val)

    do i=1,phi_max+1

      ! System now iterates -pi -> pi
      phi=2.0_dp*real(i-1,kind=dp)/real(phi_max,kind=dp)-1.0_dp

      !defining phi as $\frac{\phi}{\phi_0}$
      theta_val=((real_pi*2)*phi)
      pexp=exp(cmplx_i*theta_val)

        kexp=aexp**(k_val)

        allocate(htn(size,size), stat=istat)
        if(istat/=0) stop 'error allocating htn matrix'

        htn=(0.0_dp,0.0_dp)

        call make_htn(sys_size, cell, kexp, pexp, htn)

        call zheev_evals(htn,egn)

        loc_array(i,1) =phi
        dat_array(i,1) =phi
        grad_array(i,1)=phi
        dat_array(i,2:)=egn(:)

        call get_local(htn,int(size,kind=i64),loc_temp)

        loc_array(i,2:)=loc_temp(:)

        deallocate(htn, stat=istat)
        if(istat/=0) stop 'error deallocating htn array'

        deallocate(egn, stat=istat)
        if(istat/=0) stop 'error deallocating egn array'
  
    end do

    d_phi = (2.0_dp * real_pi) / real(phi_max, kind=dp)

    call get_grad(dat_array(:, 2:size+1), grad_array(:, 2:size+1),&
                &d_phi, (phi_max+1), size)

    allocate(grad_sum(1:phi_max+1,2),stat=istat)
    if(istat/=0) stop 'error allocating grad_sum'

    fermi=maxval(dat_array(:,filling))

    grad_sum(:,1)=grad_array(:,1)
    do i=1,phi_max+1
      grad_sum(i,2) =sum(grad_array(i,2:filling+1))
    !  dat_array(i,2:)=dat_array(i,2:)-fermi
    end do

    grad_sum(:,2:)=grad_sum(:,2:)/maxval(grad_sum(:,2:))

    call dat_write('tbtest.dat',dat_array,13)

    call dat_write('currenttest.dat',grad_sum,12)

    call dat_write('loctest.dat',loc_array,11)

    deallocate(dat_array, stat=istat)
    if(istat/=0) stop 'error deallocating dat_array array'

    deallocate(grad_array,stat=istat)
    if(istat/=0) stop 'error deallocating grad_array array'

    deallocate(grad_sum, stat=istat)
    if(istat/=0) stop 'error deallocating grad_sum array'

    deallocate(loc_array, stat=istat)
    if(istat/=0) stop 'error deallocating loc_array array'

    deallocate(loc_temp, stat=istat)
    if(istat/=0) stop 'error deallocating loc_temp array'

  end if

end program main_project
