program main_project
  !use omp_lib
  !use omp_lib_kinds
  use constants
  use tight_binding
  implicit none

  complex(kind=dp), dimension(:,:), allocatable:: htn
  integer                                      :: size,istat=0,i=1,&
                                                 &filling
  real(kind=dp)                                :: e_val, a_val, k_val,phase,&
                                                 &theta_val, phi=0.0_dp,t_val
  real(kind=dp), dimension(:), allocatable     :: egn
  real(kind=dp), dimension(:,:),allocatable    :: t_vals,dat_array,grad_array
  complex(kind=dp)                             :: kexp, aexp, pexp!,ctn
  real(kind=dp)                                :: Inought
  logical, dimension(:,:,:), allocatable       :: t_table
  integer                                      :: phi_max=1E4,k
  !character(len=100)                             :: solve
  !logical                                      :: dynamic
  !Note- as above, neighbours is set to 1 for initial testing

  call get_params(size,filling,phase,a_val,e_val,t_val)

  if(filling>size.or.filling==0) then 
    print*, 'invalid filling value, 1/2 filled current calculated'
    filling=int(real(size,kind=dp)/2.0_dp)
  end if
  
  k_val=phase*real_pi

  !Note neighbours again
  allocate(t_vals(2,1), stat=istat)
  if(istat/=0) stop 'error allocating t_vals array'

  allocate(dat_array(-phi_max:phi_max,size+1),stat=istat)
  if(istat/=0) stop 'error allocating dat_array'

  allocate(grad_array(-phi_max:phi_max,2),stat=istat)
  if(istat/=0) stop 'error allocating grad_array'

  grad_array=0.0_dp

  t_vals(1,1)=-1.0_dp*t_val
  t_vals(2,1)=0.0_dp !only considering one layer of ring

  !calculating texp
  aexp=exp(cmplx_i*a_val)
  !phi=0.0_dp!real_pi/4.0_dp

  call make_t_table(size,t_table)

  !do i=1,size
  !  print*, t_table(i,:,1)
  !end do

  do i=-phi_max,phi_max

    phi=real(i,kind=dp)*(real_pi/2.0_dp)/real(phi_max,kind=dp)

    !print*, phi
 
    !defining phi as $\frac{\phi}{\phi_0}$
    theta_val=((real_pi*2)*phi)/real(size, kind=dp)
    !Note: fudge factor
    !theta_val=theta_val*4.0_dp/3.0_dp
    pexp=exp(cmplx_i*theta_val)

    !print*, 'pexp main',pexp

    !Note- neighbours still assumed to be 1
    kexp=aexp**(k_val)
    !print*, kexp

    allocate(htn(size,size), stat=istat)
    if(istat/=0) stop 'error allocating htn matrix'

    call make_htn(size,e_val,t_vals,kexp,pexp,t_table,htn)
    !do k=1,size
      !print*, htn(1,1)
      !print*, htn(size,size) 
    !end do
    !print*, t_table(1,8,1)
    !do i=1,size
    !  print*, t_table(:,i,1)
    !end do

    call zheev_evals(htn,egn)

    !print*, phi

    dat_array(i,1) =phi
    grad_array(i,1)=phi
    dat_array(i,2:)=egn(:)

    !print*, dat_array(i,:2)

    !print*, i

    deallocate(htn, stat=istat)
    if(istat/=0) stop 'error deallocating htn array'

    deallocate(egn, stat=istat)
    if(istat/=0) stop 'error deallocating egn array'
 
  end do

  do k=2,filling+1
    call get_grad(dat_array(:,k),grad_array(:,2))
  end do
  Inought=maxval(grad_array(:,2))
  grad_array(:,2)=grad_array(:,2)/Inought

  call dat_write('tbtest.dat',dat_array,13)

  call dat_write('currenttest.dat',grad_array,12)
  
  deallocate(t_vals, stat=istat)
  if(istat/=0) stop 'error deallocating t_vals array'

  deallocate(dat_array, stat=istat)
  if(istat/=0) stop 'error deallocating dat_array array'

  deallocate(grad_array,stat=istat)
  if(istat/=0) stop 'error deallocating grad_array array'

  deallocate(t_table, stat=istat)
  if(istat/=0) stop 'error deallocating t_table'

end program main_project
