module constants
  implicit none
  public

  integer, parameter:: dp=selected_real_kind(15,300)
  integer, parameter:: i32=selected_int_kind(32)

  complex(kind=dp), parameter:: cmplx_pi=(0.0_dp,3.1415926_dp)
  complex(kind=dp), parameter:: cmplx_i=(0.0_dp,1.0_dp)
  real(kind=dp), parameter :: real_pi=4.0_dp*atan(1.0_dp)

  contains
  subroutine dat_write(filename_in , data_in , file_unit)
    implicit none
    character(len = *) , intent(IN)                    :: filename_in
    real(kind = dp) , dimension(:,:) , intent(INOUT)   :: data_in
    integer , intent(IN)                               :: file_unit

    logical                                            :: lexist
    integer                                            :: ierr=0, i=0

    inquire(file=filename_in, exist=lexist)
    if (lexist) then                                       
      print*, "old data file:", filename_in," overwritten"     
      open (unit=file_unit, file=filename_in, status="replace"&
            &, action="write" ,position="append", iostat=ierr)
      if (ierr/= 0) stop "error overwriting data file"
    else
      open (unit=file_unit, file=filename_in, status="new"&
            &, action = "write", position="append",iostat = ierr)
      if (ierr/= 0) stop "error creating data file"
    end if

    do i=1,size(data_in , 1)
      write(unit=file_unit , fmt=* , iostat=ierr) data_in(i , :)
      if(ierr/=0) stop "error writing data array to file"
    end do
    close(unit=file_unit , iostat=ierr)
    if(ierr/=0) stop "error closing data file"
  end subroutine dat_write

  subroutine zheev_evals(A, W)
    implicit none
    character                                         :: JOBZ, UPLO
    integer                                           :: N, LDA, LWORK, INFO
    integer                                           :: IERR, IOSTAT = 0

    ! arrays and matrices
    complex(kind = dp) , dimension(:,:), intent(inout)   :: A
    real(kind = dp),dimension(:),allocatable, intent(out):: W
    complex(kind = dp) , dimension(:), allocatable       :: WORK
    real(kind = dp) , dimension(:), allocatable          :: RWORK

    JOBZ = 'N'
    UPLO = 'U'

    N = size(A, 1)
    LDA  = max(1, N)
    LWORK = max(1, 2*N - 1)
    INFO = IOSTAT

    allocate(W(N), stat=IERR)
    if (IERR /= 0) stop "failed to allocate W"

    allocate(WORK(max(1,LWORK)), stat=IERR)
    if (IERR /= 0) stop "failed to allocate WORK"

    allocate(RWORK(max(1, 3*N - 2)), stat=IERR)
    if (IERR /= 0) stop "failed to allocate RWORK"

    call zheev(JOBZ , UPLO , N , A , LDA , W, WORK, LWORK, RWORK, INFO)

    deallocate(WORK,stat=ierr)
    if(ierr/=0) stop 'error deallocating WORK'

    deallocate(RWORK,stat=ierr)
    if(ierr/=0) stop 'error deallocating RWORK'

    !print*, W(1)

  end subroutine zheev_evals

end module constants

module TB_Hamiltonian
  use constants
  implicit none
  public

  contains
!---------------------------------------------------------------------------!
!SUBROUTINE: on_site                                                        !
!IN-VALS   : N_site, epsilon_val                                            !
!INOUT-VALS  : htn_row                                                      !
!PURPOSE   : To allocate the onsite energy terms of the tight binding       !
!            Hamiltonian -i.e. the terms defined by <Rln|H|Rln> within the  !
!            handwritten notes. Subroutine is very simple as no 'neighbours'!
!            are required to be considered when allocating this component of!
!            the Hamiltonian matrix.                                        !
!DATE      : 07/10/2025                                                     !
!---------------------------------------------------------------------------!
  subroutine on_site(N_site, epsilon_val, htn_row)
    implicit none
    integer, intent(in)                          :: N_site
    real(kind=dp), intent(in)                    :: epsilon_val
    complex(kind=dp), dimension(:), intent(inout):: htn_row

    htn_row=0.0_dp

    htn_row(N_site)=epsilon_val
    return
  end subroutine on_site

!---------------------------------------------------------------------------!
!SUBROUTINE: intra_cell                                                     !
!IN-VALS   : N_site, t_vals, t_table, pexp                                  !
!INOUT-VALS: htn_row                                                        !
!PURPOSE   : To allocate the intra-cell hopping component of the Hamiltonian!
!            Given that this component of the Hamiltonian encompasses       !
!            hopping, the 'truth table' for the interations must be         !
!            considered in the form of the logical t_table 2 dimensional    !
!            array. This calculates the term defined by <Rlm|H|Rln> within  !
!            the tight binding formalism.                                   !
!DATE      : 07/10/2025, updated for fwrd hopping: 19/11/25                 !
!---------------------------------------------------------------------------!
  subroutine intra_cell(N_site, t_table, t_vals, pexp, htn_row)
    implicit none

    integer, intent(in)                          :: N_site
    logical, dimension(:,:), intent(in)          :: t_table
    real(kind=dp), dimension(:,:), intent(in)    :: t_vals
    complex(kind=dp), intent(in)                 :: pexp
    complex(kind=dp), dimension(:), intent(inout):: htn_row
    !complex(kind=dp), intent(in), optional       :: ctn

    integer:: i=0

    htn_row=0.0_dp

    do i=1,size(htn_row)
    !print*, 'allocating htn'
    !determining if the sites are 'hoppable', true means that they are not
    !as it is more effective to have as a default
      if(t_table(N_site,i).eqv..true.) cycle
      htn_row(i)=t_vals(1,1)*pexp
      if(N_site==1) htn_row(size(htn_row))  =htn_row(size(htn_row))*pexp**(-2)
      !Note removed -1.0_dp, replaced by making t=-1.0_dp default
      if(mod(i,N_site)==N_site-1) htn_row(i)=htn_row(i)*pexp**(-2)

      !print*, pexp
    end do
  end subroutine intra_cell

!---------------------------------------------------------------------------!
!SUBROUTINE: inter_cell                                                     !
!IN-VALS   : N_site, t_vals, t_table, texp, pexp                            !
!INOUT-VALS: htn_row                                                        !
!PURPOSE   : To allocate the inter-cell hopping component of the Hamiltonian!
!            Given that this component of the Hamiltonian encompasses       !
!            hopping, the 'truth table' for the interactions must be        !
!            considered in the form of the logical t_table 2 dimensional    !
!            array. The array used here is populated differently to the     !
!            truth table used for the components above. This component of   !
!            the Hamiltonian calculates the term given by <RLm|H|Rln> within!
!            the tight binding formalism.                                   !
!DATE      : 07/10/2025,edited for persistent current:06/11/25              !
!---------------------------------------------------------------------------!
  subroutine inter_cell(N_site, t_vals, t_table, texp, pexp, htn_row)
    implicit none

    integer, intent(in)                          :: N_site
    real(kind=dp), dimension(:,:), intent(in)    :: t_vals
    logical, dimension(:,:), intent(in)          :: t_table
    complex(kind=dp), intent(in)                 :: texp
    complex(kind=dp), intent(in)                 :: pexp
    complex(kind=dp), dimension(:), intent(inout):: htn_row
    !complex(kind=dp), optional, intent(in)       :: ctn

    integer:: i=0

    !pexp=pexp**N_site

    htn_row=0.0_dp

    !print*, 'N_site,inter_cell',N_site
    !print*, 'pexp,inter_cell', pexp

    do i=1,size(htn_row)
      if(t_table(N_site,i).eqv..true.) cycle

        !+1 factor added as nearest neighbour & different level
        !swaps dims of t_vals to correspond to correct format 17/11/25
        htn_row(i)=t_vals(2,abs(N_site-i+1))*texp*pexp&
                 &-t_vals(2,abs(N_site-i+1))*texp**(-1)*pexp**(-1)
        !print*, htn_row(i), 'inter_cell'
        !print*, pexp

        if(i>N_site) htn_row(i)=-1.0_dp*htn_row(i)

    end do

  end subroutine inter_cell
end module TB_Hamiltonian

module tight_binding
  use constants
  use TB_Hamiltonian
  implicit none
  public

  contains
!---------------------------------------------------------------------------!
!SUBROUTINE: make_t_table                                                   !
!IN-VALS   : nbr, size                                                      !
!OUT-VALS  : t_table                                                        !
!PURPOSE   : To calculate a table which describes the interations between   !
!            various sites within the 'unit cell' that is being considered. !
!            Given that no values are being stored, it is general to all k  !
!              a values etc. it also means that there is no need to         !
!            recalculate this matrix/table for each iteration of any given  !
!            run.                                                           !
!DATE      : 08/10/2025                                                     !
!---------------------------------------------------------------------------!
  subroutine make_t_table(size,t_table)
    implicit none

    integer, intent(in)                                :: size
    logical, dimension(:,:,:), allocatable, intent(out):: t_table


    integer :: istat=0,i=0,j=0,A,B

    allocate(t_table(size,size,2), stat=istat)
    if(istat/=0) stop 'error allocating t_table:tight_binding->make_t_table'

    !allocate so that there is no hopping by default
    t_table=.true.

    !initially considering the case of nearest neighbours only
    !considering intra-cell hopping
    !15/10/2025: Considering additional species of atoms within the system
    do i=1,size
      A=mod(i,size)
      do j=1,size
        B=mod(j,size)
        !or statement accounts for the periodic boundary conditions
        if(abs(A-B)==1) t_table(i,j,1)=.false.
        if(A==0) A=size
        if(B==0) B=size
        if(abs(A-B)==1) t_table(i,j,1)=.false.
      end do
    end do

    !considering inter-cell hopping
    !this works assuming that the unit cells are in layers
    do i=1,size
      do j=1,size
        if(abs((i+1)-j)==1) t_table(i,j,2)=.false.
      end do
    end do
    
  end subroutine make_t_table

!---------------------------------------------------------------------------!
!SUBROUTINE: make_htn                                                       !
!IN-VALS   : size, e_val, t_vals, texp, pexp                                !
!INOUT-VALS: htn                                                            !
!PURPOSE   : To combine the components of the Hamiltonian above and populate!
!            the full NxN matrix for the unit cell of size N.               !
!DATE      : 08/10/2025                                                     !
!---------------------------------------------------------------------------!
subroutine make_htn(size, e_val, t_vals, texp, pexp, t_table, htn)
  implicit none

  integer, intent(in)                            :: size
  real(kind=dp), intent(in)                      :: e_val
  real(kind=dp), dimension(:,:), intent(in)      :: t_vals
  complex(kind=dp), intent(in)                   :: texp
  complex(kind=dp), intent(inout)                :: pexp
  complex(kind=dp), intent(inout), dimension(:,:):: htn
  !neighbours is implicitly 1 here (08/10/2025)
  logical, dimension(:,:,:), intent(in)          :: t_table
  !complex(kind=dp), optional, intent(in)         :: ctn

  integer                                    :: istat=0, N_site
  complex(kind=dp), dimension(:), allocatable:: on_site_row
  complex(kind=dp), dimension(:), allocatable:: intra_cell_row
  complex(kind=dp), dimension(:), allocatable:: inter_cell_row
             

  allocate(on_site_row(size), stat=istat)
  if(istat/=0) stop 'error allocating on_site'

  allocate(intra_cell_row(size), stat=istat)
  if(istat/=0) stop 'error allocating intra_cell'

  allocate(inter_cell_row(size), stat=istat)
  if(istat/=0) stop 'error allocating inter_cell'

  do N_site=1,size


      !print*, 'pexp make_htn',pexp

      call on_site(N_site, e_val, on_site_row)
      call intra_cell(N_site, t_table(:,:,1), t_vals, pexp, intra_cell_row)
      call inter_cell(N_site, t_vals, t_table(:,:,2), texp, pexp,&
                   &inter_cell_row)


    htn(:,N_site)=on_site_row+intra_cell_row+inter_cell_row

  end do

  deallocate(on_site_row, stat=istat)
  if(istat/=0) stop 'error deallocating on_site'

  deallocate(intra_cell_row, stat=istat)
  if(istat/=0) stop 'error deallocating intra_cell'

  deallocate(inter_cell_row, stat=istat)
  if(istat/=0) stop 'error deallocating inter_cell'

  return

end subroutine make_htn

  subroutine get_grad(data,grad)
    implicit none

    real(kind=dp),dimension(:),intent(in)   :: data
    real(kind=dp),dimension(:),intent(inout):: grad

    integer:: i

    !grad=0.0_dp

    do i=1,size(data)-1
      grad(i)=grad(i)+(-1.0_dp*3E8*(data(i)-data(i+1))/size(data))
    end do

    !grad(i+1)=0.0_dp

    return
  end subroutine get_grad


end module tight_binding

program main_project
  !use omp_lib
  !use omp_lib_kinds
  use constants
  use tight_binding
  implicit none

  complex(kind=dp), dimension(:,:), allocatable:: htn
  integer                                      :: size,istat=0,i=1,nthreads=8
  real(kind=dp)                                :: e_val, a_val, k_val,&
                                                 &theta_val, phi=0.0_dp
  real(kind=dp), dimension(:), allocatable     :: egn
  real(kind=dp), dimension(:,:),allocatable    :: t_vals,dat_array,grad_array
  complex(kind=dp)                             :: kexp, aexp, pexp!,ctn
  logical, dimension(:,:,:), allocatable       :: t_table
  integer                                      :: phi_max=1E4,k
  !character(len=100)                             :: solve
  !logical                                      :: dynamic
  !Note- as above, neighbours is set to 1 for initial testing

  !solve='current'
  size=20
  e_val=0.0_dp
  !a_val should be high for effective single layer ring
  a_val=1.0_dp
  k_val=0.0_dp

  !Note neighbours again
  allocate(t_vals(2,1), stat=istat)
  if(istat/=0) stop 'error allocating t_vals array'

  allocate(dat_array(-phi_max:phi_max,size+1),stat=istat)
  if(istat/=0) stop 'error allocating dat_array'

  allocate(grad_array(-phi_max:phi_max,2),stat=istat)
  if(istat/=0) stop 'error allocating grad_array'

  grad_array=0.0_dp

  t_vals(1,1)=-1.0_dp
  t_vals(2,1)=0.0_dp !only considering one layer of ring

  !calculating texp
  aexp=exp(cmplx_i*a_val)

  !phi=real_pi/4.0_dp

  call make_t_table(size,t_table)

  !do i=1,size
  !  print*, t_table(i,:,1)
  !end do

  !$OMP parallel do default(none) &
  !$OMP & private(i,phi,theta_val,pexp,htn,egn,k_val,kexp) &
  !$OMP & shared(phi_max,size,aexp,istat,e_val,t_vals,a_val,t_table) &
  !$OMP & schedule(dynamic) num_threads(nthreads) reduction(+:dat_array)

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

  !$OMP end parallel do

  do k=2,int(size/2)
    call get_grad(dat_array(:,k),grad_array(:,2))
    !print*, int(size/2)
  end do

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
