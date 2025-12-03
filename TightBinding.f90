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
subroutine make_htn(size, e_val, t_vals, kexp, pexp, t_table, htn)
  implicit none

  integer, intent(in)                            :: size
  real(kind=dp), intent(in)                      :: e_val
  real(kind=dp), dimension(:,:), intent(in)      :: t_vals
  complex(kind=dp), intent(in)                   :: kexp
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
      call intra_cell(N_site, t_table(:,:,1), t_vals, pexp, kexp,&
                   &intra_cell_row)
      call inter_cell(N_site, t_vals, t_table(:,:,2), kexp, pexp,&
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
      grad(i)=grad(i)+(3E8*(data(i)-data(i+1))/size(data))
    end do

    !grad(i+1)=0.0_dp

    return
  end subroutine get_grad


end module tight_binding

