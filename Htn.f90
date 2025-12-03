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
  subroutine intra_cell(N_site, t_table, t_vals, pexp, kexp, htn_row)
    implicit none

    integer, intent(in)                          :: N_site
    logical, dimension(:,:), intent(in)          :: t_table
    real(kind=dp), dimension(:,:), intent(in)    :: t_vals
    complex(kind=dp), intent(in)                 :: pexp,kexp
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

      !allocating hopping internal to ring    :- phi dept
      !Note removed -1.0_dp, replaced by making t=-1.0_dp default
      if(mod(i,N_site)==N_site-1) htn_row(i)=htn_row(i)*pexp**(-2)

      !print*, pexp
    end do
    !completing edge case
    if(N_site==1) htn_row(size(htn_row))=htn_row(size(htn_row))*pexp**(-2)

    !allocating hopping btwn revs of ring
    if(N_site==1) htn_row(size(htn_row))=htn_row(size(htn_row))*kexp**(-1)
    if(N_site==size(htn_row)) htn_row(1)=htn_row(1)*kexp
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

