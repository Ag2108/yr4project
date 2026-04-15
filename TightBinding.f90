module tight_binding
  use constants
  !use TB_Hamiltonian
  implicit none
  public

  type:: unit_cell
    integer                                    :: n_sites!gives number of sites in cell
    integer                                    :: nbr_range
    real(kind=dp),dimension(:),allocatable     :: epsilon
    complex(kind=dp),dimension(:,:),allocatable:: intra
    complex(kind=dp),dimension(:,:),allocatable:: inter
  end type unit_cell

  contains

!---------------------------------------------------------------------------!
!SUBROUTINE: init_sys                                                      !
!IN-VALS   : NONE                                                           !
!OUT_VALS  : CELL; contains required information for the system Hamiltonian !
!            SYS_SIZE; the number of sites within the system globally       !
!            FILLING; used to calculate how many bands are used to calculate!
!            the current; see get_grad()                                    !
!PURPOSE   : Acts as a wrapper subroutine for get_params so that less juggl-!
!            -ing is required.                                              !
!DATE      : REPLACED make_t_table 09/03/2026                               !
!---------------------------------------------------------------------------!

  subroutine init_sys(cell,sys_size,filling,phase_case,a_val,phi_max,&
                     &many_body,U_val,k_max)
    use constants
    implicit none

    type(unit_cell),intent(out)     :: cell
    integer,intent(out)             :: sys_size ! Number of cells
    integer,intent(out)             :: filling
    integer,intent(out)             :: phi_max
    integer,intent(out)             :: k_max
    real(kind=dp),intent(OUT)       :: phase_case
    real(kind=dp),intent(OUT)       :: a_val
    integer,intent(out),dimension(2):: many_body
    real(kind=dp),intent(out)       :: U_val

    integer                              :: cell_size,nbr_range
    logical,dimension(:,:),allocatable   :: intra_tbl,inter_tbl

    integer                                 :: istat=0,n,i,j
    logical,dimension(:,:),allocatable      :: tmp_map
    integer,dimension(:,:),allocatable      :: dist_inter,dist_intra,tmp_dist
    real(kind=dp),dimension(:,:),allocatable:: sites

    call get_params(sys_size,sites,filling,phase_case,a_val,&
                   &cell_size,intra_tbl,inter_tbl,nbr_range,&
                   &phi_max,many_body,U_val,k_max)

    cell%n_sites=cell_size
    cell%nbr_range=nbr_range
    allocate(cell%epsilon(cell_size),&
            &cell%intra(cell_size,cell_size),&
            &cell%inter(cell_size,cell_size),&
            &stat=istat)
    if(istat/=0)stop 'Error in allocation of cell arrays'
    
    cell%epsilon(:)=sites(:,1)

    ! Allocate a temporary map for 3*cells
    ! Assume that no more than neighbour cells are required; higher hops can be put to 0
    allocate(tmp_map(3*cell_size, 3*cell_size),tmp_dist(3*cell_size, 3*cell_size),stat=istat)
    if(istat/=0) stop 'error allocating tmp_map,tmp_dist'
    tmp_map=.true. ! Assume that no sites can initially see each other (.false.->hopping)

    ! 3 Intra-cell blocks down the diagonal
    tmp_map(1:cell_size, 1:cell_size)                             = intra_tbl  ! Left Cell
    tmp_map(cell_size+1:2*cell_size, cell_size+1:2*cell_size)     = intra_tbl ! Centre Cell
    tmp_map(2*cell_size+1:3*cell_size, 2*cell_size+1:3*cell_size) = intra_tbl ! Right Cell

    ! Inter-cell forward hops (Left -> Centre, Centre -> Right)
    tmp_map(1:cell_size, cell_size+1:2*cell_size)                 = inter_tbl
    tmp_map(cell_size+1:2*cell_size, 2*cell_size+1:3*cell_size)   = inter_tbl

    ! Inter-cell backward hops (Centre -> Left, Right -> Centre)
    ! Transpose implemented to ensure correct ordering
    tmp_map(cell_size+1:2*cell_size, 1:cell_size)                 = transpose(inter_tbl)
    tmp_map(2*cell_size+1:3*cell_size, cell_size+1:2*cell_size)   = transpose(inter_tbl)

    call get_dist(3*cell_size,tmp_map,tmp_dist)

    allocate(dist_inter(cell_size,cell_size),&
            &dist_intra(cell_size,cell_size),stat=istat)
    if(istat/=0) stop 'error with nbr hopping allocation'

    dist_intra=tmp_dist(cell_size+1:2*cell_size,cell_size+1:2*cell_size)
    dist_inter=tmp_dist(cell_size+1:2*cell_size,2*cell_size+1:3*cell_size)

    deallocate(tmp_dist,stat=istat)
    if(istat/=0) stop 'error deallocating tmp_dist'
    
    ! Initialisation of intra hopping tables
    cell%intra(:,:)=0.0_dp
    cell%inter(:,:)=0.0_dp
    ! Nbrs defined from dist
    do n=1,nbr_range
      do j=1,cell_size
        do i=1,cell_size
          if(dist_intra(i,j)==n) then
            cell%intra(i,j)=cell%intra(i,j)+sites(i,n+1)
          end if
          if(dist_inter(i,j)==n) then
            cell%inter(i,j)=cell%inter(i,j)+sites(i,cell%nbr_range+n+1)
            ! Given that only fwrd hopping is included in inter 
            !/sf(both are included in intra; see get_params)
            !cell%inter(j,i)=cell%inter(j,i)+sites(i,n+1)
          end if
        end do
      end do
    end do
    ! Note that the pexp term has not yet been added
    deallocate(sites,stat=istat)
    if(istat/=0) stop 'error deallocating sites'

  end subroutine init_sys

!---------------------------------------------------------------------------!
!SUBROUTINE: make_htn                                                       !
!IN-VALS   : size, e_val, t_vals, texp, pexp                                !
!INOUT-VALS: htn                                                            !
!PURPOSE   : To combine the components of the Hamiltonian above and populate!
!            the full NxN matrix for the unit cell of size N.               !
!DATE      : 08/10/2025; REWORKED for CELL 10/03/2026                       !
!---------------------------------------------------------------------------!
subroutine make_htn(sys_size, cell, kexp, pexp, htn)
  implicit none

  type(unit_cell),intent(in)                     :: cell
  integer, intent(in)                            :: sys_size 
  ! Gives the NUMBER OF CELLS
  complex(kind=dp), intent(in)                   :: kexp
  complex(kind=dp), intent(inout)                :: pexp
  complex(kind=dp), intent(inout), dimension(:,:):: htn

  integer:: cc_start,cc_end ! Integers used to define parts of htn which 
  integer:: nc_start,nc_end ! Integers used to define PBCs, inter hopping
  integer:: on_site         ! Integer to define placement of epsilon
  integer:: i,j

  do i=1,sys_size

    cc_start=cell%n_sites*(i-1)+1 ! Indexed from 1
    cc_end  =cell%n_sites*(i)     ! Indexed from number of sites in cell

    ! Intra cell hopping (defined in init_sys()) applied to htn
    htn(cc_start:cc_end,cc_start:cc_end)=cell%intra 
    ! Given that each cell is homogenous

    do j=1,cell%n_sites
      on_site=cc_start+j-1
      htn(on_site,on_site)=htn(on_site,on_site)+cell%epsilon(j)
    end do

    ! Inter cell hopping (defined in init_sys()) applied to htn
    if(i<sys_size) then 
      ! Always want the 'next cell'; back hopping included through transpose
      nc_start=cell%n_sites*i+1
      nc_end  =cell%n_sites*(i+1)

      ! Forward hop defined here
      htn(cc_start:cc_end,nc_start:nc_end)=&
      &htn(cc_start:cc_end,nc_start:nc_end)+cell%inter

      ! Backwards hop defined here
      htn(nc_start:nc_end,cc_start:cc_end)=&
      &htn(nc_start:nc_end,cc_start:cc_end)+transpose(cell%inter)
    end if
  end do

  ! PBCs defined here
  cc_start=(sys_size-1)*cell%n_sites+1 ! start of the last cell
  cc_end  =sys_size*cell%n_sites
  nc_start=1
  nc_end  =cell%n_sites

  ! Forward boundary hop
  htn(cc_start:cc_end,nc_start:nc_end)=&
  &htn(cc_start:cc_end,nc_start:nc_end)+cell%inter*pexp*kexp

  ! Reverse boundary hop
  htn(nc_start:nc_end,cc_start:cc_end)=&
  &htn(nc_start:nc_end,cc_start:cc_end)+transpose(cell%inter)*conjg(pexp*kexp)

end subroutine make_htn

  subroutine get_grad(data,grad,d_dat,points,bands)
    implicit none

    real(kind=dp),dimension(:,:),intent(in) :: data
    real(kind=dp),dimension(:,:),intent(out):: grad
    real(kind=dp),intent(in)                :: d_dat
    integer,intent(in)                      :: points,bands

    integer:: p,b

    grad=0.0_dp

    do b=1,bands

      grad(1,b)=(data(1,b)-data(2,b))/d_dat

      do p=2,points-1
        grad(p,b)=((data(p-1,b)-data(p+1,b)))/(2.0_dp*d_dat)
      end do

      grad(points,b)=(data(p-1,b)-data(p,b))/d_dat

    end do

    return
  end subroutine get_grad

    subroutine get_basis(n_up,n_down,n_sites,hbt_tot,d_tot)

    integer,intent(in)                                      :: n_up,n_down,n_sites
    integer(kind=i64),dimension(:,:),allocatable,intent(out):: hbt_tot
    
    
    integer(kind=i64),dimension(:),allocatable:: hbt_up,hbt_down
    integer                                   :: state_up,state_down,state_tot,i,j
    integer(kind=i64)                         :: num_states,ics
    integer(kind=i64)                         :: d_up,d_down,d_tot

    ! Defining the number of dims in UP and DOWN Hbt spaces
    d_up  =nCr_choose(n_sites,n_up)
    d_down=nCr_choose(n_sites,n_down)
    d_tot =d_up*d_down

    ! Allocating the arrays that hold the states for the Hbt spaces
    allocate(hbt_up(d_up))
    allocate(hbt_down(d_down))
    allocate(hbt_tot(d_tot,2))

    ! The MAXIMUM number of states is given as 2^L per spin state
    num_states=(2_i64**n_sites)-1

    ! Finding the integers which correspond to system filling
    state_up  =1
    state_down=1
    do ics=0_i64,num_states
      if(popcnt(ics)==n_up)then
        hbt_up(state_up)=ics
        state_up=state_up+1
      end if

      if(popcnt(ics)==n_down)then
        hbt_down(state_down)=ics
        state_down=state_down+1
      end if
    end do

    state_tot=1
    do i=1,int(d_up)
      do j=1,int(d_down)
        hbt_tot(state_tot,1)=hbt_up(i)
        hbt_tot(state_tot,2)=hbt_down(j)
        state_tot=state_tot+1
      end do
    end do

    !print*, 'states generated', d_tot

  end subroutine get_basis

  subroutine make_hmb(sys_size, cell, kexp, pexp, U_val, d_tot, hbt_tot, hmb)
    implicit none

    type(unit_cell),intent(in)                         :: cell
    integer, intent(in)                                :: sys_size 
    ! Gives the NUMBER OF CELLS
    complex(kind=dp), intent(in)                       :: kexp
    complex(kind=dp), intent(inout)                    :: pexp
    real(kind=dp),intent(in)                           :: U_val
    integer(kind=i64),intent(in)                       :: d_tot
    integer(kind=i64),dimension(d_tot,2)               :: hbt_tot
    complex(kind=dp),dimension(d_tot,d_tot),intent(out):: hmb

    integer(kind=i64):: state_up,state_down,new_up,new_down
    integer          :: col,row,i,j,k,updown
    integer          :: cc_start,nc_start,on_site    

    hmb=0.0_dp

    do col=1,int(d_tot)
      state_up  =hbt_tot(col,1)
      state_down=hbt_tot(col,2)

      ! Calculate \sum_nn_\uparrow n_\downarrow
      updown=popcnt(iand(state_up,state_down))
      hmb(col, col)=hmb(col, col)+U_val*real(updown,kind=dp)

      do i=1,sys_size

        cc_start=cell%n_sites*(i-1) ! Indexed from 0

        do j=1,cell%n_sites
          on_site=cc_start+j
          ! On-site energy for spin up state
          if(btest(state_up,on_site-1)) hmb(col,col)=hmb(col,col)+cell%epsilon(j)
          ! On-site energy for spin down state
          if(btest(state_down,on_site-1)) hmb(col,col)=hmb(col,col)+cell%epsilon(j)

          do k=1,cell%n_sites

            ! Intra-cell hopping (no AB phase applied)
            ! Test if the hop is possible based off 'range'
            if(abs(cell%intra(j,k))<1E-6.and.abs(cell%inter(j,k))<1E-6) cycle

            ! Test if the hop is legal based off Pauli
            ! UP
            if((btest(state_up, cc_start+j-1).and..not.btest(state_up, cc_start+k-1))&
              &.and.abs(cell%intra(j,k))>1E-6) then
              new_up=ibclr(state_up,cc_start+j-1)
              new_up=ibset(new_up,cc_start+k-1)
              !print*,'made it to A'

              ! Find the matrix state which corresponds to new_up
              row=get_index(d_tot,hbt_tot,new_up,state_down)
              hmb(row,col)=hmb(row,col)+fermi_sign(cc_start+j,cc_start+k,state_up)*cell%intra(j,k)

            end if

            ! DOWN
            if((btest(state_down, cc_start+j-1).and..not.btest(state_down, cc_start+k-1))&
              &.and.abs(cell%intra(j,k))>1E-6) then
              new_down=ibclr(state_down,cc_start+j-1)
              new_down=ibset(new_down,cc_start+k-1)
              !print*,'made it to B'

              ! Find the matrix state which corresponds to new_down
              row=get_index(d_tot,hbt_tot,state_up,new_down)
              hmb(row,col)=hmb(row,col)+fermi_sign(cc_start+j,cc_start+k,state_down)*cell%intra(j,k)

              !hmb(col,row)=hmb(col,row)+fermi_sign*cell%intra(k,j)
            end if

            ! Inter cell hopping defined here
            if(i<sys_size) then 
              nc_start=cell%n_sites*i

              ! UP
              if((btest(state_up, cc_start+j-1).and..not.btest(state_up, nc_start+k-1))&
                &.and.abs(cell%inter(j,k))>1E-6) then
                new_up=ibclr(state_up,cc_start+j-1)
                new_up=ibset(new_up,nc_start+k-1)
                !print*,'made it to C'

                ! Find the matrix state which corresponds to new_up
                row=get_index(d_tot,hbt_tot,new_up,state_down)
                hmb(row,col)=hmb(row,col)+fermi_sign(cc_start+j,nc_start+k,state_up)*cell%inter(j,k)
                ! Taking the transpose for reverse hopping
                ! Much like in htn, only fwrd hopping is defined by default
                hmb(col,row)=hmb(col,row)+fermi_sign(cc_start+j,nc_start+k,state_up)*cell%inter(j,k)
              end if

              ! DOWN
              if((btest(state_down, cc_start+j-1).and..not.btest(state_down, nc_start+k-1))&
                &.and.abs(cell%inter(j,k))>1E-6) then
                new_down=ibclr(state_down,cc_start+j-1)
                new_down=ibset(new_down,nc_start+k-1)
                !print*,'made it to D'

                !print*, state_up, new_down, cc_start+j-1,nc_start+k-1

                ! Find the matrix state which corresponds to new_down
                row=get_index(d_tot,hbt_tot,state_up,new_down)
                hmb(row,col)=hmb(row,col)+fermi_sign(cc_start+j,nc_start+k,state_down)*cell%inter(j,k)
                ! Taking the transpose for reverse hopping
                ! Much like in htn, only fwrd hopping is defined by default
                hmb(col,row)=hmb(col,row)+fermi_sign(cc_start+j,nc_start+k,state_down)*cell%inter(j,k)
              end if

            ! PBCs are defined here; analogous to in single particle system
            else

              cc_start=(sys_size-1)*cell%n_sites ! start of the last cell
              nc_start=0

              ! UP
              if((btest(state_up, cc_start+j-1).and..not.btest(state_up, nc_start+k-1))&
                &.and.abs(cell%inter(j,k))>1E-6) then
                new_up=ibclr(state_up,cc_start+j-1)
                new_up=ibset(new_up,nc_start+k-1)
                !print*,'made it to E'

                ! Find the matrix state which corresponds to new_up
                row=get_index(d_tot,hbt_tot,new_up,state_down)
                hmb(row,col)=hmb(row,col)+fermi_sign(cc_start+j,nc_start+k,state_up)*cell%inter(j,k)*pexp*kexp
                ! Taking the transpose for reverse hopping
                ! Much like in htn, only fwrd hopping is defined by default
                hmb(col,row)=hmb(col,row)+fermi_sign(cc_start+j,nc_start+k,state_up)*cell%inter(j,k)*conjg(pexp*kexp)
              end if

              ! DOWN
              if((btest(state_down, cc_start+j-1).and..not.btest(state_down, nc_start+k-1))&
                &.and.abs(cell%inter(j,k))>1E-6) then
                new_down=ibclr(state_down,cc_start+j-1)
                new_down=ibset(new_down,nc_start+k-1)
                !print*,'made it to F'

                ! Find the matrix state which corresponds to new_down
                row=get_index(d_tot,hbt_tot,state_up,new_down)
                hmb(row,col)=hmb(row,col)+fermi_sign(cc_start+j,nc_start+k,state_down)*cell%inter(j,k)*pexp*kexp
                ! Taking the transpose for reverse hopping
                ! Much like in htn, only fwrd hopping is defined by default
                hmb(col,row)=hmb(col,row)+fermi_sign(cc_start+j,nc_start+k,state_down)*cell%inter(j,k)*conjg(pexp*kexp)
              end if
            end if
          end do
        end do
      end do
    end do

  end subroutine make_hmb

  integer function get_index(d_tot,hbt_tot,target_up,target_down)
    implicit none

    integer(kind=i64),intent(in)                   :: d_tot
    integer(kind=i64),dimension(d_tot,2),intent(in):: hbt_tot
    integer(kind=i64),intent(in)                   :: target_up,target_down

    integer:: i,j

    get_index=0

    do i=1,int(d_tot)
      if(hbt_tot(i,1)==target_up)then
        do j=i,int(d_tot)
          if(hbt_tot(j,2)==target_down)then
            if(hbt_tot(j,1)/=target_up) exit
            get_index=j
            return
          end if
        end do
      end if
    end do

    if(get_index==0) stop 'Error finding index in get_index; see Hubbard'
  end function get_index

  real(kind=dp) function fermi_sign(j_index,k_index,state_in)
    implicit none
    integer,intent(in)          :: j_index,k_index
    integer(kind=i64),intent(in):: state_in

    integer          :: min_state,max_state
    integer(kind=i64):: check_sites

    min_state=min(j_index-1,k_index-1)
    max_state=max(j_index-1,k_index-1)
    check_sites=(ishft(1_i64,int(max_state,kind=i64)))-&
               &(ishft(1_i64,int(min_state,kind=i64)+1_i64))

    ! Use the bin parity poppar() to return 1,-1 as required
    fermi_sign=1.0_dp-2.0_dp*real(poppar(iand(state_in,check_sites)),kind=dp)
    return
  end function fermi_sign

end module tight_binding

