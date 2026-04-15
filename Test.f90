program test
  use constants
  use tight_binding
  implicit none

  integer(kind=i64),dimension(:,:),allocatable:: hbt
  integer(kind=i64)                           :: tot

  call get_basis(4,4,8,hbt,tot)

!  contains
!   subroutine get_basis(n_up,n_down,n_sites,hbt_tot,d_tot)

!     integer,intent(in)                                      :: n_up,n_down,n_sites
!     integer(kind=i64),dimension(:,:),allocatable,intent(out):: hbt_tot
    
    
!     integer(kind=i64),dimension(:),allocatable:: hbt_up,hbt_down
!     integer                                   :: state_up,state_down,state_tot,i,j
!     integer(kind=i64)                         :: num_states,ics
!     integer(kind=i64)                         :: d_up,d_down,d_tot

!     ! Defining the number of dims in UP and DOWN Hbt spaces
!     d_up  =nCr_choose(n_sites,n_up)
!     d_down=nCr_choose(n_sites,n_down)
!     d_tot =d_up*d_down

!     ! Allocating the arrays that hold the states for the Hbt spaces
!     allocate(hbt_up(d_up))
!     allocate(hbt_down(d_down))
!     allocate(hbt_tot(d_tot,2))

!     ! The MAXIMUM number of states is given as 2^L per spin state
!     num_states=(2_i64**n_sites)-1

!     ! Finding the integers which correspond to system filling
!     state_up  =1
!     state_down=1
!     do ics=0_i64,num_states
!       if(popcnt(ics)==n_up)then
!         hbt_up(state_up)=ics
!         state_up=state_up+1
!       end if

!       if(popcnt(ics)==n_down)then
!         hbt_down(state_down)=ics
!         state_down=state_down+1
!       end if
!     end do

!     state_tot=1
!     do i=1,int(d_up)
!       do j=1,int(d_down)
!         hbt_tot(state_tot,1)=hbt_up(i)
!         hbt_tot(state_tot,2)=hbt_down(j)
!         state_tot=state_tot+1
!       end do
!     end do

!   end subroutine get_basis

!   subroutine make_hmb(sys_size, cell, kexp, pexp, U_val, d_tot, hbt_tot, hmb)
!     implicit none

!     type(unit_cell),intent(in)                         :: cell
!     integer, intent(in)                                :: sys_size 
!     ! Gives the NUMBER OF CELLS
!     complex(kind=dp), intent(in)                       :: kexp
!     complex(kind=dp), intent(inout)                    :: pexp
!     real(kind=dp),intent(in)                           :: U_val
!     integer(kind=i64),intent(in)                       :: d_tot
!     integer(kind=i64),dimension(d_tot,2)               :: hbt_tot
!     complex(kind=dp),dimension(d_tot,d_tot),intent(out):: hmb

!     real(kind=dp)    :: fermi_sign
!     integer(kind=i64):: state_up,state_down,new_up,new_down,check_sites
!     integer          :: col,row,i,j,k,updown
!     integer          :: cc_start,nc_start,on_site,min_state,max_state     

!     hmb=0.0_dp

!     do col=1,int(d_tot)
!       state_up  =hbt_tot(col,1)
!       state_down=hbt_tot(col,2)

!       updown=popcnt(iand(state_up,state_down))
!       hmb(col, col)=hmb(col, col)+U_val*real(updown,kind=dp)

!       do i=1,sys_size

!         cc_start=cell%n_sites*(i-1)+1 ! Indexed from 1

!         do j=1,cell%n_sites
!           on_site=cc_start+j-1
!           ! On-site energy for spin up state
!           if(btest(state_up,on_site-1)) hmb(col,col)=hmb(col,col)+cell%epsilon(j)
!           ! On-site energy for spin down state
!           if(btest(state_down,on_site-1)) hmb(col,col)=hmb(col,col)+cell%epsilon(j)

!           do k=1,cell%n_sites

!             ! Intra-cell hopping (no AB phase applied)
!             ! Test if the hop is possible based off 'range'
!             if(abs(cell%intra(j,k))<1E-6.or.abs(cell%inter(j,k))<1E-6) cycle

!             ! Test if the hop is legal based off Pauli
!             if((btest(state_up, cc_start+j-1).and..not.btest(state_up, cc_start+k-1))&
!               &.and.abs(cell%intra(j,k))>1E-6) then
!               new_up=ibclr(state_up,cc_start+j-1)
!               new_up=ibset(new_up,cc_start+k-1)

!               min_state=min(cc_start+j-1,cc_start+k-1)
!               max_state=max(cc_start+j-1,cc_start+k-1)
!               check_sites=(2_i64**(max_state) - 1_i64) - (2_i64**(min_state + 1) - 1_i64)

!               ! Use the bin parity poppar() to return 1,-1 as required
!               fermi_sign=1.0_dp-2.0_dp*real(poppar(iand(state_up,check_sites)),kind=dp)

!               ! Find the matrix state which corresponds to new_up
!               row=get_index(d_tot,hbt_tot,new_up,state_down)
!               hmb(row,col)=hmb(row,col)+fermi_sign*cell%intra(j,k)
!               ! Taking the transpose for reverse hopping
!               ! Much like in htn, only fwrd hopping is defined by default
!               hmb(col,row)=hmb(col,row)+fermi_sign*cell%intra(k,j)
!             end if

!             if((btest(state_down, cc_start+j-1).and..not.btest(state_down, cc_start+k-1))&
!               &.and.abs(cell%intra(j,k))>1E-6) then
!               new_down=ibclr(state_down,cc_start+j-1)
!               new_down=ibset(new_down,cc_start+k-1)

!               min_state=min(cc_start+j-1,cc_start+k-1)
!               max_state=max(cc_start+j-1,cc_start+k-1)
!               check_sites=(2_i64**(max_state) - 1_i64) - (2_i64**(min_state + 1) - 1_i64)

!               ! Use the bin parity poppar() to return 1,-1 as required
!               ! Combination of poppar, iand represents finding the parity of the shared sites within states
!               fermi_sign=1.0_dp-2.0_dp*real(poppar(iand(state_down,check_sites)),kind=dp)

!               ! Find the matrix state which corresponds to new_down
!               row=get_index(d_tot,hbt_tot,state_up,new_down)
!               hmb(row,col)=hmb(row,col)+fermi_sign*cell%intra(j,k)

!               hmb(col,row)=hmb(col,row)+fermi_sign*cell%intra(k,j)
!             end if

!             ! Inter cell hopping defined here
!             if(i<sys_size) then 
!               nc_start=cell%n_sites*i+1

!               if((btest(state_up, cc_start+j-1).and..not.btest(state_up, nc_start+k-1))&
!                 &.and.abs(cell%inter(j,k))>1E-6) then
!                 new_up=ibclr(state_up,cc_start+j-1)
!                 new_up=ibset(new_up,nc_start+k-1)

!                 min_state=min(cc_start+j-1,nc_start+k-1)
!                 max_state=max(cc_start+j-1,nc_start+k-1)
!                 check_sites=(2_i64**(max_state) - 1_i64) - (2_i64**(min_state + 1) - 1_i64)

!                 ! Use the bin parity poppar() to return 1,-1 as required
!                 fermi_sign=1.0_dp-2.0_dp*real(poppar(iand(state_up,check_sites)),kind=dp)

!                 ! Find the matrix state which corresponds to new_up
!                 row=get_index(d_tot,hbt_tot,new_up,state_down)
!                 hmb(row,col)=hmb(row,col)+fermi_sign*cell%inter(j,k)
!                 ! Taking the transpose for reverse hopping
!                 ! Much like in htn, only fwrd hopping is defined by default
!                 hmb(col,row)=hmb(col,row)+fermi_sign*cell%inter(k,j)
!               end if

!               if((btest(state_up, cc_start+j-1).and..not.btest(state_up, cc_start+k-1))&
!                 &.and.abs(cell%inter(j,k))>1E-6) then
!                 new_down=ibclr(state_down,cc_start+j-1)
!                 new_down=ibset(new_down,nc_start+k-1)

!                 min_state=min(cc_start+j-1,nc_start+k-1)
!                 max_state=max(cc_start+j-1,nc_start+k-1)
!                 check_sites=(2_i64**(max_state) - 1_i64) - (2_i64**(min_state + 1) - 1_i64)

!                 ! Use the bin parity poppar() to return 1,-1 as required
!                 fermi_sign=1.0_dp-2.0_dp*real(poppar(iand(state_up,check_sites)),kind=dp)

!                 ! Find the matrix state which corresponds to new_down
!                 row=get_index(d_tot,hbt_tot,state_up,new_down)
!                 hmb(row,col)=hmb(row,col)+fermi_sign*cell%inter(j,k)
!                 ! Taking the transpose for reverse hopping
!                 ! Much like in htn, only fwrd hopping is defined by default
!                 hmb(col,row)=hmb(col,row)+fermi_sign*cell%inter(k,j)
!               end if

!             ! PBCs are defined here; analogous to in single particle system
!             else

!               cc_start=(sys_size-1)*cell%n_sites+1 ! start of the last cell
!               nc_start=1

!               if((btest(state_up, cc_start+j-1).and..not.btest(state_up, nc_start+k-1))&
!                 &.and.abs(cell%inter(j,k))>1E-6) then
!                 new_up=ibclr(state_up,cc_start+j-1)
!                 new_up=ibset(new_up,nc_start+k-1)

!                 min_state=min(cc_start+j-1,nc_start+k-1)
!                 max_state=max(cc_start+j-1,nc_start+k-1)
!                 check_sites=(2_i64**(max_state) - 1_i64) - (2_i64**(min_state + 1) - 1_i64)

!                 ! Use the bin parity poppar() to return 1,-1 as required
!                 fermi_sign=1.0_dp-2.0_dp*real(poppar(iand(state_up,check_sites)),kind=dp)

!                 ! Find the matrix state which corresponds to new_up
!                 row=get_index(d_tot,hbt_tot,new_up,state_down)
!                 hmb(row,col)=hmb(row,col)+fermi_sign*cell%inter(j,k)*pexp*kexp
!                 ! Taking the transpose for reverse hopping
!                 ! Much like in htn, only fwrd hopping is defined by default
!                 hmb(col,row)=hmb(col,row)+fermi_sign*cell%inter(k,j)*conjg(pexp*kexp)
!               end if

!               if((btest(state_up, cc_start+j-1).and..not.btest(state_up, cc_start+k-1))&
!                 &.and.abs(cell%inter(j,k))>1E-6) then
!                 new_up=ibclr(state_up,cc_start+j-1)
!                 new_up=ibset(new_up,cc_start+k-1)

!                 min_state=min(cc_start+j-1,cc_start+k-1)
!                 max_state=max(cc_start+j-1,cc_start+k-1)
!                 check_sites=(2_i64**(max_state) - 1_i64) - (2_i64**(min_state + 1) - 1_i64)

!                 ! Use the bin parity poppar() to return 1,-1 as required
!                 fermi_sign=1.0_dp-2.0_dp*real(poppar(iand(state_up,check_sites)),kind=dp)

!                 ! Find the matrix state which corresponds to new_down
!                 row=get_index(d_tot,hbt_tot,state_up,new_down)
!                 hmb(row,col)=hmb(row,col)+fermi_sign*cell%inter(j,k)*pexp*kexp
!                 ! Taking the transpose for reverse hopping
!                 ! Much like in htn, only fwrd hopping is defined by default
!                 hmb(col,row)=hmb(col,row)+fermi_sign*cell%inter(k,j)*conjg(pexp*kexp)
!               end if
!             end if
!           end do
!         end do
!       end do
!     end do

!   end subroutine make_hmb

!   integer function get_index(d_tot,hbt_tot,target_up,target_down)
!     implicit none

!     integer(kind=i64),intent(in)                   :: d_tot
!     integer(kind=i64),dimension(d_tot,2),intent(in):: hbt_tot
!     integer(kind=i64),intent(in)                   :: target_up,target_down

!     integer:: i,j

!     get_index=0

!     do i=1,int(d_tot)
!       if(hbt_tot(i,1)==target_up)then
!         do j=i,int(d_tot)
!           if(hbt_tot(j,1)/=target_up) exit
!           if(hbt_tot(j,2)==target_down)then
!             get_index=j
!             return
!           end if
!         end do
!       end if
!     end do

!     if(get_index==0) stop 'Error finding index in get_index; see Hubbard'
!   end function get_index

 end program test