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

  subroutine get_dist(size, nn_map, dist)
    implicit none
    integer, intent(in)                  :: size
    logical, dimension(:,:), intent(in)  :: nn_map
    integer, dimension(:,:), intent(out) :: dist
    integer :: i, j, k
    integer, parameter :: INF = 999999 ! Represents no connection

    ! Initialising the Distance Matrix
    do i = 1, size
      do j = 1, size
        if (i == j) then
          dist(i, j) = 0           ! Distance to itself is 0
        else if (nn_map(i, j).eqv..false.) then
          dist(i, j) = 1           ! 1st NN(1 hop)
        else
          dist(i, j) = INF         ! Unknown distance
        end if
      end do
    end do

    ! Floyd-Warshall 
    do k = 1, size
      do i = 1, size
        do j = 1, size

          if (dist(i, k) + dist(k, j) < dist(i, j)) then
            dist(i, j) = dist(i, k) + dist(k, j)

          end if
        end do
      end do
    end do

  end subroutine get_dist

  subroutine get_params(sys_size,sites,filling,phase_case,a_val,&
                      &cell_size,nbr_map,nbr_inter,nbr_range,phi_max,k_max)
    implicit none
    integer,intent(OUT)                                 :: sys_size,filling
    integer,intent(OUT)                                 :: cell_size,nbr_range,phi_max,k_max
    real(kind=dp),intent(OUT)                           :: phase_case,a_val
    real(kind=dp),dimension(:,:),allocatable,intent(OUT):: sites
    logical,dimension(:,:),allocatable,intent(out)      :: nbr_map,nbr_inter

    character(len=100)                    :: ifl
    character(len=100)                    :: line,param
    integer                               :: node1,node2,i,cell_num,n
    integer                               :: file_unit=15,istat=0,dimer
    logical                               :: lexist
    real(kind=dp),dimension(:),allocatable:: tmp_sites
    real(kind=dp)                         :: AAH

    !reading input filename from make file
    read(*,*) ifl
    print*, 'file input read as:',ifl

    inquire(file=ifl,exist=lexist)
    !params retain default initial vals if no input file is provided
    if(.not.lexist) return !'no input, default params used

!---------------------------------------------------------------------------!
!                    INITIALISING DEFAULT VALUES                            !
!---------------------------------------------------------------------------!

    phi_max   =1E4
    k_max     =1
    sys_size  =1    !/Number of sites
    filling   =10    !/Number of energy bands filled with electrons
    phase_case=0.0_dp!To be multiplied by pi to determine phase/antiphase etc
    a_val     =1.0_dp!Atomic Spacing/Amstrongs
    cell_size =8!number of sites within a unit cell
    nbr_range =1!hopping neighbours
    AAH       =0.0_dp
    dimer     =1

    allocate(nbr_map(cell_size,cell_size),stat=istat)
    if(istat/=0) stop 'error allocating default nbr_map'
    allocate(nbr_inter(cell_size,cell_size),stat=istat)
    if(istat/=0) stop 'error allocating default nbr_inter'
    allocate(sites(cell_size,nbr_range+1),stat=istat)
    if(istat/=0) stop 'error allocating default sites'
    
    nbr_inter=.true.
    nbr_map=.true.
    nbr_map(1,cell_size)=.false.
    nbr_map(cell_size,1)=.false.
    do i=2,cell_size-1
      nbr_map(i,i+1)=.false.
      nbr_map(i,i-1)=.false.
    end do

    ! Initial assumptions for hopping assume homogeneous on-site and hopping terms
    sites(:,1) =0.0_dp
    sites(:,2:)=1.0_dp

    open (unit=file_unit,file=ifl,status="old",action="read",iostat=istat)
    if (istat/=0) STOP 'error openning input file'

    !DO LOOP over lines in the input file
    do while(istat==0)

      !READS the current line into the line character variable
      read(unit=file_unit,fmt='(a)',iostat=istat) line

      !verifies that the line isn't a comment;removes leading blank space
      line=adjustl(line)
      if (line(1:1) == "#") cycle

      param=trim(line(1:index(line,'=')-1))

      !Verify the keywords
      !TRIM() removes the trailing whitespace
      select case(trim(param))
        
      case('PHIMAX')
        read(line(index(line,'=')+1:),*,iostat=istat) phi_max
        if(istat/=0) stop 'invalid PHIMAX'

      case('KMAX')
        read(line(index(line,'=')+1:),*,iostat=istat) k_max
        if(istat/=0) stop 'invalid KMAX'

      case('SIZE')
        read(line(index(line,'=')+1:),*,iostat=istat) sys_size
        if(istat/=0) stop 'invalid NUMBER OF CELLS'

      case('CELLSIZE')
        print*, 'RANGE should already be defined to get desired system'
        read(line(index(line,'=')+1:),*,iostat=istat) cell_size
        if(istat/=0) stop 'invalid CELLSIZE'

        deallocate(sites,stat=istat)
        if(istat/=0) stop 'error switching size of nbr_map from default'
        allocate(sites(cell_size,2*nbr_range+1),stat=istat)
        sites=0.0_dp

        if(istat/=0) stop 'error allocating nbr_map with CELLSIZE'
        deallocate(nbr_map,stat=istat)
        if(istat/=0) stop 'error switching size of nbr_map from default'
        allocate(nbr_map(cell_size,cell_size),stat=istat)
        if(istat/=0) stop 'error allocating nbr_map with CELLSIZE'
        nbr_map=.true.

        deallocate(nbr_inter,stat=istat)
        if(istat/=0) stop 'error switching size of nbr_inter from default'
        allocate(nbr_inter(cell_size,cell_size),stat=istat)
        if(istat/=0) stop 'error allocating nbr_map with CELLSIZE'
        nbr_inter=.true. !default assumption that none of the sites cannot see each other

      case('SITE')
        allocate(tmp_sites(2*nbr_range+1),stat=istat)
        if(istat/=0) stop 'error allocating tmp_sites'

        print*, 'RANGE and CELLSIZE MUST be defined to get desired system.'
        read(line(index(line,'=')+1:),*,iostat=istat) cell_num, tmp_sites(:)
        if(istat/=0) stop 'invalid SITE'

        if (cell_num > 0 .and. cell_num <= cell_size) then
            sites(cell_num,:)=tmp_sites(:)
        else
            print*, 'WARNING: SITE index ', cell_num, ' is out of bounds! Skipping.'
        end if

        deallocate(tmp_sites,stat=istat)
        if(istat/=0) stop 'error deallocating tmp_sites'

      case('FILLING')
        read(line(index(line,'=')+1:),*,iostat=istat) filling
        if(istat/=0) stop 'invalid FILLING'

      case('PHASE')
        read(line(index(line,'=')+1:),*,iostat=istat) phase_case
        if(istat/=0) stop 'invalid PHASE'

      case('SPACING')
        read(line(index(line,'=')+1:),*,iostat=istat) a_val
        if(istat/=0) stop 'invalid SPACING'

      case('RANGE')
        read(line(index(line,'=')+1:),*,iostat=istat) nbr_range
        if(istat/=0) stop 'invalid RANGE'

      case('CONNECT')
        !Read the two integers separated by a space
        read(line(index(line,'=')+1:),*,iostat=istat) node1, node2
        if(istat/=0) stop 'invalid CONNECT parameters'
        ! Check bounds to prevent array crashing
        if (node1 > 0 .and. node1 <= cell_size .and. node2 > 0 .and. node2 <= cell_size) then
          nbr_map(node1, node2) = .false.
          nbr_map(node2, node1) = .false.
        else
          print*, 'WARNING: CONNECT out of bounds: ', node1, node2
        end if

      case('CONNECTN')
        read(line(index(line,'=')+1:),*,iostat=istat) node1, node2
        if(istat/=0) stop 'invalid CONNECT parameters'
        if (node1 > 0 .and. node1 <= cell_size .and. node2 > 0 .and. node2 <= cell_size) then
          ! NOTE: This is NOT symmetric! ONLY forward hop has been recorded
          nbr_inter(node1, node2) = .false.
        else
          print*, 'WARNING: CONNECT_NEXT out of bounds: ', node1, node2
        end if

      !AAH case must come last
      case('AAH')
        read(line(index(line,'=')+1:),*,iostat=istat) AAH,dimer
        if(istat/=0) stop 'invalid AAH parameter'

        n=1

        do i=1,size(sites,1)
          sites(i,1)=sites(i,1)*cos(2.0_dp*real_pi*AAH*n+phase_case*real_pi)
          if(mod(i,dimer)==0) n=n+1
        end do

      case('EXIT')
        print*,'input file read'
        exit
      case default
        print*, 'WARNING:', trim(param), ' is an invalid parameter'
      end select
    end do

    close (unit = file_unit , iostat = istat)
    if (istat /= 0) STOP 'error closing input file'

!---------------------------------------------------------------------------!
!                   PRINTING VALS WHICH WILL BE USED                        !
!---------------------------------------------------------------------------!
    
    print*,'!--------------------------------------------------------------!'
    print*,'!PARAMETERS USED:                                              !'
    print*,'!--------------------------------------------------------------!'
    print*,'                  PHIMAX  = ',phi_max
    print*,'                  KMAX    = ',k_max
    print*,'                  SIZE    = ',sys_size
    print*,'                  CELLSIZE= ',cell_size
    print*,'                  FILLING = ',filling
    print*,'                  PHASE   = ',phase_case
    print*,'                  SPACING = ',a_val
    print*,'                  RANGE   = ',nbr_range
    print*,'                  AAH     = ',AAH,dimer
    print*,'!--------------------------------------------------------------!'
    do i=1,size(sites,1)
      print*,'               SITE ',i,' = ' ,sites(i,:)
    end do
    print*,'!--------------------------------------------------------------!'
    do i=1,size(nbr_map,1)
      print*,'                  CONNECT = ',nbr_map(i,:)
    end do
    print*,'!--------------------------------------------------------------!'
    do i=1,size(nbr_inter,1)
      print*,'              CONNECT_NBR = ',nbr_inter(i,:)
    end do
    print*,'!--------------------------------------------------------------!'
  end subroutine get_params


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

    JOBZ = 'V'
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

    if(INFO/=0) STOP 'Zheev error'

    deallocate(WORK,stat=ierr)
    if(ierr/=0) stop 'error deallocating WORK'

    deallocate(RWORK,stat=ierr)
    if(ierr/=0) stop 'error deallocating RWORK'

    !print*, W(1)

  end subroutine zheev_evals

end module constants

