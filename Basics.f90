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

  subroutine get_params(size,filling,phase_case,a_val,epsn_val,t_val)
    implicit none
    integer,intent(OUT)      :: size,filling
    real(kind=dp),intent(OUT):: phase_case,epsn_val,t_val,a_val

    character(len=100):: ifl
    character(len=100):: line,param
    integer           :: file_unit=15,istat=0
    logical           :: lexist

    !reading input filename from make file
    read(*,*) ifl
    print*, 'file input read as:',ifl

    inquire(file=ifl,exist=lexist)
    !params retain default initial vals if no input file is provided
    if(.not.lexist) return !'no input, default params used

!---------------------------------------------------------------------------!
!                    INITIALISING DEFAULT VALUES                            !
!---------------------------------------------------------------------------!

    size      =20    !/Number of sites
    filling   =10    !/Number of energy bands filled with electrons
    phase_case=0.0_dp!To be multiplied by pi to determine phase/antiphase etc
    a_val     =1.0_dp!Atomic Spacing/Amstrongs
    epsn_val  =0.0_dp!On site energy/eV
    t_val     =-1.0_dp!Hopping value (nearest neighbour)

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
      case('SIZE')
        read(line(index(line,'=')+1:),*,iostat=istat) size
        if(istat/=0) stop 'invalid NUMBER OF SITES'
      case('FILLING')
        read(line(index(line,'=')+1:),*,iostat=istat) filling
        if(istat/=0) stop 'invalid FILLING'
      case('PHASE')
        read(line(index(line,'=')+1:),*,iostat=istat) phase_case
        if(istat/=0) stop 'invalid PHASE'
      case('SPACING')
        read(line(index(line,'=')+1:),*,iostat=istat) a_val
        if(istat/=0) stop 'invalid SPACING'
      case('ON-SITE')
        read(line(index(line,'=')+1:),*,iostat=istat) epsn_val
        if(istat/=0) stop 'invalid ON-SITE'
      case('HOPPING')
        read(line(index(line,'=')+1:),*,iostat=istat) t_val
        if(istat/=0) stop 'invalid HOPPING'
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
    print*,'                  SIZE   = ',size
    print*,'                  FILLING= ',filling
    print*,'                  PHASE  = ',phase_case
    print*,'                  SPACING= ',a_val
    print*,'                  ON-SITE= ',epsn_val
    print*,'                  HOPPING= ',t_val
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

