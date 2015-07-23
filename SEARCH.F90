!FLO SEARCH CREATED 2015
!RENSSELAER POLYTECHNIC INSTITUTE
!BY QUADIS EVANS
!FOR BYSTROFF LAB GROUP 
!TO BE INTEGRATED INTO INTERACTIVE ROSETTA


program main
  type distancetype
		integer ::  cabin, cbbin
		character(len=4) :: rightanchor
		integer :: rightresiduenumber 
		character(len=4) ::  leftanchor
		integer :: leftresiduenumber
		real :: caseparation, cbseparation, cacross, cbcross
		character(len=4) :: pdbid 
		character(len=1) :: chainid
	end type distancetype

  character(len=8) :: args
  integer :: argnum
  real :: casep, cbsep, camaxdev, cbmaxdev
  type(distancetype), allocatable :: bins(:)
  integer :: totalelements
  integer, parameter :: binmax = 200

  argnum = iargc()

	if(argnum /= 4) stop 
  !there will be four arguments, casep, cbsep, max deviation ca, max deviation cb
	call getarg(1,args)
  write(*,*) args
  read(args(1:5), fmt='(F8.3)') casep
  call getarg(2,args)
  write(*,*) args
  read(args(1:5), fmt='(F8.3)') cbsep
  call getarg(3,args)
  write(*,*) args
  read(args(1:5), fmt='(F8.3)') camaxdev
  call getarg(4,args)
  write(*,*) args
  read(args(1:5), fmt='(F8.3)') cbmaxdev
  
  !write(*,*) casep, cbsep, camaxdev, cbmaxdev
	!write error message for arguments
  call buildlist(bins, totalelements, casep, camaxdev, cbsep, cbmaxdev)
  !call buildcbbinlist(bins, totalelements, cbsep, cbmaxdev)
  call printlist(totalelements, bins)
  !generate the list of cabins to vet for distances
  contains
  
  subroutine buildlist(bins, totalelements, casep, cadev, cbsep, cbdev)
    !take in separation and deviation and pick out the cabins

    real, intent(in) :: casep, cadev, cbsep, cbdev
    integer, intent(inout) :: totalelements
    integer :: cabin, cbbin
    type(distancetype), allocatable, intent(inout) :: bins(:)
    integer, dimension(202) :: line
    integer :: ostat,cstat,fstat
    integer :: ostat2,cstat2,fstat2
    integer :: i, n, l
    integer :: alow, ahigh, blow, bhigh, first, last, start, finish, astart
    
    l = 1
    totalelements = 0
    start = 0
    finish = 0
    !calculate the lower and upper bounds for bin selection
    cabin = ceiling(casep * 10 - 0.5)
    alow = ceiling(casep*10 - ceiling(cadev*10 - 0.5) - 0.5)
    ahigh = ceiling(casep*10 + ceiling(cadev*10 - 0.5) - 0.5)

    cbbin = ceiling(cbsep * 10 - 0.5)
    blow = ceiling(cbsep*10 - ceiling(cbdev*10 - 0.5) - 0.5)
    bhigh = ceiling(cbsep*10 + ceiling(cbdev*10 - 0.5) - 0.5)

    !conditionals
    if(alow < 1) then
      alow = 1
    endif
    
    if(ahigh > binmax) then
      ahigh = binmax
    endif

    if(blow < 1) then
      blow = 1
    endif

    if(bhigh > binmax) then
      bhigh = binmax
    endif

    !open the doubly sorted file
    open(unit=14, file="doublysorted.txt", form="formatted", iostat=ostat, access="direct", action="read", status="old", recl=62)

    !open the abinsize file
    open(unit=19, file="lookuptable.txt", form="formatted", iostat=ostat2, access="direct", action="read", status="old", recl=1616)

    !count the number of bin elements for allocation and find the start and end records

    do i=alow,ahigh
      read(unit=19, fmt='( 2(I8),200(I8) )', iostat=fstat2, rec=i) line
      start = line(blow + 2 - 1) + 1
      finish = line(bhigh + 2)
      totalelements = totalelements + (finish - start)
    enddo

    write(*,*) totalelements
    !allocate the bin array
    if(totalelements > 100) then
      totalelements = 100
    endif

    allocate(bins(totalelements))

    do i=alow,ahigh
      read(unit=19, fmt='( 2(I8),200(I8) )', iostat=fstat2, rec=i) line
      
      start = line(blow + 2 - 1) + 1
      finish = line(bhigh + 2) - 1
      astart = line(2)
        
      first = astart + start
      last =  astart + finish

      do j=first, last
        if(l == 100) exit
        read(unit=14, fmt='( 2(I4), 2(A4, I4), 4F8.3, A4, A1)', iostat=fstat, rec=j) bins(l)
        write(*,*) bins(l)
        l = l + 1
      enddo

    enddo

    close(unit=14, iostat=cstat)
    close(unit=19, iostat=cstat2)

    !export to build cbbin filtered list

  end subroutine buildlist

  subroutine printlist(totalelements, onedimen)
  !just to check to se if the list is properly allocating
    integer, intent(in) :: totalelements
    type(distancetype), intent(in) :: onedimen(totalelements)
    integer :: i

    do i=1,totalelements - 1
      write(*,*) onedimen(i)
    enddo 
  
  end subroutine printlist

end program main
