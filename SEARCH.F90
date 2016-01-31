!FLO SEARCH CREATED 2015
!RENSSELAER POLYTECHNIC INSTITUTE
!BY QUADIS EVANS
!FOR BYSTROFF LAB GROUP 
!TO BE INTEGRATED INTO INTERACTIVE ROSETTA

module search
  contains
  subroutine foo(casep, cbsep, camaxdev, cbmaxdev, limit, length, lendev)

    real, intent(in) :: casep, cbsep, camaxdev, cbmaxdev
    integer, intent(in) :: limit, length, lendev
    integer :: totalelements
    integer :: ostat
    character(len=*), parameter :: filename="binfile.txt"

    !there will be four arguments, casep, cbsep, max deviation ca, max deviation cb
    write(*,*) casep
		
    write(*,*) cbsep
		
    write(*,*) camaxdev
		
    write(*,*) cbmaxdev
		
    write(*,*) limit

    write(*,*) length

    write(*,*) lendev

    write(*,*) ""
		
    open(unit=29, file=filename, form="formatted", iostat=ostat, access="sequential", action="readwrite", status="replace")

    call buildlist(totalelements, casep, camaxdev, cbsep, cbmaxdev, limit)
    !call buildcbbinlist(bins, totalelements, cbsep, cbmaxdev)
    call printlist(totalelements,length, lendev)
    !generate the list of cabins to vet for distances 
  end subroutine foo

  subroutine buildlist(totalelements, casep, cadev, cbsep, cbdev, limit)
    !take in separation and deviation and pick out the cabins
    type distancetype
      integer::  cabin, cbbin
      character(len=4) :: rightanchor
      integer :: rightresiduenumber 
      character(len=4) ::  leftanchor
      integer :: leftresiduenumber
      real :: caseparation, cbseparation, cacross, cbcross
      character(len=4) :: pdbid 
      character(len=1) :: chainid
      integer :: length
    end type distancetype

    integer, parameter :: binmax = 200
    type(distancetype) :: bin
    real, intent(in) :: casep, cadev, cbsep, cbdev
    integer, intent(inout) :: totalelements
    integer, intent(in) :: limit
    integer :: cabin, cbbin
    integer, dimension(binmax + 2) :: line
    integer :: ostat,cstat,fstat
    integer :: ostat2,cstat2,fstat2
    integer :: fstat3
    integer :: i, l
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
    if(alow < 1 ) then !set less than bottom
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

    !write(*,*) cabin, alow, ahigh, cbbin, blow, bhigh
    !open the doubly sorted file
    open(unit=14, file="doublysorted.txt", form="formatted", iostat=ostat, access="direct", action="read", status="old", recl=66)

    !open the abinsize file
    open(unit=19, file="lookuptable.txt", form="formatted", iostat=ostat2, access="direct", action="read", status="old", recl=1616)

    !count the number of bin elements for allocation and find the start and end records

    do i=alow,ahigh
      read(unit=19, fmt='( 2(I8),200(I8) )', iostat=fstat2, rec=i) line
      start = line(blow + 2 - 1) + 1
      finish = line(bhigh + 2)
      totalelements = totalelements + (finish - start)
    enddo

    write(*,*) ""

    do i=alow,ahigh
      read(unit=19, fmt='( 2(I8),200(I8) )', iostat=fstat2, rec=i) line
      start = line(blow + 2 - 1) + 1
      finish = line(bhigh + 2) - 1
      astart = line(2)

      first = astart + start
      last = astart + finish

      do j=first, last
        if(l == totalelements) then
          exit
        else if(l == limit) then
          totalelements = limit
          exit
        endif
        read(unit=14, fmt="( 2(I4), 2(A4, I4), 4F8.3, A4, A1, I4)", iostat=fstat, rec=j) bin
        write(unit=29, fmt="( 2(I4), 2(A4, I4), 4F8.3, A4, A1, I4)", iostat=fstat3) bin
        l = l + 1
      enddo

    enddo

    close(unit=14, iostat=cstat)
    close(unit=19, iostat=cstat2)

  end subroutine buildlist

  subroutine printlist(totalelements, length, lendev)
  !just to check to se if the list is properly allocating
    type distancetype
      integer :: cabin, cbbin
      character(len=4) :: rightanchor
      integer :: rightresiduenumber 
      character(len=4) :: leftanchor
      integer :: leftresiduenumber
      real :: caseparation, cbseparation, cacross, cbcross
      character(len=4) :: pdbid 
      character(len=1) :: chainid
      integer :: length
    end type distancetype

    integer, intent(in) :: totalelements, length, lendev
    integer :: cstat,fstat
    type(distancetype) :: bin
    integer :: i

    rewind 29

    do i=1,totalelements !change to read until eof
      read(unit=29, fmt="( 2(I4), 2(A4, I4), 4F8.3, A4, A1, I4)", iostat=fstat) bin
      if(bin%length < length - lendev) cycle
      if(bin%length > length + lendev) cycle
      write(*,*) bin
    enddo 

    close(unit=29, iostat=cstat, status='delete')

  end subroutine printlist
end module search




program main
	integer :: argnum,indexnum

	argnum = iargc()
	if(argnum /= 6) then
		write(*,*) "please specify all variables for input: 6 variables, see README"
	call getarg(1,allpdbs)
	indexnum = index(allpdbs,txt)
	if(indexnum == 0) then
		write(*,*) "the file", allpdbs, &
		" could not be opened, proceed 	with a txt file"
		stop
	endif
  use search
  implicit none
end program main
!run through the final set of distances, get all backbone coordinates, store in new file
