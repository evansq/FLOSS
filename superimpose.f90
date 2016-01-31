!!!LOOP SUPERIMPOSE
!WRITE A PROGRAM WHICH SUPERIMPOSES TWO LOOPS, PROGRAM IS RUN ON ALL PROTEINS IN LIST(MOVING) CALLED IN REFERENCE TO INQUIRY LOOP(FIXED)***
!PSEUDO
!	IN MAIN:
!		1. RECIEVE IN FIXED AND MOVING LOOP ANCHOR COORDINATES
!		2. CONSTRUCT TRANSFORMATION MATRIX (CONSTRUCT TRANSLATION VECTOR AND ROTATION VECTOR)
!		3. SEND TO MODULE SUPERIMPOSE
!			MODULE SUPERIMPOSE HAS A FUNCTION DEFINETRANSFORMATION, WHERE THE COORDINATES OF THE MOVING LOOP ANCHOR ARE USED TO 
!			DEFINE THE TRANSFORMATION IN 3 DIMENSIONS
!		4. SEND TO FUNCTION SUPERIMPOSE, WHICH APPLIES TRANSFORMATION VECTOR TO ALL COORDINATES OF MOVING LOOP AND STORES IN A
!			 TEMPORARY ACCESS FILE WHICH GETS REWRITTEN ON CALL OF EVERY NEW INQUIRY. 
module superimpose
	type residuetype
		character(len=4) :: res
		integer :: resnum
		integer, dimension(3) :: cacoord, cbcoord 
	end type residuetype

	type looptype
		character(len=10) :: pdbid,chainid
		integer :: length
		type (residue), dimension(length) :: res
	end type looptype

end module superimpose


program main
	type (residuetype) :: residue
	type (looptype) :: fixed
	integer :: counter
	type (looptype), dimension(:), allocatable :: moving !allocatable to the list [count] 
	character(len=10) :: listof				!name of the list of loops in file argument 1
	
	integer :: i											!loop integers
	integer :: indexnum, fstat, ostat, cstat		!file handling variables
	character(len=68) :: listitem			!line of list file
	character(len=) :: line						!line of pdb

	integer, dimension(3,3) :: idrotvec, rotvec
	integer, dimension(3) :: idtransvec, transvec
!	process the commandline arguments for 
!		1. listofloops (moving loop)
!		2. inquiry loop pdbid
!		3. inquiry loop chain id
!		4. inquiry loop length
!		5. inquiry loop residue name anchor1
!		6. inquiry loop residue number anchor1
	
	counter = 0
	argnum = iargc()

	if(argnum /= 6) stop
	call getarg(1,listof)
	call getarg(2,fixed%pdbid)
	call getarg(3,fixed%chainid)
	call getarg(4,fixed%length)
	call getarg(5,residue%res)
	call getarg(6,residue%resnum)

	indexnum = index(listof,'.txt')
	if(indexnum == 0) then
		write(*,*) "the file", listof, &
		" could not be opened, check filesource error 1"
		stop
	endif

	open(unit=15, file=listof, iostat=ostat, & 
	access="sequential", 	action="read", status="old")

	do
		read(unit=15, fmt= '( A )', iostat=fstat) listitem
		if(fstat /= 0) exit
		counter = counter + 1
	enddo

	rewind 15
	allocate(moving(counter))

	do
		read(unit=15, fmt = '( A )', iostat=fstat) listitem
		if(fstat /= 0) exit
		!populate the list of possible loops, [write a file key]
	enddo
	
!	populate 		






end program main




