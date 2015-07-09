!FLOSS CREATED 2015
!RENSSELAER POLYTECHNIC INSTITUTE
!BY QUADIS EVANS
!FOR BYSTROFF LAB GROUP 
!TO BE INTEGRATED INTO INTERACTIVE ROSETTA
!PURPOSE: SORT PDB FILES BY LOOP ANCHOR CA and CB SEPARATION A
!COMPARE THE SPEED AND ACCURACY OF THIS METHOD TO THE ALTERNATE PROGRAM: PDB FILES SORTED BY THE CROSS-DISTANCES BETWEEN THE ANCHORING RESIDUES CA/CB. 

!NOTES: 
!Length Expansion and Contraction
!FOR THE NEXT STEP: COMPUTE THE SEPARATION DISTANCES BETWEEN THE CA ATOMS
!SORT BY CA, BIN SIZES ARE .1A , COUNT THE NUMBER IN EACH BIN, WRITE THE INDEX (X=X+SIZE OF BIN),SIZEOFBIN+1
!THE PROGRAM WILL START FOR EACHPDB FILE, LOCATING THE START OF THE ATOM ANNOTATIONS AND DECLARING THAT VARIABLE RECSTART. THEN, TAKING LOOPING THROUGH THE AVAILABLE ATOMS ONE AT A TIME, INCREMENTING RECSTART. WITHIN THE INCREMENTING RECSTART LOOP, INCREMENT CURRENT REC TO DEFINE THE SECOND "ANCHOR" ATOM CURRENTREC. WHEN READING GETS TO THE END OF THE FILE, THE REC START SHOULD INCREMENT AND RESTART THE PROCESS. 
!
!CALL THE PROGRAM ONCE ON A FILE WITH A LIST OF PDB FILES, THEN CALL EACH SUBROUTINE ON EACH PDB FILE INDEPENDENTLY, KEEPING NECESSARY GLOBAL VARIABLES TO CHARACTERIZE THE OVERAL STRUCTURE
!rewrite captured info as defined types residuetype(int resnumber, char chainid, char pdbid, char aa, real(3) cacoord,cbcoord) and dist(int resnum1, int resnum2, real caca,cbcb,cacb,cbca, char pdbid, char chainid) as globals (construct, declare as allocatable)

program main
!the type definitions will remain global because i want them to
	type residuetype
		integer :: residuenumber
		character(len=1) :: chainid
		character(len=4) :: pdbid, aminoacid
		real,dimension(3) :: cacoordinate, cbcoordinate
	end type residuetype

	type distancetype
		character(len=4) :: pdbid
		character(len=1) :: chainid
		integer :: rightresiduenumber, leftresiduenumber
		real :: caseparation, cbseparation, cacross, cbcross
	end type distancetype
!the type declarations will remain global because i want them to
!return and make these allocatable, but for now, initialize them to a very large number
!return and initialize these elements to 0
!
! ex: arraydist(:)%[element]) = " " ... in loop
!
	call build

	contains
	
	subroutine build
		type (residuetype), dimension(1000) :: arrayres
		type (distancetype), dimension(:), allocatable :: arraydist

		!for the current subroutine
		character(len=4) :: pdbid
		character(len=1) :: chainid
		character(len=200) :: pdbfilename, allpdbs
		integer :: ostat,cstat,fstat
		integer :: argnum, indexnum, totaldist

		!for subroutine computeseparation

		!first element in abinsize is cabinsize, second is cbbinsize
		integer, dimension(500000) :: abinsize
		integer, dimension(10) :: restotal
		integer :: res, filecount, totalsize
		!an allocatable array to store the information

		!initialization of necessary variables for subsequent use
		totalsize = 0
		res = 1
		filecount = 0
		totaldist = 0

		!process the commandline arguments
		argnum = iargc()
		if(argnum /= 1) stop 
		call getarg(1,allpdbs)
		indexnum = index(allpdbs,'.txt')
		if(indexnum == 0) then
			write(*,*) "the file", allpdbs, " could not be opened, proceed 	with a txt file"
			stop
		endif
	
		!read through the command line file for each pdbid and chainid
		open(unit=15, file=allpdbs, iostat=ostat, access="sequential", 	action="read", status="old")

		do
			filecount = filecount + 1	
			read(unit=15,fmt="( A4,1X,A1 )", iostat=fstat) pdbid, chainid
			pdbfilename = pdbid//'.pdb' !add in a pathway specifier
			if(fstat /= 0) exit
			call readinpoints(pdbfilename, chainid, arrayres, restotal, res, filecount)
		enddo
		close(unit=15, iostat=cstat)
		call computeseparations(arrayres, arraydist, filecount, restotal, totalsize, abinsize, totaldist)
		call cabinstructure(totalsize, abinsize, arraydist, totaldist)

	end subroutine build

	subroutine readinpoints(pdbfilename, chain2, arraya, restotal, res, filecount)
	!readinpoints takes in the pdb files,reads them for the xyz location of the ca(1,2)/cb(1,2) and passes them to computeseparations to determine the values for each index
		character(len=8), intent(in) ::pdbfilename
		character(len=1), intent(in) :: chain2
		type (residuetype), intent(inout), dimension(1000) :: arraya
		integer, dimension(10), intent(inout) :: restotal
		integer, intent(inout) :: res, filecount
		
		integer :: n, ostat, fstat, cstat
		character(len=4) :: atomtype
		character(len=6) :: rectype
		character(len=3) :: residue
		character(len=78) :: line	
		character(len=1) :: chain
		character(len=4) :: pdbid
	
		pdbid = pdbfilename(1:4)
		n = 0

		open(unit=11, file=pdbfilename, iostat=ostat, access="sequential", action="readwrite", status="old")

		do 
			!for all atoms in file | read the first anchor
			read(unit=11, fmt="( A )", iostat=fstat) line
			if (fstat /= 0) exit !condition for end of file
			read(line(1:6), fmt="( A )") rectype
			if (rectype /= "ATOM  ") cycle
			read(line(13:16), fmt="( A )") atomtype
			if (atomtype /= " CA ") cycle !check for CA !add space to the other side, capture 4
			read(line(18:20), fmt="( A )") residue
			if(residue == "GLY ") cycle !ignore the glycines
			read(line(22:22), fmt="( A )") chain
			if(chain /= chain2) exit
			arraya(res+n)%aminoacid = residue
			read(line(23:26), fmt="( I4 )") arraya(res+n)%residuenumber
			read(line(31:54) , fmt="( 3F8.3 )") arraya(res+n)%cacoordinate(1:3)
			arraya(res+n)%pdbid = pdbid
			arraya(res+n)%chainid = chain
			!FINDS THE CB COORDINATES AHEAD OF THE CURRENT LINE
			do
				read(unit=11, fmt="( A )", iostat=fstat) line
				read(line(13:16), fmt="( A )") atomtype
				if (atomtype == " CB ") then
					read(line(31:54) , fmt="( 3F8.3 )") arraya(res+n)%cbcoordinate(1:3)
					exit
				endif
			enddo 
			n = n + 1 
		enddo
		res = res + n
		restotal(filecount) = n
		close(unit=11, iostat=cstat)	
	end subroutine readinpoints	

	subroutine computeseparations(arraya, arrayb, filecount, restotal, totalsize, abinsize, totaldist)
		!access the xyz coordinates for the ca(1,2) and compute the distance between them in 3d space.EACH SEPARATION IS CALCULATED USING THE DISTANCE FORUMULA FOR A 3 COORDINATE SYSTEM
		type (residuetype), intent(inout), dimension(1000) :: arraya
		type (distancetype), intent(inout), allocatable :: arrayb(:)
		integer, intent(inout):: filecount, totalsize, totaldist
		integer, dimension(10), intent(inout) :: restotal
		integer, dimension(500000), intent(inout) :: abinsize
		integer :: cabin, cbbin
		integer :: i, j, n, t, m, h
		integer :: ostat,cstat,fstat, numstat
		character(len=4) :: id
		real :: distca, distcb

		n = 0
		t = 1
		h = 1
		m = 1
		id = arraya(1)%pdbid
		totaldist = (sum(restotal(1:filecount)**2)/2)
		allocate(arrayb(totaldist))

		open(unit=29, file="prelim.txt", form="formatted", iostat=ostat, access="sequential", action="readwrite", status="replace")
		!this wont work if there is no available cb, determine for glycines

		do m=1,filecount - 1
			write(*,*) m
			h = t
			t = t + restotal(m)
			do i = h , t - 2
				do j = i + 2, t
					!CACA
					!arrayd(i,j,1) = sqrt(((arraya(i,1)-arraya(j,1))**2)+((arraya(i,2)-arraya(j,2))**2)+((arraya(i,3)-arraya(j,3))**2))
			    distca = sqrt(sum((arraya(i)%cacoordinate(1:3)-arraya(j)%cacoordinate(1:3))**2))
					cabin = int(distca * 10)
					if(cabin > 500) cycle
					!CBCB
					distcb = sqrt(sum((arraya(i)%cbcoordinate(1:3)-arraya(j)%cbcoordinate(1:3))**2))
					cbbin = int(distcb * 10)
		
					if(cbbin > 500) cycle !make this a parameter
					
					n = n + 1

					arrayb(n)%cbseparation = distcb
					arrayb(n)%caseparation = distca

					!CACB
					arrayb(n)%cacross = sqrt(sum((arraya(i)%cacoordinate(1:3)-arraya(j)%cbcoordinate(1:3))**2))

					!CBCA
					arrayb(n)%cbcross = sqrt(sum((arraya(i)%cbcoordinate(1:3)-arraya(j)%cacoordinate(1:3))**2))

					!save the residue numbers (left and right)
					arrayb(n)%pdbid = arraya(i)%pdbid
					arrayb(n)%chainid = arraya(i)%chainid
					arrayb(n)%rightresiduenumber = arraya(j)%residuenumber
					arrayb(n)%leftresiduenumber = arraya(i)%residuenumber
					!WRITE THE CABIN NUMBER, CBBIN NUMBER, ATOM1NUMBER, ATOM2NUMBER, DISTANCES, AND PDBID ONLY LESS THAN 50Angstroms
					call countbins(cabin, abinsize, totalsize)
					!write(unit=29, fmt='( I4, I4, A4, A1, I4, I4, 4F8.3 )', iostat=fstat) cabin, cbbin, arrayb(n)				
					!write (*,*) arrayb(n)%pdbid, n
				enddo
			enddo
		enddo
		write(*,*) n
	end subroutine computeseparations





	subroutine countbins(cabin, abinsize, totalsize)
		integer, intent(inout) ::  cabin, totalsize
		integer, dimension(500000), intent(inout) :: abinsize
	!count bins will populate an array which holds integers corresponding to the number of elements which each bin will need to have capacity for
		abinsize(cabin) = abinsize(cabin) + 1 
		totalsize = totalsize + 1
		!write(*,*) cabin, abinsize(cabin), totalsize
	end subroutine countbins

!_______________integrated up to about this line_______________________

subroutine cabinstructure(totalsize, abinsize, arrayb, totaldist)
	!use this subroutine organizes the ca bin structure, filled with zeroes
		integer, intent(inout) :: totalsize
		integer, dimension(:), intent(inout) :: abinsize
		integer, dimension(size(abinsize,1)) :: abinsizecounter
		integer, intent(in):: totaldist
		type (distancetype), intent(inout), dimension(totaldist):: arrayb
		

		character(len=*), parameter :: filename="diracloop.txt"		
		character(len=53) :: line

		integer :: ostat,fstat,cstat
		integer :: ostat2,fstat2,cstat2
		integer :: lstat
		integer :: i, j, t, last
		integer :: cabin, cbbin

		t = 1 !give this counter a name
		abinsizecounter = abinsize

		!using the cabin structure, print to a file the records in bins organized by cabinstructure. include in the record the cabin number, cbbin number, cacb, cbca, anchoring residue by chain, orthogonal cordinates, and pdbid
		 
		!loop through the bins to find the locations where the beginning of the cabin record is
		open(unit=19, file=filename, form="formatted", iostat=ostat, access="direct", action="readwrite", status="replace", recl=57)
		
		!loop through the array structure arraydist (arrayb), cl=alculate the cabin from the casep and cbbin from cbsep
		do i=1, size(arrayb,1)
			!read(unit=29, fmt='( I4, I4, A4, A1, I4, I4, 4F8.3 )', iostat=fstat2) cabin, cbbin, binelement
			cabin = int(arrayb(i)%caseparation * 10)
			cbbin = int(arrayb(i)%cbseparation * 10)
			if(cabin == 0) cycle
			t = sum(abinsize(1:cabin)) - abinsize(cabin) + abinsizecounter(cabin) 
			write(unit=19, fmt='( I4, I4, A4, A1, I4, I4, 4F8.3 )', iostat=fstat, rec=t) cabin, cbbin, arrayb(i)
			abinsizecounter(cabin) = abinsizecounter(cabin) - 1
			write(*,*) cabin, abinsizecounter(cabin), t
		enddo
		close(unit=19, iostat=cstat)
	end subroutine cabinstructure

!the next thing you need to do is to set up the cbbinstructure the same way you did with the cabinstructure. repeat the last two subroutines but for the second filter

end program main
