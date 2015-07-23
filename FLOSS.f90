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
		integer ::  cabin, cbbin
		character(len=4) :: rightanchor
		integer :: rightresiduenumber 
		character(len=4) ::  leftanchor
		integer :: leftresiduenumber
		real :: caseparation, cbseparation, cacross, cbcross
		character(len=4) :: pdbid 
		character(len=1) :: chainid
	end type distancetype

  integer, parameter :: top = 200, bottom = 2
		
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
		integer :: argnum, indexnum

		!for subroutine computeseparation

		!first element in abinsize is cabinsize, second is cbbinsize
		integer, dimension(top) :: abinsize
		integer, dimension(10) :: restotal
		integer :: res, filecount, totalsize
		!an allocatable array to store the information

		!initialization of necessary variables for subsequent use
		totalsize = 0
		res = 1
		filecount = 0

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
		call computeseparations(arrayres, arraydist, filecount, restotal, totalsize, abinsize)
		call cabinstructure(totalsize, abinsize, arraydist)

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

	subroutine computeseparations(arraya, arrayb, filecount, restotal, totalsize, abinsize)
		!access the xyz coordinates for the ca(1,2) and compute the distance between them in 3d space.EACH SEPARATION IS CALCULATED USING THE DISTANCE FORUMULA FOR A 3 COORDINATE SYSTEM
		type (residuetype), intent(inout), dimension(1000) :: arraya
		type (distancetype), intent(inout), allocatable :: arrayb(:)
		integer, intent(inout):: filecount, totalsize
		integer, dimension(10), intent(inout) :: restotal
		integer, dimension(top), intent(inout) :: abinsize
		integer :: cabin, cbbin, lastbin
		integer :: i, j, n, t, m, h
		integer :: ostat,cstat,fstat, numstat
		character(len=4) :: id
		real :: distca, distcb

		n = 0
		t = 1
		m = 1
		id = arraya(1)%pdbid
    lastbin = 0

    do i=1,size(abinsize,1)
      abinsize(i) = 0
    enddo

		allocate(arrayb(sum(restotal(1:filecount-1)**2/2) - sum(restotal(1:filecount-1)/2)))
		!this wont work if there is no available cb, determine for glycines

		do m=1, filecount-1
			h = t
			t = t + restotal(m)
			do i = h + 2 , t - 1
				do j = h, i - 2
					!CACA
			    distca = sqrt(sum((arraya(i)%cacoordinate(1:3)-arraya(j)%cacoordinate(1:3))**2))
					cabin = int(distca * 10)
					if(cabin > top) cycle
          if(cabin < bottom) cycle
          
					!CBCB
					distcb = sqrt(sum((arraya(i)%cbcoordinate(1:3)-arraya(j)%cbcoordinate(1:3))**2))
					cbbin = int(distcb * 10)
		
					if(cbbin > top) cycle !make this a parameter
          if(cbbin < bottom) cycle
					
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
					arrayb(n)%rightanchor = arraya(j)%aminoacid
					arrayb(n)%leftanchor = arraya(i)%aminoacid
					arrayb(n)%cabin = cabin
					arrayb(n)%cbbin = cbbin
					!WRITE THE CABIN NUMBER, CBBIN NUMBER, ATOM1NUMBER, ATOM2NUMBER, DISTANCES, AND PDBID ONLY LESS THAN 50Angstroms
          if(abinsize(cabin) == 100) then
            n = n - 1
            cycle
          endif

					call countbins(cabin, abinsize, totalsize)

				enddo
			enddo
		enddo
    !write(*,*) restotal(1:filecount-1)**2/2
	end subroutine computeseparations





	subroutine countbins(cabin, abinsize, totalsize)
		integer, intent(inout) ::  cabin, totalsize
		integer, dimension(top), intent(inout) :: abinsize
	!count bins will populate an array which holds integers corresponding to the number of elements which each bin will need to have capacity for
		abinsize(cabin) = abinsize(cabin) + 1 
		totalsize = totalsize + 1
	end subroutine countbins

!_______________integrated up to about this line_______________________

	subroutine cabinstructure(totalsize, abinsize, arrayb)
	!use this subroutine organizes the ca bin structure, filled with zeroes
		integer, intent(inout) :: totalsize
		integer, dimension(:), intent(inout) :: abinsize
		integer, dimension(size(abinsize,1)) :: abinsizecounter
		type (distancetype), intent(inout), dimension(totalsize) :: arrayb
		type (distancetype), dimension(totalsize) :: arrayc
		

		character(len=*), parameter :: infile="singlysorted.txt"		

		integer :: ostat,fstat,cstat
		integer :: ostat2,fstat2,cstat2
		integer :: lstat
		integer :: i, j, t, n
		integer :: cabin, cbbin

		t = 1 !give this counter a name
		abinsizecounter = abinsize

		!using the cabin structure, print to a file the records in bins organized by cabinstructure. include in the record the cabin number, cbbin number, cacb, cbca, anchoring residue by chain, orthogonal cordinates, and pdbid
		 
		!loop through the bins to find the locations where the beginning of the cabin record is

		open(unit=19, file=infile, form="formatted", iostat=ostat, access="direct", action="readwrite", status="replace", recl=62)
		!loop through the array structure arraydist (arrayb), calculate the cabin from the casep and cbbin from cbsep
		do i=1, totalsize
			cabin = int(arrayb(i)%caseparation * 10)
			cbbin = int(arrayb(i)%cbseparation * 10)
			t = sum(abinsize(1:cabin-1)) + abinsizecounter(cabin) 
      write(unit=19, fmt='( 2(I4), 2(A4, I4), 4F8.3, A4, A1)', iostat=fstat, rec=t) arrayb(i)
			abinsizecounter(cabin) = abinsizecounter(cabin) - 1
		enddo
		close(unit=19,iostat=cstat)

    call cbbinstructure(abinsize,totalsize)

	end subroutine cabinstructure



	subroutine cbbinstructure(abinsize, totalsize)
    integer, intent(inout) :: abinsize(:)
    integer, intent(inout) :: totalsize

    integer, dimension(top + 1) :: init
    integer, dimension(top) :: bbinsize, bbinsizecounter, bbinsizefinal
    integer :: i, j, linenum, t, m, cabin, cbbin, width
    integer :: ostat, fstat, cstat, ostat2, fstat2, cstat2

    type(distancetype) :: line, record

    character(len=*), parameter :: infile="singlysorted.txt"
    character(len=*), parameter :: outfile="doublysorted.txt"
    character(len=*), parameter :: tabfile="lookuptable.txt"

    
    init(1:top) = 0
    i = 1
    t = 0
    linenum = 1

    open(unit=19, file=infile, form="formatted", iostat=ostat, access="direct", action="readwrite", status="old", recl=62)
    open(unit=39, file=outfile, form="formatted", iostat=ostat2, access="direct", action="readwrite", status="replace", recl=62)
    open(unit=37, file=tabfile, form="formatted", iostat=ostat, access="direct", action="readwrite", status="replace", recl=1616)

    !go through and make the counts you need

    do m=1, top
      write(unit=37, fmt='( I8,200(I8) )', iostat=fstat, rec=m) 0, linenum, init
    enddo

    do i=1, top

      do m=1, top
        bbinsize(m) = 0
        bbinsizecounter(m) = 0
        bbinsizefinal(m) = 0  
      enddo

      width = abinsize(i)

      if(width == 0) then
        write(unit=37, fmt='( 2(I8),200(I8) )', iostat=fstat, rec=i) width, linenum, bbinsize
        cycle
      endif
      !within the width, build the bin, write the file, move to the next bin
      !write an error here for if i /= line%cabin
      do j=linenum, linenum + width - 1
        read(unit=19, fmt='( 2(I4), 2(A4, I4), 4F8.3, A4, A1)' , iostat = fstat, rec=j) line
        cbbin = line%cbbin
        bbinsize(cbbin) = bbinsize(cbbin) + 1
      enddo

        bbinsizecounter = bbinsize

      do k=1, size(bbinsize,1)
        bbinsizefinal(k) = sum(bbinsize(1:k))
        if(bbinsize(k) == 0) then
          bbinsizefinal(k) = 0
        endif
      enddo

      write(unit=37, fmt='( 2(I8),200(I8) )', iostat=fstat, rec=i) width, linenum, bbinsizefinal

      do j = linenum, linenum + width - 1
        read(unit=19, fmt='( 2(I4), 2(A4, I4), 4F8.3, A4, A1)' , iostat=fstat, rec=j) line
        cbbin = line%cbbin
        t = sum(bbinsize(1:cbbin)) - bbinsizecounter(cbbin)
        t = t + linenum
        write(unit=39, fmt='( 2(I4), 2(A4, I4), 4F8.3, A4, A1)', iostat=fstat2, rec=t) line
        bbinsizecounter(cbbin) = bbinsizecounter(cbbin) - 1
      enddo

      linenum = linenum + width
    enddo

    close(unit=37, iostat=cstat)
    close(unit=19, iostat=cstat)
    close(unit=39, iostat=cstat2)

  end subroutine cbbinstructure

!the cbbin values for each bin need to be saved, as does the cabinsize values in two separate files. save the cabin file as a cabin number, the cbbin number and the count. save the cabin file as it is structured, cabin and count

end program main


