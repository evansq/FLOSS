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
module floss

	type residuetype 																			!type definition which stores the information for individual residues
		integer :: residuenumber
		character(len=1) :: chainid
		character(len=4) :: pdbid, aminoacid
		real,dimension(3) :: cacoordinate, cbcoordinate
	end type residuetype

	type distancetype																			!type definition which stores the information for paired residue distances
		integer ::  cabin, cbbin
		character(len=4) :: rightanchor
		integer :: rightresiduenumber 
		character(len=4) ::  leftanchor
		integer :: leftresiduenumber
		real :: caseparation, cbseparation, cacross, cbcross
		character(len=4) :: pdbid 
		character(len=1) :: chainid
		integer :: length
	end type distancetype

	contains
		subroutine getpdbfile(pdbfilename)
			integer :: i
			character(len=8), intent(in) :: pdbfilename

			call execute_command_line("wget http://www.rcsb.org/pdb/files/"//pdbfilename//" -nH", exitstat=i)
		end subroutine getpdbfile

		subroutine filecounter(pdbfilename, chainid, top, bottom, filemax, filecount, restotal, abinsize)
		!!!!!----IN SUBROUTINE FILECOUNTER,CALL THE SUBROUTINE TO WGET THE FILE, EXECUTE ON WGOTTEN FILE----
		character(len=8), intent(in) :: pdbfilename
		character(len=1), intent(in) :: chainid
		integer, intent(in) :: top, filemax, filecount, bottom
		integer, dimension(filemax), intent(inout) :: restotal
		integer, dimension(top), intent(inout) :: abinsize
		real,dimension(3) :: rightanchor, leftanchor				!coordinates of left and right anchoring residues
		real :: distca																			!real value for computed distance
		integer :: cabin																		!c alpha bin number
		integer :: n, ostat, fstat, cstat		
		integer :: i, j, k																	!n: counter for number of residues in file being used
		character(len=4) :: atomtype												!atom type	
		character(len=6) :: rectype													!record type
		character(len=3) :: residue													!residue name (3 letter)
		character(len=78) :: line														!pdb record
		character(len=1) :: chain														!chain id at record

		n = 0

		call getpdbfile(pdbfilename)
		open(unit=11, file=pdbfilename, iostat=ostat, access="sequential", action="readwrite", status="old")

		do
			read(unit=11, fmt="( A )", iostat=fstat) line
			if (fstat /= 0) exit !condition for end of file
			read(line(1:6), fmt="( A )") rectype
			if (rectype == "ATOM  ") exit
		enddo

		do
			read(unit=11, fmt="( A )", iostat=fstat) line
			if (fstat /= 0) exit !condition for end of file
			read(line(1:6), fmt="( A )") rectype
			if (rectype == "TER   ") exit
			if (rectype /= "ATOM  ") cycle			
			read(line(22:22), fmt="( A )") chain
			if(chain /= chainid) exit
			read(line(13:16), fmt="( A )") atomtype
			if (atomtype == " CA ") then
				n = n + 1
				read(line(18:20), fmt="( A )") residue
				if(residue == "GLY ") then
					n = n - 1
				endif
			endif
		enddo
		restotal(filecount) = n
		close(unit=11, iostat=cstat, status="delete")
	end subroutine filecounter




	subroutine firstpass(pdbfilename, chainid, top, bottom, filemax, &
	totalsize, filenumber, filecount, arrayres, restotal, &
	abinsize, abinsizecounter)
	!takes in the pdb files, reads them for the xyz location of the ca(1,2)/cb(1,2) and passes them to bincounts to populate the array of bincounts to allocate arraydist
	!!!!!-------FIRSTPASS NEEDS TO WGET THE FILE AGAIN TO ACCESS---
		character(len=8), intent(in) :: pdbfilename
		character(len=1), intent(in) :: chainid
		integer, intent(in) :: top, bottom, filemax, totalsize, filenumber, filecount
		type (residuetype), intent(inout), allocatable :: arrayres(:)
		integer, intent(inout) :: restotal(filemax)
		integer, intent(inout) :: abinsize(top)
		integer, intent(inout) :: abinsizecounter(top)

		integer :: n, ostat, fstat, cstat, alloc
		character(len=4) :: atomtype
		character(len=6) :: rectype
		character(len=3) :: residue
		character(len=78) :: line	
		character(len=1) :: chain
		character(len=4) :: pdbid

		allocate(arrayres(restotal(filenumber)))
		pdbid = pdbfilename(1:4)
		n = 1

		call getpdbfile(pdbfilename)

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
			if(chain /= chainid) exit
			arrayres(n)%aminoacid = residue
			read(line(23:26), fmt="( I4 )") arrayres(n)%residuenumber
			read(line(31:54) , fmt="( 3F8.3 )") arrayres(n)%cacoordinate(1:3)
			arrayres(n)%pdbid = pdbid
			arrayres(n)%chainid = chain
			!FINDS THE CB COORDINATES AHEAD OF THE CURRENT LINE
			do
				read(unit=11, fmt="( A )", iostat=fstat) line
				read(line(13:16), fmt="( A )") atomtype
				if (atomtype == " CB ") then
					read(line(31:54) , fmt="( 3F8.3 )") arrayres(n)%cbcoordinate(1:3)
					exit
				endif
			enddo 
			n = n + 1 
		enddo

		call bincounts(top, bottom, filemax, totalsize, filenumber, &
		filecount, arrayres, restotal, abinsize, abinsizecounter)

		close(unit=11, iostat=cstat, status="delete")
		deallocate(arrayres)
	end subroutine firstpass

	subroutine secondpass(pdbfilename, chainid, top, bottom, filemax, &
	totalsize, filenumber, filecount, arrayres, arraydist, restotal, &
	abinsize, abinsizecounter)
	!takes in the pdb files,reads them for the xyz location of the ca(1,2)/cb(1,2) and passes them to computeseparations to populate arraydist and build the singlysorted file
	!!!!-------SECONDPASS NEEDS TO WGET THE FILE AS WELL, BUT IT IS THE LAST TIME THIS HAPPENS
		character(len=8), intent(in) :: pdbfilename
		character(len=1), intent(in) :: chainid
		integer, intent(in) :: top, bottom, filemax, totalsize, filenumber, filecount
		type (residuetype), intent(inout), allocatable :: arrayres(:)
		type (distancetype),  intent(inout), allocatable :: arraydist(:)
		integer, intent(inout) :: restotal(filemax)
		integer, intent(inout) :: abinsize(top)
		integer, intent(inout) :: abinsizecounter(top)

		integer :: n, ostat, fstat, cstat, alloc
		character(len=4) :: atomtype
		character(len=6) :: rectype
		character(len=3) :: residue
		character(len=78) :: line	
		character(len=1) :: chain
		character(len=4) :: pdbid

		allocate(arrayres(restotal(filenumber)))
		pdbid = pdbfilename(1:4)
		n = 1

		call getpdbfile(pdbfilename)
		open(unit=11, file=pdbfilename, iostat=ostat, access="sequential", &
		action="readwrite", status="old")

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
			if(chain /= chainid) exit
			arrayres(n)%aminoacid = residue
			read(line(23:26), fmt="( I4 )") arrayres(n)%residuenumber
			read(line(31:54) , fmt="( 3F8.3 )") arrayres(n)%cacoordinate(1:3)
			arrayres(n)%pdbid = pdbid
			arrayres(n)%chainid = chain
			!FINDS THE CB COORDINATES AHEAD OF THE CURRENT LINE
			do
				read(unit=11, fmt="( A )", iostat=fstat) line
				read(line(13:16), fmt="( A )") atomtype
				if (atomtype == " CB ") then
					read(line(31:54) , fmt="( 3F8.3 )") arrayres(n)%cbcoordinate(1:3)
					exit
				endif
			enddo 
			n = n + 1 
		enddo

		call computeseparations(top, bottom, filemax, totalsize, &
		filenumber, filecount, arrayres, arraydist, restotal, &
		abinsize, abinsizecounter)

		close(unit=11, iostat=cstat, status="delete")
		deallocate(arrayres)
	end subroutine secondpass

	subroutine bincounts(top, bottom, filemax, totalsize, filenumber, &
	filecount, arrayres, restotal, abinsize, abinsizecounter)
		!access the xyz coordinates for the ca(1,2) and compute the distance between them in 3d space. EACH SEPARATION IS CALCULATED USING THE DISTANCE FORUMULA FOR A 3 COORDINATE SYSTEM
		integer, intent(in) :: filemax
		integer, intent(in) :: filenumber, totalsize, filecount
		integer, intent(in) :: top, bottom
		integer, dimension(filemax), intent(inout) :: restotal
		type (residuetype), intent(in) :: arrayres(restotal(filenumber))
		integer, dimension(top), intent(inout) :: abinsize
		integer, dimension(top), intent(inout) :: abinsizecounter
		integer, dimension(top) :: temp !temporary bin counter
		integer :: cabin, cbbin, counter
		integer :: i, j
		integer :: ostat,cstat,fstat
		real :: distca, distcb

		counter = 0
		temp(1:top) = 0

		do i = 1, restotal(filenumber) - 2
			do j = i + 2, restotal(filenumber)
				!CACA
			  distca = sqrt(sum((arrayres(i)%cacoordinate(1:3)-arrayres(j)%cacoordinate(1:3))**2))
				cabin = ceiling((distca * 10) - 0.5)

				if(distca*10 > real(top)) cycle
		    if(distca*10 < real(bottom)) cycle
		    
				!CBCB
				distcb = sqrt(sum((arrayres(i)%cbcoordinate(1:3)-arrayres(j)%cbcoordinate(1:3))**2))
				cbbin = ceiling((distcb * 10) - 0.5)

				if(distcb*10 .gt. real(top)) cycle
				if(distcb*10 .lt. real(bottom)) cycle

				abinsize(cabin) = abinsize(cabin) + 1
				enddo
			enddo
		end subroutine bincounts

		subroutine computeseparations(top, bottom, filemax, totalsize, &
		filenumber, filecount, arrayres, arraydist, restotal, &
		abinsize, abinsizecounter)
			integer, intent(in) :: filemax
			integer, intent(in) :: filenumber, totalsize, filecount
			integer, intent(in) :: top, bottom
			integer, intent(inout) :: restotal(filecount)
			type (residuetype), intent(in) :: arrayres(restotal(filenumber))
			type (distancetype), intent(inout), allocatable :: arraydist(:)
			integer, dimension(top), intent(inout) :: abinsize
			integer, dimension(top), intent(inout) :: abinsizecounter
			integer, dimension(top) :: temp 																	!temporary bin counter
			integer :: cabin, cbbin, counter
			integer :: i, j
			integer :: ostat,cstat,fstat
			real :: distca, distcb

			counter = 0
			allocate(arraydist(restotal(filenumber)*((restotal(filenumber) - 1))/2))

			do i = 1, restotal(filenumber) - 2
				do j = i + 2, restotal(filenumber)
					!CACA
			  	distca = sqrt(sum((arrayres(i)%cacoordinate(1:3)-arrayres(j)%cacoordinate(1:3))**2))
					cabin = ceiling((distca * 10) - 0.5)

					if(distca*10 > real(top)) cycle
			    if(distca*10 < real(bottom)) cycle
		    
					!CBCB
					distcb = sqrt(sum((arrayres(i)%cbcoordinate(1:3)-arrayres(j)%cbcoordinate(1:3))**2))
					cbbin = ceiling((distcb * 10) - 0.5)
					if(distcb*10 .lt. real(bottom)) cycle
					if(distcb*10 .gt. real(top)) cycle

					counter = counter + 1
					arraydist(counter)%cbseparation = distcb
					arraydist(counter)%caseparation = distca
	
					!CACB
					arraydist(counter)%cacross = sqrt(sum((arrayres(i)%cacoordinate(1:3)-arrayres(j)%cbcoordinate(1:3))**2))
	
					!CBCA
					arraydist(counter)%cbcross = sqrt(sum((arrayres(i)%cbcoordinate(1:3)-arrayres(j)%cacoordinate(1:3))**2))
	
					!save the residue numbers (left and right)
					arraydist(counter)%pdbid = arrayres(i)%pdbid
					arraydist(counter)%chainid = arrayres(i)%chainid
					arraydist(counter)%rightresiduenumber = arrayres(j)%residuenumber
					arraydist(counter)%leftresiduenumber = arrayres(i)%residuenumber
					arraydist(counter)%rightanchor = arrayres(j)%aminoacid
					arraydist(counter)%leftanchor = arrayres(i)%aminoacid
					arraydist(counter)%cabin = cabin
					arraydist(counter)%cbbin = cbbin
					arraydist(counter)%length = arrayres(j)%residuenumber - arrayres(i)%residuenumber
					!WRITE THE CABIN NUMBER, CBBIN NUMBER, ATOM1NUMBER, ATOM2NUMBER, DISTANCES, AND PDBID ONLY LESS THAN 20Angstroms
				enddo				
			enddo
		  call cabinstructure(counter, totalsize, top, abinsize, arraydist, abinsizecounter)
			deallocate(arraydist)
	end subroutine computeseparations


	subroutine cabinstructure(counter, totalsize, top, abinsize, arraydist, abinsizecounter)
	!use this subroutine organizes the ca bin structure, filled with zeroes
		integer, intent(in) :: counter, top, totalsize
		integer, dimension(top), intent(inout) :: abinsize
		integer, dimension(top), intent(inout) :: abinsizecounter
		type (distancetype), intent(inout), dimension(totalsize) :: arraydist
		type (distancetype), dimension(totalsize) :: arrayc
		

		character(len=*), parameter :: infile="singlysorted.txt"		

		integer :: ostat,fstat,cstat, alloc
		integer :: ostat2,fstat2,cstat2
		integer :: lstat
		integer :: i, j, t, n
		integer :: cabin, cbbin

		t = 1 !give this counter a name

		!using the cabin structure, print to a file the records in bins organized by cabinstructure. include in the record the cabin number, cbbin number, cacb, cbca, anchoring residue by chain, orthogonal cordinates, and pdbid
		 
		!loop through the bins to find the locations where the beginning of the cabin record is

		open(unit=19, file=infile, form="formatted", iostat=ostat, access="direct", action="readwrite", status="old", recl=66)
		!loop through the array structure arraydist (arraydist), calculate the cabin from the casep and cbbin from cbsep
		do i=1, counter
			cabin = arraydist(i)%cabin
			t = sum(abinsize(1:cabin-1)) + abinsizecounter(cabin)
      write(unit=19, fmt="( 2(I4), 2(A4, I4), 4F8.3, A4, A1, I4)", iostat=fstat, rec=t) arraydist(i)
			abinsizecounter(cabin) = abinsizecounter(cabin) - 1
		enddo
		close(unit=19,iostat=cstat)
	end subroutine cabinstructure


	subroutine cbbinstructure(totalsize, top, abinsize)
    integer, intent(inout) :: abinsize(:)
    integer, intent(inout) :: totalsize
		integer, intent(in) :: top

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

    open(unit=19, file=infile, form="formatted", iostat=ostat, access="direct", action="read", status="old", recl=66)
    open(unit=39, file=outfile, form="formatted", iostat=ostat2, access="direct", action="readwrite", status="replace", recl=66)
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
        read(unit=19, fmt="( 2(I4), 2(A4, I4), 4F8.3, A4, A1, I4 )" , iostat = fstat, rec=j) line
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

			

      write(unit=37, fmt='( 2(I8),200(I8) )', iostat=fstat, rec=i) &
			width, linenum, bbinsizefinal

      do j = linenum, linenum + width - 1
        read(unit=19, fmt="( 2(I4), 2(A4, I4), 4F8.3, A4, A1, I4)" , &
				iostat=fstat, rec=j) line
        cbbin = line%cbbin
        t = sum(bbinsize(1:cbbin)) - bbinsizecounter(cbbin)
        t = t + linenum
        write(unit=39, fmt="( 2(I4), 2(A4, I4), 4F8.3, A4, A1, I4)", iostat=fstat2, rec=t) line
        bbinsizecounter(cbbin) = bbinsizecounter(cbbin) - 1
      enddo

      linenum = linenum + width
    enddo

    close(unit=37, iostat=cstat)
    close(unit=19, iostat=cstat)
    close(unit=39, iostat=cstat2)

  end subroutine cbbinstructure

!the cbbin values for each bin need to be saved, as does the cabinsize values in two separate files. save the cabin file as a cabin number, the cbbin number and the count. save the cabin file as it is structured, cabin and count

!WRITE IN THE COMMENTS FOR EACH VARIABLE SO TH
end module floss


program main
	use floss

	character(len=4) :: pdbid												!4 character pdb id
	character(len=1) :: chainid											!1 character chain id
	character(len=*), parameter :: txt = '.txt' 		!txt extension parameter
	real, parameter :: maximum = 20.0 							!maximum angstrom distance
	real, parameter :: resolution = 0.10						!detail
	integer, parameter :: top = int(maximum/resolution) 	!maximum number of calpha bins
	real, parameter :: minimum = 	4.0								!minimum angstrom distance
	integer, parameter :: bottom = int(minimum/resolution) !minimum number of allowable cabin
	integer, parameter :: filemax = 1000						!maximum number of files to be read
	character(len=200) :: pdbfilename, allpdbs			!pdbfilename: stores the name of the pdb file !capture the file name argument 
	integer :: ostat,cstat,fstat 										!file iostat carriers o:open, c:close, f:in use file
	integer :: argnum, indexnum 										!argnum: iargc carrier, the number of arguments given in the command line !
	integer, dimension(top) :: abinsize 						!array structure stores the cabin with its binsize
	integer, dimension(top) :: abinsizecounter			!counter for each cabin
	integer ::filecount, totalsize 									!filecount: total number of files opened, totalsize: total number of distances calculated
	integer :: filenumber														!filenumber for indexing through restotal array to define size of arrayres
	integer, dimension(filemax) :: restotal 				!array structure stores the total number of residues in each file in order of processing
	type (residuetype), dimension(:), allocatable :: arrayres		!array structure stores the residues of each pdbfile
	type (distancetype), dimension(:), allocatable :: arraydist !array stores the distances between each residue pair in a single pdb file

	!initialization of necessary variables
	restotal(1:filemax) = 0
	abinsize(1:top) = 0
	abinsizecounter = abinsize
	filecount = 0
	
	!process the commandline arguments
	argnum = iargc()
	if(argnum /= 1) stop 
	call getarg(1,allpdbs)
	indexnum = index(allpdbs,txt)
	if(indexnum == 0) then
		write(*,*) "the file", allpdbs, &
		" could not be opened, proceed 	with a txt file"
		stop
	endif

	!read through the command line file for each pdbid and chainid
	!!!!!!!REFORMAT TO READ THE SELECTFILE!!!!!!!!!!!!!
	open(unit=15, file=allpdbs, iostat=ostat, & 
	access="sequential", 	action="read", status="old")

	do i=1, filemax
		filecount = filecount + 1	
		read(unit=15,fmt="( A4,A1 )", iostat=fstat) pdbid, chainid
		pdbfilename = pdbid//'.pdb' !ADD IN A PATHWAY SPECIFIER

		if(fstat /= 0) then
			filecount = filecount - 1
			exit
		endif

		call filecounter(pdbfilename, chainid, top, bottom, filemax, filecount, restotal, abinsize)
		write(*,*) "file ", filecount, " completed, with ", &
		restotal(filecount), "residues"
	enddo
	

	totalsize = sum((restotal(1:filecount)*(restotal(1:filecount) - 1)/2))

	rewind 15

	do	i=1,filecount
		filenumber = i
		read(unit=15,fmt="( A4,A1 )", iostat=fstat) pdbid, chainid
		pdbfilename = pdbid//'.pdb'
		call firstpass(pdbfilename, chainid, top, bottom, filemax, &
		totalsize, filenumber, filecount, arrayres, restotal, &
		abinsize, abinsizecounter)
	enddo

	abinsizecounter = abinsize
	rewind 15

	do	i=1,filecount
		filenumber = i				
		read(unit=15,fmt="( A4,A1 )", iostat=fstat) pdbid, chainid
		pdbfilename = pdbid//'.pdb'
		call secondpass(pdbfilename, chainid, top, bottom, filemax, &
		totalsize, filenumber, filecount, arrayres, arraydist, restotal, &
		abinsize, abinsizecounter)
		write(*,*) totalsize, top, sum(abinsizecounter(1:top)), sum(abinsize(1:top))
	enddo

	call cbbinstructure(totalsize, top, abinsize)
	close(unit=15, iostat=cstat)
!----WRITE A CHECK TO IDENTIFY WHEN BINS ARE FULL. IF ALL BINS ARE FULL, EXIT PROGRAM---!
!---ERASE DOWNLOADED FILES AFTER EACH SUBROUTINE FINISHES---!
end program main









					
