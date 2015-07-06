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
!
program main
	character(len=*), parameter :: filename="diracloop.txt"
	character(len=200) :: pdbfilename, allpdbs
	!1 = pdbfilename, 2=diracloop.txt, 3=allpdbs, 4=prelim
	integer :: ostat,cstat,fstat
	integer :: argnum, indexnum, linenum, cabin, cbbin
	integer :: bincount, totalsize, totalres

	!first element in abinsize is cabinsize, second is cbbinsize
	integer, dimension(1000) :: abinsize, bbinsize, binloc
	integer, parameter :: r8=selected_real_kind(4,7)
	!an allocatable array to store the information
	integer, allocatable :: bins(:)
	real(kind=r8), dimension(1,3):: arraya, arrayb !make allocatable to correspond to the variable number of residues in each file
	real(kind=r8), dimension(600,600,4):: arrayd

	bincount = 0
	totalsize = 0
	linenum = 1

	!commenting here represents the switch from reading in the pdb fileindividually to reading from a file with a list of pdb file names
	argnum = iargc()
	if(argnum /= 1) stop 
	call getarg(1,allpdbs)
	!call getarg(1,pdbfilename)
	indexnum = index(allpdbs,'.txt')
	!indexnum = index(pdbfilename,'.pdb')
	if(indexnum == 0) then
		write(*,*) "the file", allpdbs, " could not be opened, proceed with a txt file"
		stop
	endif

	!do for all pdb files in text list of pdb filenames
	
	open(unit=15, file=allpdbs, iostat=ostat, access="sequential", action="read", status="old")
	!FIRST PASS THROUGH, THIS PASS ORGANIZES THE STRUCTURE FOR THE ARRAY BEFORE WRITING.
	do
		read(unit=15,fmt="( A )", iostat=fstat) pdbfilename
		!write(*,*) pdbfilename
		if(fstat /= 0) exit
		call readinpoints(pdbfilename)
		call computeseparations(totalres)
	enddo
	rewind 15
	call cabinstructure
	close(unit=15, iostat=cstat)
	contains
	
	subroutine readinpoints(filename)
	!readinpoints takes in the pdb files,reads them for the xyz location of the ca(1,2)/cb(1,2) and passes them to computeseparations to determine the values for each index
		character(len=8), intent(in) ::filename
		integer :: n, ostat, fstat, cstat
		character(len=4) :: atomtype
		character(len=6) :: rectype
		character(len=3) :: residue
		character(len=78) :: line	

		n = 1
		totalres = 0

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
			read(line(31:54) , fmt="( 3F8.3 )") arraya(n,1:3)
			n = n + 1 !that is so dumb, make n the residue number, get residue number and the chain ID because it may not start from 1
		enddo

		totalres = n
		rewind 11 !go back to the top of the file		
		n = 1

		do 
			!for all atoms in file | read the first anchor
			read(unit=11, fmt="( A )", iostat=fstat) line
			if (fstat /= 0) exit !condition for end of file
			read(line(1:4), fmt="( A )") rectype
			if (rectype /= "ATOM  ") cycle
			read(line(13:16), fmt="( A )") atomtype
			if (atomtype /= " CB") cycle !check for CA
			if(residue == "GLY ") cycle !ignore the glycines
			read(line(31:54) , fmt="( 3F8.3 )") arrayb(n,1:3)
			n = n + 1
		enddo
	
		close(unit=11, iostat=cstat)	
	end subroutine readinpoints
	






	subroutine computeseparations(totalres)
		!rewrite as dists 
		!access the xyz coordinates for the ca(1,2) and compute the distance between them in 3d space.EACH SEPARATION IS CALCULATED USING THE DISTANCE FORUMULA FOR A 3 COORDINATE SYSTEM
		integer :: i,j, num
		integer, intent(in) :: totalres
		integer :: ostat,cstat,fstat, numstat
		character(len=4) :: pdbid

		pdbid = pdbfilename(1:4)

		open(unit=29, file="prelim.txt", form="formatted", iostat=ostat, access="sequential", action="readwrite", status="replace")
		!this wont work if there is no available cb, determine for glycines
		
		do i=1, totalres
			do j=1, i-2
				if (i == j) cycle
				!CACA
				!arrayd(i,j,1) = sqrt(((arraya(i,1)-arraya(j,1))**2)+((arraya(i,2)-arraya(j,2))**2)+((arraya(i,3)-arraya(j,3))**2))
        arrayd(i,j,1) = sqrt(sum((arraya(i,1:3)-arraya(j,1:3))**2))!bin (by .1A) !IMPLIED LOOP
				cabin = int(arrayd(i,j,1)*10)

				!CBCB
				arrayd(i,j,2) = sqrt(((arrayb(i,1)-arrayb(j,1))**2)+((arrayb(i,2)-arrayb(j,2))**2)+((arrayb(i,3)-arrayb(j,3))**2))
				!which bin (by .1A)
				cbbin = int(arrayd(i,j,2)*10)

				!CACB
				arrayd(i,j,3) = sqrt(((arraya(i,1)-arrayb(j,1))**2)+((arraya(i,2)-arrayb(j,2))**2)+((arraya(i,3)-arrayb(j,3))**2))

				!CBCA
				arrayd(i,j,4) = sqrt(((arrayb(i,1)-arraya(j,1))**2)+((arrayb(i,2)-arraya(j,2))**2)+((arrayb(i,3)-arraya(j,3))**2))

				!WRITE THE CABIN NUMBER, CBBIN NUMBER, ATOM1NUMBER, ATOM2NUMBER, DISTANCES, AND PDBID ONLY LESS THAN 50Angstroms
				if(cabin<500 .and. cabin > 18) then
					if(cbbin<500 .and. cbbin > 18) then
						call countbins(cabin,cbbin)
						write(unit=29, fmt='( I4, I4, I4, I4, 4F8.3, A5 )', iostat=fstat, advance="yes") cabin, cbbin, i, j, arrayd(i,j,1:4), pdbid
					endif
				endif	
			enddo
		enddo
	end subroutine computeseparations





	subroutine countbins(cabin, cbbin)
		integer, intent(in) ::  cabin, cbbin
	!count bins will populate an array which holds integers corresponding to the number of elements which each bin will need to have capacity for
		abinsize(cabin) = abinsize(cabin) + 1 
		bbinsize(cbbin) = bbinsize(cbbin) + 1
		totalsize = totalsize + 1
		!write(*,*) cabin, abinsize(cabin), cbbin, bbinsize(cbbin), totalsize
	end subroutine countbins
	
subroutine cabinstructure
	!use this subroutine organizes the ca bin structure, filled with zeroes
		character(len=53) :: line
		integer :: ostat,fstat,cstat
		integer :: ostat2,fstat2,cstat2
		integer :: lstat
		integer :: anchor1,anchor2
		integer :: i, j, t, last
		real(kind=r8), dimension(4) :: dist
		
		last = 0
		t = 1 !give this counter a name

		allocate (bins(totalsize)) 
		!populate array with default value of 0
		do i=1,size(bins,1)
			bins(i) = 0
		enddo		

	!using the cabin structure, print to a file the records in bins organized by cabinstructure. include in the record the cabin number, cbbin number, cacb, cbca, anchoring residue by chain, orthogonal cordinates, and pdbid
		 
		!loop through the bins to find the locations where the beginning of the cabin record is
		open(unit=19, file=filename, form="formatted", iostat=ostat, access="direct", action="readwrite", status="replace", recl=57)

		rewind 29

		do i=1, size(abinsize,1)
			if(abinsize(i) == 0) cycle
			bins(t:t+abinsize(i) - 1) = i
			t = t + abinsize(i) - 1
		enddo
		
		do
			read(unit=29, fmt='( I4, I4, I4, I4, 4F8.3, A5 )', iostat=fstat2) cabin, cbbin, anchor1, anchor2, dist(1:4), pdbid
			if(fstat2 /= 0) exit
			do j=1,size(bins,1)
				if(cabin == bins(j)) then !only works for exct matches	
					last = j + abinsize(cabin) - 1
					write(unit=19, fmt='( I4, I4, I4, I4, 4F8.3, A5 )', iostat=fstat, rec=last) cabin, cbbin, anchor1, anchor2,dist(1:4), pdbid
					abinsize(cabin) = abinsize(cabin) - 1
					write(*,*) cabin, abinsize(cabin), last, totalsize
					exit
				endif
			enddo
			
		enddo
		close(unit=29, iostat=cstat2)
		close(unit=19, iostat=cstat)
	end subroutine cabinstructure

!the next thing you need to do is to set up the cbbinstructure the same way you did with the cabinstructure. repeat the last two subroutines but for the second filter

!save the chain id information
!that is so dumb, make n the residue number, get residue number and the chain ID (see following comment)
!rewrite captured info as defined types residuetype(int resnumber, char chainid, char pdbid, char aa, real(3) cacoord,cbcoord) and dist(int resnum1, int resnum2, real caca,cbcb,cacb,cbca, char pdbid, char chainid) as globals (construct, declare as allocatable)
!rewrite the initial allpdb file as a file of pdbid and chain Id to be subsequently written to a file of pathways (use trim and string concatenation)
end program main
