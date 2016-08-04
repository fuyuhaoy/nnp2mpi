! 2016-8-3
! fix some bugs

program maketrain
  implicit none
  include '/opt/openmpi-1.6.5/include/mpif.h'
  !include '/opt/intel/impi/4.1.3.049/include64/mpif.h'
  !include '/opt/intel/impi/4.1.0.024/include64/mpif.h'
  
  integer i,j,k
  integer nat
  parameter (nat=136) ! max number of atoms in structure for dimensioning of arrays
  integer nn ! number of input nodes
  parameter (nn=56)
  integer counter
  integer ncut
  integer ntest,ntrain,nrej
  integer natom ! real number of atoms in structure
  integer nexample
  integer ntot

  real*8 lattice(3,3)
  real*8 xyz(nat,3)
  real*8,  dimension(:,:), allocatable :: xxyyzz 
  !real*8 xxyyzz(27*nat,3)
  real*8 energy,e,diffe,de,ej
  real*8 zran
  real*8 cutoff
  real*8 xx(nn,1,nat)
  real*8 dxdy(nn,nat,nat,3) ! derivatives of symmetry functions with respect to x,y and z
  real*8 strs(3,3,nn,nat)
  real*8 bohr2a,a2bohr
  parameter (bohr2a=0.529177d0)
  parameter (a2bohr=1.8897268d0)
  real*8 ethres 
  real*8 volume, volumemax
  real*8 tempvec(3) 
  
  character*15 cdummy
  
  logical ldebug
  logical lforce
  !lforce=.true.
  
  ! --------------------fyuhaoy--------------------
  integer :: ios ! status when read some files

  ! related with mpi
  integer :: myid, numprocs, structures_per_procs
  integer :: ierr, status(MPI_STATUS_SIZE), sender
  
  integer :: quotient, remainder    ! used for calculation of counts array
  integer, allocatable :: counts(:) ! number of structure for evergy process
  integer, dimension(2) :: displs   ! range of structure for every process
  
  integer :: structures ! total number of structures
  character(len=15) :: fdebug, ftrain, ftest, fmonitor
  character(len=1000) :: fdebug2, ftrain2, ftest2, fmonitor2 ! used for string concatenation

  ! judge whether process is finished
  logical, allocatable :: pstatus(:)
  logical :: subpstatus
  
  logical :: tmps ! monitor status of mpi message transmission

  integer, allocatable :: ntrain2(:), ntest2(:), nrej2(:) ! for all porcesses
  character(len=1024) :: buffer ! buffer of IO for train.data and test.data files
  integer :: size ! actual size of buffer during reading
  
  ! ------------------end fyuhaoy------------------
  
  ! initializations
  counter=0
  !ldebug=.true.
  ldebug=.false.
  ncut=1
  ntrain=0
  ntest=0
  nrej=0
  !ethres=-7.7960 ! energy per atom in Ry = 498.944 for 64 
  !ethres=-7.8125 ! energy per atom in Ry = 500 for 64
  !ethres=-7.765 ! energy per atom in Ry = 497 for 64
  ethres=0.0 !  in fact no energy threshold 
  !ethres=-106.290 ! energy per atom in eV
  !ethres=-500.d0 !for 64 atom cell in Ry
  !volumemax=140.d0 ! Bohr3 per atom, caution with clusters in large boxes
  volumemax=400.d0 ! in fact no volume threshold 
  natom=0
  nexample=0
  ntot=0

  ! --------------------fyuhaoy--------------------
  ios=0
  structures=0
  quotient=0
  remainder=0
  displs=0
  fdebug2=""
  ftrain2=""
  ftest2=""
  fmonitor2=""
  buffer=""
  size=0
  
  call MPI_INIT (ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD, myid, ierr)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, numprocs, ierr)

  allocate(pstatus(numprocs))
  pstatus=.false.
  subpstatus=.false.

  allocate(ntrain2(numprocs))
  allocate(ntest2(numprocs))
  allocate(nrej2(numprocs))
  ntrain2=0
  ntest2=0
  nrej2=0
  
  if (myid .eq. 0) then
     ! count atoms in structures firstly
     open(20,file='input.data',form='formatted',status='old')
     open(21,file='temp.data',form='formatted',status='replace')
     
     ! count nunber of atoms in every structure, and output to temp.data file
     do while (ios .eq. 0)
        read(20,*,IOSTAT=ios) cdummy
        if (ios .eq. 0) then
           nexample=nexample+1 ! counter of lines
           if (cdummy(1:4) .eq. 'Next') then
              structures=structures+1
              if (nexample .gt. 1) then
                 write(21,*) natom-4
              end if
              natom=0
           else
              natom=natom+1
           end if
        end if
     end do
     write(21,*) natom-4 ! the last structure

     close(20)
     close(21)
     ios=0 ! the status recover to 0

     ! check whether the number of atoms are beyond the upper limit;
     ! search the maximum value of atoms in all structures
     open(21,file='temp.data',form='formatted',status='old')
     do while (ios .eq. 0)
        read(20,*,IOSTAT=ios) natom
        if (ios .eq. 0) then
           if (natom .gt. nat) then
              write(*,*)'Redimension nat!',natom
              stop
           end if
        end if
     end do
     
     close(21)
     ios=0 ! the status recover to 0
   
  end if ! myid=0
  
  call MPI_Bcast(structures, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
  !call MPI_Barrier(MPI_COMM_WORLD, ierr)

  ! ensure the array of distribution
  quotient=structures/numprocs
  remainder=mod(structures,numprocs)

  allocate(counts(numprocs))
  do i=1,numprocs-1
     counts(i)=quotient
  end do
  counts(numprocs)=quotient+remainder

  if (myid .eq. 0) then
     displs=[1,counts(myid+1)]
  else
     displs=[sum(counts(1:myid))+1,sum(counts(1:myid+1))]
  end if
  
  call structure2process(numprocs, myid, structures, counts, displs)

  ! parallel calculation of symmetry function 
  open(20,file='input.data',form='formatted',status='old')
  open(24,file='temp.data',form='formatted',status='old')

  ! every process hold its own file. i.e. debug.data.xx (xx=myid)
  if (myid .lt. 10) then
     write(fdebug, '(a11,i1)') "debug.data.", myid
     write(ftrain, '(a11,i1)') "train.data.", myid
     write(ftest, '(a10,i1)') "test.data.", myid
     write(fmonitor, '(a13,i1)') "monitor.data.", myid
  elseif (myid .lt. 100) then
     write(fdebug, '(a11,i2)') "debug.data.", myid
     write(ftrain, '(a11,i2)') "train.data.", myid
     write(ftest, '(a10,i2)') "test.data.", myid
     write(fmonitor, '(a13,i2)') "monitor.data.", myid
  elseif (myid .lt. 1000) then
     write(fdebug, '(a11,i3)') "debug.data.", myid
     write(ftrain, '(a11,i3)') "train.data.", myid
     write(ftest, '(a10,i3)') "test.data.", myid
     write(fmonitor, '(a13,i3)') "monitor.data.", myid
  else
     write(*,*) "Warning! too much processes(>=1000). The program will exit."
     stop
  end if
     
  open(22,file=fdebug,form='formatted',status='replace')
  open(21,file=ftrain,form='formatted',status='replace')
  open(23,file=ftest,form='formatted',status='replace')
  open(33,file=fmonitor)

  ntot=0
  counter=0
  ! start to calculation of symmetry function
  do while (ios .eq. 0)
     read(20,*,IOSTAT=ios) cdummy
     read(24,*,IOSTAT=ios) natom
     
     if (ios .eq. 0) then
        if (cdummy(1:4).eq.'Next') then
           ntot=ntot+1       ! total structure
           counter=counter+1 ! vaild structure
        end if

        if ((ntot .ge. displs(1)) .and. (ntot .le. displs(2))) then ! only address these process own's structures
           do i=1,3
              read(20,*) lattice(i,1),lattice(i,2),lattice(i,3)
           end do
           do i=1,natom
              read(20,*) xyz(i,1),xyz(i,2),xyz(i,3)
           end do
           read(20,*)energy
           ! we have to use the energy per atom if systems of different size are used
           ! to have larger numbers, we use eV
           energy=energy/dble(natom) ! still in Ry
           ! calculation of the volume
           volume=0.0d0
           tempvec(1)=lattice(1,2)*lattice(2,3)-lattice(1,3)*lattice(2,2)
           tempvec(2)=lattice(1,3)*lattice(2,1)-lattice(1,1)*lattice(2,3)
           tempvec(3)=lattice(1,1)*lattice(2,2)-lattice(1,2)*lattice(2,1)
           volume=tempvec(1)*lattice(3,1)+tempvec(2)*lattice(3,2)&
                +tempvec(3)*lattice(3,3)
           volume=volume/dble(natom)

           if (ldebug) then
              write(22,*)'lattice vectors in Bohr'
              do i=1,3
                 write(22,'(3f14.6)')lattice(i,1),lattice(i,2),lattice(i,3)
              enddo
              write(22,*)'Atomic positions in Bohr'
              do i=1,natom
                 write(22,'(3f14.8)')xyz(i,1),xyz(i,2),xyz(i,3)
              enddo
           endif

           if ((energy .gt. ethres) .or. (volume .gt. volumemax)) then
              if (energy .gt. ethres) then
                 write(*,'(a,i6,a,i6,a,f14.8,a,i6)') &
                      'Point ',ntot,' : ',natom,' atoms, E=',energy,' rejected '
              end if
              if (volume .gt. volumemax) then
                 write(*,'(a,i6,a,i6,a,f14.8,a,i6)') &
                      'Point ',ntot,' : ',natom,' atoms, V=',volume,' rejected '
              end if
              
              counter=counter-1
              nrej=nrej+1
           else ! vaild structure
              ! convert from Bohr to Angstrom
              do i=1,3
                 lattice(i,1)=lattice(i,1)*bohr2a
                 lattice(i,2)=lattice(i,2)*bohr2a
                 lattice(i,3)=lattice(i,3)*bohr2a
              end do
              do i=1,natom
                 xyz(i,1)=xyz(i,1)*bohr2a
                 xyz(i,2)=xyz(i,2)*bohr2a
                 xyz(i,3)=xyz(i,3)*bohr2a
              end do
              if (ldebug) then
                 write(22,*)'lattice vectors in Angstrom'
                 do i=1,3
                    write(22,'(3f14.6)')lattice(i,1),lattice(i,2),lattice(i,3)
                 end do
                 write(22,*)'Atomic positions in Angstrom'
                 do i=1,natom
                    write(22,'(3f14.8)')xyz(i,1),xyz(i,2),xyz(i,3)
                 end do
              endif

              allocate (xxyyzz(natom,3))
              do i=1,natom
                 xxyyzz(i,1)=xyz(i,1)
                 xxyyzz(i,2)=xyz(i,2)
                 xxyyzz(i,3)=xyz(i,3)
              end do
              ! one determination of all function values per geometry
              ! functions for all atoms are calculated at once
              ! write(*,*) lattice
              call functions(counter,natom,nn,lattice,xxyyzz,&
                   xx,dxdy,strs,.false.,.true.,ldebug,e,volume)
              deallocate (xxyyzz)
              
              ! check if write format is sufficient to separate numbers
              do i=1,natom
                 do j=1,9
                    if (xx(j,1,i) .gt. 10000.d0) then
                       write(*,*)'Warning too large xx ',j,i,xx(j,1,i)
                       stop
                    elseif (xx(j,1,i) .lt. -9999.d0) then
                       write(*,*) 'Warning too low xx ',j,i,xx(j,1,i)
                       stop
                    end if
                 end do
              end do
              
              call getran(zran)

              de=-14.9817443422860
              ej=e/(2*natom*1312.22) 
      
              if (zran .gt. 0.2d0)then ! train.data
              !if (.true.)then ! train.data
                 ntrain=ntrain+1
                 write(21,*)natom
                 diffe=energy-de-ej
                 do i=1,natom
                    write(21,'(600f18.11)') (xx(j,1,i),j=1,nn)
                 end do
                 write(21,'(f14.8,x,a)') energy,"    myid"
                 write(33,'(i8,x,2f20.12,i8)') ntot,ej,diffe,myid
              else ! test.data
                 diffe=energy-de-ej
                 ntest=ntest+1
                 write(23,*) natom
                 do i=1,natom
                    write(23,'(600f18.11)') (xx(j,1,i),j=1,nn)
                 end do
                 write(23,'(f14.8,x,a)') energy,"    myid"
                 write(33,'(i8,x,2f20.12,i8)') ntot,ej,diffe,myid
              endif
              
           end if ! vaild structure
        else
           do i=1,natom+4
              read(20,*)
           end do
        end if ! process's own structure
     end if ! status of reading
  end do ! read every line

  close(20)
  close(21)
  close(22)
  close(23)
  close(24)
  ios=0 ! the status recover to 0
  
  if ((ntrain+ntest) .ne. counts(myid+1)) then
     if (myid .lt. 10) then
        write(*,'(a, i1, a)')'Error, ntest+ntrain.ne.counter in ', myid, " process."
     elseif (myid .lt. 100) then
        write(*,'(a, i2, a)')'Error, ntest+ntrain.ne.counter in ', myid, " process."
     elseif (myid .lt. 100) then
        write(*,'(a, i3, a)')'Error, ntest+ntrain.ne.counter in ', myid, " process."
     else
        write(*,*) "Warning! too much processes(>=1000). The program will exit."
     end if
  endif

  subpstatus=.true.
  call MPI_Gather(subpstatus, 1, MPI_LOGICAL, pstatus, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  
  write(*,'(a, i3, a)')'----------------process ', myid, '----------------'
  write(6,*)counts(myid+1),' structures read from input.data'
  write(6,*)ntrain,' training points'
  write(6,*)ntest,' test points'
  write(6,*)nrej ,' rejected high energy points'

  call MPI_Gather(ntrain, 1, MPI_INT, ntrain2, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
  call MPI_Gather(ntest, 1, MPI_INT, ntest2, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
  call MPI_Gather(nrej, 1, MPI_INT, nrej2, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
  
  if (myid .eq. 0) then
     
     do i=0,numprocs-1

        if (i .lt. 10) then
           write(fdebug, '(a11,i1)') "debug.data.", i
           write(ftrain, '(a11,i1)') "train.data.", i
           write(ftest, '(a10,i1)') "test.data.", i
           write(fmonitor, '(a13,i1)') "monitor.data.", i
        elseif (i .lt. 100) then
           write(fdebug, '(a11,i2)') "debug.data.", i
           write(ftrain, '(a11,i2)') "train.data.", i
           write(ftest, '(a10,i2)') "test.data.", i
           write(fmonitor, '(a13,i2)') "monitor.data.", i
        elseif (i .lt. 1000) then
           write(fdebug, '(a11,i3)') "debug.data.", i
           write(ftrain, '(a11,i3)') "train.data.", i
           write(ftest, '(a10,i3)') "test.data.", i
           write(fmonitor, '(a13,i3)') "monitor.data.", i
        else
           write(*,*) "Warning! too much processes(>=1000). The program will exit."
           stop
        end if

        fdebug2=trim(fdebug2)//" "//fdebug
        ftrain2=trim(ftrain2)//" "//ftrain
        ftest2=trim(ftest2)//" "//ftest
        fmonitor2=trim(fmonitor2)//" "//fmonitor
        
     end do

     ! ensure all processes are finished
     tmps=.false.
     do while(.not.tmps)
        do i=1,numprocs
           if (pstatus(i)) then
              tmps=.true.
           else
              tmps=.false.
           end if
        end do
        call sleep(3) ! sleep 3s
     end do

     ! linux shell command
     call system("cat "//trim(fdebug2)//" > debug.data")
     call system("cat "//trim(ftrain2)//" > train.data")
     call system("cat "//trim(ftest2)//" > test.data")
     call system("cat "//trim(fmonitor2)//" > monitor.data")

     ! remove temporary file i.e. debug.data.xx
     call system("rm debug.data.*")
     call system("rm train.data.*")
     call system("rm test.data.*")
     call system("rm monitor.data.*")
     
     ! method 1
     ! train.data
     !if (sum(ntrain2) .gt. 0) then
     !   do i=1,sum(ntrain2)
     !      if (i .lt. 10) then
     !         write(string, '(a23,i1,a13)') "sed -i '0,/myid/s/myid/", i, "/' train.data"
     !      elseif (i .lt. 100) then
     !         write(string, '(a23,i2,a13)') "sed -i '0,/myid/s/myid/", i, "/' train.data"
     !      elseif (i .lt. 1000) then
     !         write(string, '(a23,i3,a13)') "sed -i '0,/myid/s/myid/", i, "/' train.data"
     !      else
     !         write(*,*) "Warning! too much processes(>=1000). The program will exit."
     !         stop
     !      end if
     !      call system(trim(string))
     !   end do
     !end if
     ! test.data
     !if (sum(ntest2) .gt. 0) then
     !   do i=1,sum(ntest2)
     !      if (i .lt. 10) then
     !         write(string, '(a23,i1,a12)') "sed -i '0,/myid/s/myid/", i, "/' test.data"
     !      elseif (i .lt. 100) then
     !         write(string, '(a23,i2,a12)') "sed -i '0,/myid/s/myid/", i, "/' test.data"
     !      elseif (i .lt. 1000) then
     !         write(string, '(a23,i3,a12)') "sed -i '0,/myid/s/myid/", i, "/' test.data"
     !      else
     !         write(*,*) "Warning! too much processes(>=1000). The program will exit."
     !         stop
     !      end if
     !      call system(trim(string))
     !   end do
     !end if

     ! method 2
     ! send pstatus to task 1, its aim is to for replacing (myid -> xx) test.data on this task.
     if (numprocs .gt. 1) then
        call MPI_Send(pstatus, numprocs, MPI_LOGICAL, 1, 2000, MPI_COMM_WORLD, ierr)
     end if
     
     ! replace myid to xx in train.data
     call system("mv train.data tmptrain.data")
     
     ! train.data
     open(200,file='tmptrain.data',form='formatted',status='old')
     open(201,file='train.data',form='formatted',status='replace')
     !open(204, file='monitor.data',form='formatted',position='append')
     
     ! Note that the status: 0:ok ; -1:end-of-file; -2:end-of-line
     counter=0
     do while (.not. IS_IOSTAT_END(ios))
        read(200,'(a)',ADVANCE='NO',IOSTAT=ios,SIZE=size)buffer
        if (ios .eq. 0) then
           write(201,'(a)',ADVANCE='NO')buffer(1:size)
        else if (ios .eq. -2) then ! end-of-line
           if (buffer(size-3:size) .eq. "myid") then
              counter=counter+1
              write(201,'(a,i8)',ADVANCE='NO')buffer(1:size-4),counter
           else
              write(201,'(a)',ADVANCE='NO')buffer(1:size)
           end if
           write(201,*) ! linefeed
        else if (ios .ne. -1) then
           write(*,'(a30,i4,i8,a27)')"Error! unknown status. IOSTAT=",ios,counter,"th structure in train.data."
           write(204,'(a30,i4,i8,a27)')"Error! unknown status. IOSTAT=",ios,counter,"th structure in train.data."
           stop ! exit program
        end if
     end do
     
     close(200)
     close(201)
     ios=0 ! the status recover to 0
     call system("rm tmptrain.data")
     
     ! replace myid to xx in test.data
     if (numprocs .eq. 1) then ! run the code when only using one process.
        call system("mv test.data tmptest.data")
        
        open(202,file='tmptest.data',form='formatted',status='old')
        open(203,file='test.data',form='formatted',status='replace')
        ! test.data
        counter=0
        do while (.not. IS_IOSTAT_END(ios))
           read(202,'(a)',ADVANCE='NO',IOSTAT=ios,SIZE=size)buffer
           if (ios .eq. 0) then
              write(203,'(a)',ADVANCE='NO')buffer(1:size)
           else if (ios .eq. -2) then ! end-of-line
              if (buffer(size-3:size) .eq. "myid") then
                 counter=counter+1
                 write(203,'(a,i8)',ADVANCE='NO')buffer(1:size-4),counter
              else
                 write(203,'(a)',ADVANCE='NO')buffer(1:size)
              end if
              write(203,*) ! linefeed
           else if (ios .ne. -1) then
              write(*,'(a30,i4,i8,a26)')"Error! unknown status. IOSTAT=",ios,counter,"th structure in test.data."
              write(204,'(a30,i4,i8,a26)')"Error! unknown status. IOSTAT=",ios,counter,"th structure in test.data."
              stop ! exit program          
           end if
        end do
        
        close(202)
        close(203)
        ios=0 ! the status recover to 0
        call system("rm tmptest.data")
     end if
     
     close(204)
     
     write(*,*)'>>>>>>>>>>>>>>>>>>>> total >>>>>>>>>>>>>>>>>>>>'
     write(6,*)sum(counts),' structures read from input.data'
     write(6,*)sum(ntrain2),' training points'
     write(6,*)sum(ntest2),' test points'
     write(6,*)sum(nrej2) ,' rejected high energy points'

  end if

  ! replace myid to xx in test.data
  if (myid .eq. 1) then
     call MPI_Recv(pstatus, numprocs, MPI_LOGICAL, 0, 2000, MPI_COMM_WORLD, status, ierr)

     ! ensure all processes are finished
     tmps=.false.
     do while(.not.tmps)
        do i=1,numprocs
           if (pstatus(i)) then
              tmps=.true.
           else
              tmps=.false.
           end if
        end do
        call sleep(3) ! sleep 3s
     end do
     
     call system("mv test.data tmptest.data")
     
     open(202,file='tmptest.data',form='formatted',status='old')
     open(203,file='test.data',form='formatted',status='replace')
     !open(205, file='monitor.data',form='formatted',position='append')
     ! test.data
     counter=0
     do while (.not. IS_IOSTAT_END(ios))
        read(202,'(a)',ADVANCE='NO',IOSTAT=ios,SIZE=size)buffer
        if (ios .eq. 0) then
           write(203,'(a)',ADVANCE='NO')buffer(1:size)
        else if (ios .eq. -2) then ! end-of-line
           if (buffer(size-3:size) .eq. "myid") then
              counter=counter+1
              write(203,'(a,i8)',ADVANCE='NO')buffer(1:size-4),counter
           else
              write(203,'(a)',ADVANCE='NO')buffer(1:size)
           end if
           write(203,*) ! linefeed
        else if (ios .ne. -1) then
           write(*,'(a30,i4,i8,a26)')"Error! unknown status. IOSTAT=",ios,counter,"th structure in test.data."
           write(205,'(a30,i4,i8,a26)')"Error! unknown status. IOSTAT=",ios,counter,"th structure in test.data."
           stop ! exit program          
        end if
     end do
     
     close(202)
     close(203)
     ios=0 ! the status recover to 0
     call system("rm tmptest.data")
     
     close(205)
     
  end if
  
  call MPI_FINALIZE(ierr)
  stop  
end program maketrain

!********************************************************
! calculate the range of structures by considering the number of atoms.
! numprocs: number of processes
! myid: id of process
! structures: number of structures
! displs: range of structures for every process
subroutine structure2process(numprocs, myid, structures, counts, displs)
  integer, intent(in) :: numprocs, myid, structures
  integer, intent(out) :: counts(numprocs), displs(2)
  integer :: atom(structures), tnatom
  integer :: tmpdispls(numprocs,2)
  integer :: i, p, natom
  counts=0
  displs=0
  atom=0
  natom=0
  p=0
  tmpdispls=0

  open(21,file='temp.data',form='formatted',status='old')
  do i=1,structures
     read(21,*) atom(i)
  end do
  close(21)
  
  ! total number of atoms in all structures
  tnatom=sum(atom(:))
  ! ensure the array of distribution
  quotient=tnatom/numprocs
  
  !do i=1,structures
  !   natom=natom+atom(i)
  !   if ((natom .le. quotient) .and. &
  !        (natom+atom(i+1) .gt. quotient)) then ! reaches the boundary
  !      p=p+1
  !      if (p .eq. 1) then
  !         tmpdispls(p,:)=[1,i]
  !      elseif (p .lt. numprocs) then
  !         tmpdispls(p,:)=[tmpdispls(p-1,2)+1,i]
  !      end if
  !      natom=0
  !      counts(p)=tmpdispls(p,2)+1-tmpdispls(p,1)
  !   end if
  !end do
  !tmpdispls(numprocs,:)=[tmpdispls(numprocs-1,2)+1,structures]
  !counts(numprocs)=tmpdispls(numprocs,2)+1-tmpdispls(numprocs,1)
  !displs=tmpdispls(myid+1,:)

  do i=1,structures-1
     natom=natom+atom(i)
     if ((natom .le. quotient) .and. &
          (natom+atom(i+1) .gt. quotient)) then ! reaches the boundary
        p=p+1
        if (p .eq. 1) then
           tmpdispls(p,:)=[1,i]
           counts(p)=tmpdispls(p,2)+1-tmpdispls(p,1)
        elseif (p .lt. numprocs) then
           tmpdispls(p,:)=[tmpdispls(p-1,2)+1,i]
           counts(p)=tmpdispls(p,2)+1-tmpdispls(p,1)
        end if
        natom=0
        
     end if
  end do
  tmpdispls(numprocs,:)=[tmpdispls(numprocs-1,2)+1,structures]
  counts(numprocs)=tmpdispls(numprocs,2)+1-tmpdispls(numprocs,1)
  displs=tmpdispls(myid+1,:)
  
  
end subroutine structure2process

!********************************************************

      subroutine getran(zran)
      implicit none

      integer time
      integer iseed
      real*8 zran
      real*8 ran0
      logical ldebug

      save iseed
! caution: iseed is changed in ran0 and this is necessary!

      if(iseed.eq.0) iseed=time()
      zran=ran0(iseed)

      return
      end

!*******************************************************
! Numerical Recipes

      FUNCTION ran0(idum)
      INTEGER idum,IA,IM,IQ,IR,MASK
      REAL*8 ran0,AM ! *8 is necessary for run on Opteron with pgf90
!      REAL ran0,AM
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,&
        MASK=123459876)
      INTEGER k
      idum=ieor(idum,MASK)
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      ran0=AM*idum
      idum=ieor(idum,MASK)
      return
      END

!********************************************************
      subroutine functions(counter,nat,nn,lattice,xyz,xx,dxdy,strs,&
                lforce,lperiodic,ldebug,e,volume)

! This subroutine is called once for each atomic geometry (configuration of the system)
! and calculates the set of functions xx for each atom in the system.
! The derivatives of the xx with respect to the x, y and z coordinate (dxdy) is also calculated.
! dxdy is the derivative dxx/dr, to be combined outside this routine with dE/dxx to dE/dr
      implicit none

      integer nat
      integer listdim
      parameter (listdim=50000)
      integer nn,nnn
      integer lsta(2,nat)
      integer lstc(listdim)
      real*8 lstb(listdim,4)
      integer i,j,k,h,m,n,ii,kk,jj
      integer counter
      integer ndummy,mdummy
      integer ntype

      real*8 cutoff
      real*8 cutoffg
      real*8 cshift
      real*8 rshift
      real*8 lattice(3,3)
      real*8 xyz(nat,3)
      real*8 xx(nn,1,nat) ! second index is for geometry
      real*8 dxx(nn,nat,nat)
      real*8 rij,rik,rih,rjk,rjh,rkh,fij
      real*8 alpha
      real*8 fcutij,fcutik,fcutih
      real*8 fcutjk,fcutjh,fcutkh
      real*8 pi
      real*8 costheta,theta
      real*8 norm1(3) ! norm vector of plane
      real*8 norm2(3)
      real*8 vec12(3),vec13(3) ! auxiliary vectors
      real*8 expxyz
      real*8 dexpxyzdxi
      real*8 dexpxyzdyi
      real*8 dexpxyzdzi
      real*8 dexpxyzdxj
      real*8 dexpxyzdyj
      real*8 dexpxyzdzj
      real*8 dexpxyzdxk
      real*8 dexpxyzdyk
      real*8 dexpxyzdzk
      real*8 temp1,temp2,temp3,temp4
      real*8 dangvec(3)
      real*8 dangvector

      real*8 dgdr(nn,nat,3)
      real*8 dxdy(nn,nat,nat,3)
      real*8 ddxxyy(nn,nat,nat,3) ! derivatives of symmetry functions with respect to x,y and z
! nn is number of symmetry functions
! 1. nat is the index of the atom the symmetry functions refer to
! 2. nat is the index for the atom with the coordinate we want the derivative for
! 3 is for x y z derivative

! temporary derivatives
      real*8 drijdx, drijdy, drijdz
      real*8 drijdxi, drijdyi, drijdzi
      real*8 drijdxj, drijdyj, drijdzj
      real*8 drijdxk, drijdyk, drijdzk
      real*8 drikdx, drikdy, drikdz
      real*8 drikdxi, drikdyi, drikdzi
      real*8 drikdxj, drikdyj, drikdzj
      real*8 drikdxk, drikdyk, drikdzk
      real*8 drihdx, drihdy, drihdz
      real*8 drihdxi, drihdyi, drihdzi
      real*8 drihdxh, drihdyh, drihdzh
      real*8 drihdxj, drihdyj, drihdzj
      real*8 drihdxk, drihdyk, drihdzk
      real*8 drjkdx, drjkdy, drjkdz
      real*8 drjkdxi, drjkdyi, drjkdzi
      real*8 drjkdxj, drjkdyj, drjkdzj
      real*8 drjkdxk, drjkdyk, drjkdzk
      real*8 dcutoffdx, dcutoffdy, dcutoffdz
      real*8 dfcutijdx,dfcutijdy,dfcutijdz
      real*8 dfcutijdxi,dfcutijdyi,dfcutijdzi
      real*8 dfcutijdxj,dfcutijdyj,dfcutijdzj
      real*8 dfcutijdxk,dfcutijdyk,dfcutijdzk
      real*8 dfcutikdx,dfcutikdy,dfcutikdz
      real*8 dfcutikdxi,dfcutikdyi,dfcutikdzi
      real*8 dfcutikdxj,dfcutikdyj,dfcutikdzj
      real*8 dfcutikdxk,dfcutikdyk,dfcutikdzk
      real*8 dfcutjkdx,dfcutjkdy,dfcutjkdz
      real*8 dfcutjkdxi,dfcutjkdyi,dfcutjkdzi
      real*8 dfcutjkdxj,dfcutjkdyj,dfcutjkdzj
      real*8 dfcutjkdxk,dfcutjkdyk,dfcutjkdzk
      real*8 dcosthetadx,dcosthetady,dcosthetadz
      real*8 dcosthetadxi,dcosthetadyi,dcosthetadzi
      real*8 dcosthetadxj,dcosthetadyj,dcosthetadzj
      real*8 dcosthetadxk,dcosthetadyk,dcosthetadzk
      real*8 f,g
      real*8 dfdx,dfdy,dfdz
      real*8 dfdxi,dfdyi,dfdzi
      real*8 dfdxj,dfdyj,dfdzj
      real*8 dfdxk,dfdyk,dfdzk
      real*8 dgdx,dgdy,dgdz
      real*8 dgdxi,dgdyi,dgdzi
      real*8 dgdxj,dgdyj,dgdzj
      real*8 dgdxk,dgdyk,dgdzk
      real*8 dthetadxi,dthetadyi,dthetadzi
      real*8 dthetadxj,dthetadyj,dthetadzj
      real*8 dthetadxk,dthetadyk,dthetadzk
      real*8 e,volume
! stress
      real*8 deltaxj,deltayj,deltazj
      real*8 deltaxk,deltayk,deltazk
      real*8 strs(3,3,nn,nat)

      logical ldebug
      logical ltest
      logical lforce
      logical lperiodic

!      ltest=.true.
      ltest=.false.
!      lforce=.true.

!      write(*,*)'functions.f starts'
      if(ltest) write(*,*)'functions.f starts'

! parameters
!      cutoffg=6.d0 ! global (largest) cutoffg
      pi=3.141592654d0

      open(66,file='neural.in',form='formatted',status='old')
      rewind(66)
      read(66,*)cutoffg
      close(66)

! translate atoms back into unit cell (needed for proper neighbor list)
      if(lperiodic)then
      call translate(nat,lattice,xyz,ldebug)
      endif

! calculate neighbor list
      call neighbor(nat,listdim,lsta,lstb,lstc,&
          cutoffg,lattice,xyz,lperiodic,ldebug,e,volume)

! initialization
!      xx(:,:,:)=0.0d0
!      dxdy(:,:,:,:)=0.0d0
      do j=1,nn
       do i=1,nat
        xx(j,1,i)=0.0d0   ! XXj value for atom i
        do m=1,nat
          dxdy(j,i,m,1)=0.0d0 ! x derivative of XXj for atom i
          dxdy(j,i,m,2)=0.0d0 ! y derivative of XXj for atom i
          dxdy(j,i,m,3)=0.0d0 ! z derivative of XXj for atom i
        enddo ! m
       enddo ! i
      enddo ! j

! read neural.in

      open(66,file='neural.in',form='formatted',status='old')
      rewind(66)
      read(66,*)cutoffg
      read(66,*)nnn
      if(ltest)write(*,*)'nnn from neural.in ',nnn
      if(nnn.ne.nn)then
        write(*,*)'nnn in neural.in does not fit to nn(0) in input'
        stop
      endif

!! initialization
      do ii=1,nnn
      do i=1,nat
       do j=1,3
        do k=1,3
        strs(j,k,ii,i)=0.0d0
        enddo
       enddo
      enddo
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! loop over all specified symmetry functions
      do ii=1,nnn
        read(66,*)ntype,alpha,cutoff
        if(ltest)write(*,*)'card is ',ntype,alpha,cutoff
        if(cutoff.gt.cutoffg) then
          write(*,*)'Error, cutoff is larger than cutoffg for neighbor'
          stop
        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TWO BODY TERM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(ntype.eq.1) then ! simple pair term

      do i=1,nat
        if(ltest) write(*,*)'calculating pair term for atom ',i

      do j=lsta(1,i),lsta(2,i) ! loop over all neighbors of atom i
        n=lstc(j)
        rij=lstb(j,4)
        if(rij.le.0.5d0) then
          write(*,*)'Error rij .le. 0.5 Angstrom',ntype
!          stop
        endif

        if(rij.le.cutoff) then ! rij < cutoff

        deltaxj=-1.d0*(xyz(i,1)-lstb(j,1))
        deltayj=-1.d0*(xyz(i,2)-lstb(j,2))
        deltazj=-1.d0*(xyz(i,3)-lstb(j,3))
        drijdxi=-deltaxj/rij
        drijdyi=-deltayj/rij
        drijdzi=-deltazj/rij
        drijdxj=-1.d0*drijdxi
        drijdyj=-1.d0*drijdyi
        drijdzj=-1.d0*drijdzi

        fcutij=0.5d0*(dcos(pi*rij/cutoff)+1.d0)
        temp1=0.5d0*(-dsin(pi*rij/cutoff))*(pi/cutoff)
        dfcutijdxi=temp1*drijdxi
        dfcutijdyi=temp1*drijdyi
        dfcutijdzi=temp1*drijdzi
        dfcutijdxj=-1.d0*dfcutijdxi
        dfcutijdyj=-1.d0*dfcutijdyi
        dfcutijdzj=-1.d0*dfcutijdzi

! summation of the Gaussians of all neighbors j within cutoff
        xx(ii,1,i)=xx(ii,1,i)+dexp(-1.d0*alpha*rij**2)*fcutij

!        if(lforce) then
! dxx/dx
        temp1=-2.d0*alpha*rij*dexp(-1.d0*alpha*rij**2)*fcutij
        temp2= dexp(-1.d0*alpha*rij**2)
        dxdy(ii,i,i,1)=dxdy(ii,i,i,1)+&
               (drijdxi*temp1&
              + temp2*dfcutijdxi)
        dxdy(ii,i,n,1)=dxdy(ii,i,n,1)+&
               (drijdxj*temp1&
              + temp2*dfcutijdxj)
        strs(1,1,ii,i)=strs(1,1,ii,i)+deltaxj*&
               (drijdxj*temp1&
              + temp2*dfcutijdxj)
        strs(2,1,ii,i)=strs(2,1,ii,i)+deltayj*&
               (drijdxj*temp1&
              + temp2*dfcutijdxj)
        strs(3,1,ii,i)=strs(3,1,ii,i)+deltazj*&
               (drijdxj*temp1&
              + temp2*dfcutijdxj)
! dxx/dy
        dxdy(ii,i,i,2)=dxdy(ii,i,i,2)+&
               (drijdyi*temp1&
              + temp2*dfcutijdyi)
        dxdy(ii,i,n,2)=dxdy(ii,i,n,2)+&
               (drijdyj*temp1&
              + temp2*dfcutijdyj)
        strs(1,2,ii,i)=strs(1,2,ii,i)+deltaxj*&
               (drijdyj*temp1&
              + temp2*dfcutijdyj)
        strs(2,2,ii,i)=strs(2,2,ii,i)+deltayj*&
               (drijdyj*temp1&
              + temp2*dfcutijdyj)
        strs(3,2,ii,i)=strs(3,2,ii,i)+deltazj*&
               (drijdyj*temp1&
              + temp2*dfcutijdyj)
! dxx/dz
        dxdy(ii,i,i,3)=dxdy(ii,i,i,3)+&
               (drijdzi*temp1&
              + temp2*dfcutijdzi)
        dxdy(ii,i,n,3)=dxdy(ii,i,n,3)+&
               (drijdzj*temp1&
              + temp2*dfcutijdzj)
        strs(1,3,ii,i)=strs(1,3,ii,i)+deltaxj*&
               (drijdzj*temp1&
              + temp2*dfcutijdzj)
        strs(2,3,ii,i)=strs(2,3,ii,i)+deltayj*&
               (drijdzj*temp1&
              + temp2*dfcutijdzj)
        strs(3,3,ii,i)=strs(3,3,ii,i)+deltazj*&
               (drijdzj*temp1&
              + temp2*dfcutijdzj)
!        endif ! lforce
        endif ! rij .le. cutoff
        enddo ! j
!       write(*,*)dxdy(ii,i,i,1),dxdy(ii,i,i,2),dxdy(ii,i,i,3)
      enddo ! i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THREE BODY TERM WITH GLOBAL CUTOFF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      elseif(ntype.eq.2) then ! three-body term with global cutoff
      do i=1,nat
        if(ltest) write(*,*)'calculating three body term for atom ',i

      do j=lsta(1,i),lsta(2,i) ! loop over all neighbors of this atom
        n=lstc(j)
        rij=lstb(j,4)

        if(rij.le.cutoff)then

        deltaxj=-1.d0*(xyz(i,1)-lstb(j,1))
        deltayj=-1.d0*(xyz(i,2)-lstb(j,2))
        deltazj=-1.d0*(xyz(i,3)-lstb(j,3))
        drijdxi=-deltaxj/rij
        drijdyi=-deltayj/rij
        drijdzi=-deltazj/rij
        drijdxj=-1.d0*drijdxi
        drijdyj=-1.d0*drijdyi
        drijdzj=-1.d0*drijdzi
        drijdxk=0.0d0
        drijdyk=0.0d0
        drijdzk=0.0d0

        fcutij=0.5d0*(dcos(pi*rij/cutoff)+1.d0)
        temp1=0.5d0*(-dsin(pi*rij/cutoff))*(pi/cutoff)
        dfcutijdxi=temp1*drijdxi
        dfcutijdyi=temp1*drijdyi
        dfcutijdzi=temp1*drijdzi
        dfcutijdxj=-1.d0*dfcutijdxi
        dfcutijdyj=-1.d0*dfcutijdyi
        dfcutijdzj=-1.d0*dfcutijdzi
        dfcutijdxk=0.0d0
        dfcutijdyk=0.0d0
        dfcutijdzk=0.0d0

        do k=lsta(1,i),lsta(2,i) ! loop over all neighbors of this atom
          m=lstc(k)
!          if(j.ne.k) then ! neighbors j and k must be different
          if(k.gt.j) then ! neighbors j and k must be different

          rik=lstb(k,4)

          if(rik.le.cutoff)then

          deltaxk=-1.d0*(xyz(i,1)-lstb(k,1))
          deltayk=-1.d0*(xyz(i,2)-lstb(k,2))
          deltazk=-1.d0*(xyz(i,3)-lstb(k,3))
          drikdxi=-deltaxk/rik
          drikdyi=-deltayk/rik
          drikdzi=-deltazk/rik
          drikdxk=-1.d0*drikdxi
          drikdyk=-1.d0*drikdyi
          drikdzk=-1.d0*drikdzi
          drikdxj=0.0d0
          drikdyj=0.0d0
          drikdzj=0.0d0

          fcutik=0.5d0*(dcos(pi*rik/cutoff)+1.d0)
          temp1=0.5d0*(-dsin(pi*rik/cutoff))*(pi/cutoff)
          dfcutikdxi=temp1*drikdxi
          dfcutikdyi=temp1*drikdyi
          dfcutikdzi=temp1*drikdzi
          dfcutikdxj=0.0d0
          dfcutikdyj=0.0d0
          dfcutikdzj=0.0d0
          dfcutikdxk=-1.d0*dfcutikdxi
          dfcutikdyk=-1.d0*dfcutikdyi
          dfcutikdzk=-1.d0*dfcutikdzi

          rjk=(lstb(j,1)-lstb(k,1))**2 +&
              (lstb(j,2)-lstb(k,2))**2 +&
              (lstb(j,3)-lstb(k,3))**2
          rjk=dsqrt(rjk)

          if(rjk.le.0.5d0) then
            write(*,*)'Error rjk .le. 0.5',ntype
!            stop
          endif

          if(rjk.le.cutoff) then
          drjkdxj=(lstb(j,1)-lstb(k,1))/rjk
          drjkdyj=(lstb(j,2)-lstb(k,2))/rjk
          drjkdzj=(lstb(j,3)-lstb(k,3))/rjk
          drjkdxk=-1.d0*drjkdxj
          drjkdyk=-1.d0*drjkdyj
          drjkdzk=-1.d0*drjkdzj
          drjkdxi=0.0d0
          drjkdyi=0.0d0
          drjkdzi=0.0d0

          fcutjk=0.5d0*(dcos(pi*rjk/cutoff)+1.d0)
          temp1=0.5d0*(-dsin(pi*rjk/cutoff))*(pi/cutoff)
          dfcutjkdxj=temp1*drjkdxj
          dfcutjkdyj=temp1*drjkdyj
          dfcutjkdzj=temp1*drjkdzj
          dfcutjkdxk=-1.d0*dfcutjkdxj
          dfcutjkdyk=-1.d0*dfcutjkdyj
          dfcutjkdzk=-1.d0*dfcutjkdzj
          dfcutjkdxi=0.0d0
          dfcutjkdyi=0.0d0
          dfcutjkdzi=0.0d0

! costheta=(rjk**2 - rij**2 -rik**2)/(-2.d0*rij*rik)
            f=rjk**2 - rij**2 -rik**2
            g=-2.d0*rij*rik
            costheta=f/g
            costheta=costheta+1.d0 ! avoid negative values

! calculate the derivatives of costheta
! (f/g)' = (f'g - fg')/g^2

            if(lforce) then
               dfdxi=-2.d0*rij*drijdxi - 2.d0*rik*drikdxi
               dfdyi=-2.d0*rij*drijdyi - 2.d0*rik*drikdyi
               dfdzi=-2.d0*rij*drijdzi - 2.d0*rik*drikdzi

               dfdxj=2.d0*rjk*drjkdxj - 2.d0*rij*drijdxj
               dfdyj=2.d0*rjk*drjkdyj - 2.d0*rij*drijdyj
               dfdzj=2.d0*rjk*drjkdzj - 2.d0*rij*drijdzj

               dfdxk=2.d0*rjk*drjkdxk - 2.d0*rik*drikdxk
               dfdyk=2.d0*rjk*drjkdyk - 2.d0*rik*drikdyk
               dfdzk=2.d0*rjk*drjkdzk - 2.d0*rik*drikdzk

               dgdxi=-2.d0*(drijdxi*rik + rij*drikdxi)
               dgdyi=-2.d0*(drijdyi*rik + rij*drikdyi)
               dgdzi=-2.d0*(drijdzi*rik + rij*drikdzi)

               dgdxj=-2.d0*drijdxj*rik
               dgdyj=-2.d0*drijdyj*rik
               dgdzj=-2.d0*drijdzj*rik

               dgdxk=-2.d0*rij*drikdxk
               dgdyk=-2.d0*rij*drikdyk
               dgdzk=-2.d0*rij*drikdzk

               dcosthetadxi=(dfdxi*g - f*dgdxi)/g**2
               dcosthetadyi=(dfdyi*g - f*dgdyi)/g**2
               dcosthetadzi=(dfdzi*g - f*dgdzi)/g**2
               dcosthetadxj=(dfdxj*g - f*dgdxj)/g**2
               dcosthetadyj=(dfdyj*g - f*dgdyj)/g**2
               dcosthetadzj=(dfdzj*g - f*dgdzj)/g**2
               dcosthetadxk=(dfdxk*g - f*dgdxk)/g**2
               dcosthetadyk=(dfdyk*g - f*dgdyk)/g**2
               dcosthetadzk=(dfdzk*g - f*dgdzk)/g**2
            endif

! calculation of exponential term
      expxyz=dexp(-alpha*(rij**2+rik**2+rjk**2))
      temp1=-alpha*2.0d0*expxyz
      dexpxyzdxi=(rij*drijdxi+rik*drikdxi+rjk*drjkdxi)*temp1
      dexpxyzdyi=(rij*drijdyi+rik*drikdyi+rjk*drjkdyi)*temp1
      dexpxyzdzi=(rij*drijdzi+rik*drikdzi+rjk*drjkdzi)*temp1
      dexpxyzdxj=(rij*drijdxj+rik*drikdxj+rjk*drjkdxj)*temp1
      dexpxyzdyj=(rij*drijdyj+rik*drikdyj+rjk*drjkdyj)*temp1
      dexpxyzdzj=(rij*drijdzj+rik*drikdzj+rjk*drjkdzj)*temp1
      dexpxyzdxk=(rij*drijdxk+rik*drikdxk+rjk*drjkdxk)*temp1
      dexpxyzdyk=(rij*drijdyk+rik*drikdyk+rjk*drjkdyk)*temp1
      dexpxyzdzk=(rij*drijdzk+rik*drikdzk+rjk*drjkdzk)*temp1

      xx(ii,1,i)=xx(ii,1,i)+costheta*expxyz*fcutij*fcutik*fcutjk

! (fgh)' = f'gh + fg'h + fgh'
! (fghi)'= f'ghi + fg'hi + fgh'i + fghi'

! for the derivatives fcutjk is a constant

      if(lforce) then
! dxxii/dx_i
        temp1=(+dcosthetadxi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxi)
        temp2=(+dcosthetadxj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxj)
        temp3=(+dcosthetadxk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxk)
        dxdy(ii,i,i,1)=dxdy(ii,i,i,1)+temp1
        dxdy(ii,i,n,1)=dxdy(ii,i,n,1)+temp2
        dxdy(ii,i,m,1)=dxdy(ii,i,m,1)+temp3
        strs(1,1,ii,i)=strs(1,1,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,1,ii,i)=strs(2,1,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,1,ii,i)=strs(3,1,ii,i)+deltazj*temp2&
              +deltazk*temp3
! dxxii/dy_i
        temp1=(+dcosthetadyi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyi)
        temp2=(+dcosthetadyj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyj)
        temp3=(+dcosthetadyk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyk)
        dxdy(ii,i,i,2)=dxdy(ii,i,i,2)+temp1
        dxdy(ii,i,n,2)=dxdy(ii,i,n,2)+temp2
        dxdy(ii,i,m,2)=dxdy(ii,i,m,2)+temp3
        strs(1,2,ii,i)=strs(1,2,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,2,ii,i)=strs(2,2,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,2,ii,i)=strs(3,2,ii,i)+deltazj*temp2&
              +deltazk*temp3
!! dxxii/dz_i
        temp1=(+dcosthetadzi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzi)
        temp2=(+dcosthetadzj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzj)
        temp3=(+dcosthetadzk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzk)
        dxdy(ii,i,i,3)=dxdy(ii,i,i,3)+temp1
        dxdy(ii,i,n,3)=dxdy(ii,i,n,3)+temp2
        dxdy(ii,i,m,3)=dxdy(ii,i,m,3)+temp3
        strs(1,3,ii,i)=strs(1,3,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,3,ii,i)=strs(2,3,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,3,ii,i)=strs(3,3,ii,i)+deltazj*temp2&
              +deltazk*temp3
      endif

            endif ! rjk .le. cutoff
            endif ! rik .le. cutoff
          endif ! j .ne. k
        enddo ! k
        endif ! rij .le. cutoff
      enddo ! j
      enddo ! i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TWO BODY TERM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      elseif(ntype.eq.3) then ! Gaussian sphere (pair term)
! caution: alpha has the meaning of the shift in the input file
      rshift=alpha
      alpha=4.d0

      do i=1,nat

      do j=lsta(1,i),lsta(2,i) ! loop over all neighbors of atom i
        n=lstc(j)
        rij=lstb(j,4)

        if(rij.le.cutoff) then ! rij < cutoff
        deltaxj=-1.d0*(xyz(i,1)-lstb(j,1))
        deltayj=-1.d0*(xyz(i,2)-lstb(j,2))
        deltazj=-1.d0*(xyz(i,3)-lstb(j,3))
        drijdxi=-deltaxj/rij
        drijdyi=-deltayj/rij
        drijdzi=-deltazj/rij
        drijdxj=-1.d0*drijdxi
        drijdyj=-1.d0*drijdyi
        drijdzj=-1.d0*drijdzi

        fcutij=0.5d0*(dcos(pi*rij/cutoff)+1.d0)
        temp1=0.5d0*(-dsin(pi*rij/cutoff))*(pi/cutoff)
        dfcutijdxi=temp1*drijdxi
        dfcutijdyi=temp1*drijdyi
        dfcutijdzi=temp1*drijdzi
        dfcutijdxj=-1.d0*dfcutijdxi
        dfcutijdyj=-1.d0*dfcutijdyi
        dfcutijdzj=-1.d0*dfcutijdzi

! summation of the Gaussians of all neighbors j within cutoff
        xx(ii,1,i)=xx(ii,1,i)+dexp(-1.d0*alpha*(rij-rshift)**2)*fcutij

        if(lforce)then
        temp1=-2.d0*alpha*(rij-rshift)
        temp2=dexp(-1.d0*alpha*(rij-rshift)**2)
! dxx/dx
        dxdy(ii,i,i,1)=dxdy(ii,i,i,1)+&
               (temp1*drijdxi*&
                temp2*fcutij&
              + temp2*dfcutijdxi)
        dxdy(ii,i,n,1)=dxdy(ii,i,n,1)+&
               (temp1*drijdxj*&
                temp2*fcutij&
              + temp2*dfcutijdxj)
        strs(1,1,ii,i)=strs(1,1,ii,i)+deltaxj*&
               (temp1*drijdxj*&
                temp2*fcutij&
              + temp2*dfcutijdxj)
        strs(2,1,ii,i)=strs(2,1,ii,i)+deltayj*&
               (temp1*drijdxj*&
                temp2*fcutij&
              + temp2*dfcutijdxj)
        strs(3,1,ii,i)=strs(3,1,ii,i)+deltazj*&
               (temp1*drijdxj*&
                temp2*fcutij&
              + temp2*dfcutijdxj)
! dxx/dy
        dxdy(ii,i,i,2)=dxdy(ii,i,i,2)+&
               (temp1*&
                drijdyi*temp2*fcutij&
              + temp2*dfcutijdyi)
        dxdy(ii,i,n,2)=dxdy(ii,i,n,2)+&
               (temp1*&
                drijdyj*temp2*fcutij&
              + temp2*dfcutijdyj)
        strs(1,2,ii,i)=strs(1,2,ii,i)+deltaxj*&
               (temp1*&
                drijdyj*temp2*fcutij&
              + temp2*dfcutijdyj)
        strs(2,2,ii,i)=strs(2,2,ii,i)+deltayj*&
               (temp1*&
                drijdyj*temp2*fcutij&
              + temp2*dfcutijdyj)
        strs(3,2,ii,i)=strs(3,2,ii,i)+deltazj*&
               (temp1*&
                drijdyj*temp2*fcutij&
              + temp2*dfcutijdyj)
! dxx/dz
        dxdy(ii,i,i,3)=dxdy(ii,i,i,3)+&
               (temp1*&
                drijdzi*temp2*fcutij&
              + temp2*dfcutijdzi)
        dxdy(ii,i,n,3)=dxdy(ii,i,n,3)+&
               (temp1*&
                drijdzj*temp2*fcutij&
              + temp2*dfcutijdzj)
        strs(1,3,ii,i)=strs(1,3,ii,i)+deltaxj*&
               (temp1*&
                drijdzj*temp2*fcutij&
              + temp2*dfcutijdzj)
        strs(2,3,ii,i)=strs(2,3,ii,i)+deltayj*&
               (temp1*&
                drijdzj*temp2*fcutij&
              + temp2*dfcutijdzj)
        strs(3,3,ii,i)=strs(3,3,ii,i)+deltazj*&
               (temp1*&
                drijdzj*temp2*fcutij&
              + temp2*dfcutijdzj)
        endif

        endif ! rij .le. cutoff
        enddo ! j
      enddo ! i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THREE BODY TERM WITH GLOBAL CUTOFF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      elseif(ntype.eq.4) then ! three-body term with global cutoff
      do i=1,nat
        if(ltest) write(*,*)'calculating three body term for atom ',i

      do j=lsta(1,i),lsta(2,i) ! loop over all neighbors of this atom
        n=lstc(j)
        rij=lstb(j,4)

        if(rij.le.cutoff)then

        deltaxj=-1.d0*(xyz(i,1)-lstb(j,1))
        deltayj=-1.d0*(xyz(i,2)-lstb(j,2))
        deltazj=-1.d0*(xyz(i,3)-lstb(j,3))
        drijdxi=-deltaxj/rij
        drijdyi=-deltayj/rij
        drijdzi=-deltazj/rij
        drijdxj=-1.d0*drijdxi
        drijdyj=-1.d0*drijdyi
        drijdzj=-1.d0*drijdzi
        drijdxk=0.0d0
        drijdyk=0.0d0
        drijdzk=0.0d0

        fcutij=0.5d0*(dcos(pi*rij/cutoff)+1.d0)
        temp1=0.5d0*(-dsin(pi*rij/cutoff))*(pi/cutoff)
        dfcutijdxi=temp1*drijdxi
        dfcutijdyi=temp1*drijdyi
        dfcutijdzi=temp1*drijdzi
        dfcutijdxj=-1.d0*dfcutijdxi
        dfcutijdyj=-1.d0*dfcutijdyi
        dfcutijdzj=-1.d0*dfcutijdzi
        dfcutijdxk=0.0d0
        dfcutijdyk=0.0d0
        dfcutijdzk=0.0d0

        do k=lsta(1,i),lsta(2,i) ! loop over all neighbors of this atom
          m=lstc(k)

!          if(j.ne.k) then ! neighbors j and k must be different
          if(k.gt.j) then ! neighbors j and k must be different

          rik=lstb(k,4)

          if(rik.le.cutoff)then

          deltaxk=-1.d0*(xyz(i,1)-lstb(k,1))
          deltayk=-1.d0*(xyz(i,2)-lstb(k,2))
          deltazk=-1.d0*(xyz(i,3)-lstb(k,3))
          drikdxi=(xyz(i,1)-lstb(k,1))/rik
          drikdyi=(xyz(i,2)-lstb(k,2))/rik
          drikdzi=(xyz(i,3)-lstb(k,3))/rik
          drikdxk=-1.d0*drikdxi
          drikdyk=-1.d0*drikdyi
          drikdzk=-1.d0*drikdzi
          drikdxj=0.0d0
          drikdyj=0.0d0
          drikdzj=0.0d0

          fcutik=0.5d0*(dcos(pi*rik/cutoff)+1.d0)
          temp1=0.5d0*(-dsin(pi*rik/cutoff))*(pi/cutoff)
          dfcutikdxi=temp1*drikdxi
          dfcutikdyi=temp1*drikdyi
          dfcutikdzi=temp1*drikdzi
          dfcutikdxj=0.0d0
          dfcutikdyj=0.0d0
          dfcutikdzj=0.0d0
          dfcutikdxk=-1.d0*dfcutikdxi
          dfcutikdyk=-1.d0*dfcutikdyi
          dfcutikdzk=-1.d0*dfcutikdzi

          rjk=(lstb(j,1)-lstb(k,1))**2 +&
              (lstb(j,2)-lstb(k,2))**2 +&
              (lstb(j,3)-lstb(k,3))**2
          rjk=dsqrt(rjk)

          if(rjk.le.cutoff)then

          drjkdxj=(lstb(j,1)-lstb(k,1))/rjk
          drjkdyj=(lstb(j,2)-lstb(k,2))/rjk
          drjkdzj=(lstb(j,3)-lstb(k,3))/rjk
          drjkdxk=-1.d0*drjkdxj
          drjkdyk=-1.d0*drjkdyj
          drjkdzk=-1.d0*drjkdzj
          drjkdxi=0.0d0
          drjkdyi=0.0d0
          drjkdzi=0.0d0

          fcutjk=0.5d0*(dcos(pi*rjk/cutoff)+1.d0)
          temp1=0.5d0*(-dsin(pi*rjk/cutoff))*(pi/cutoff)
          dfcutjkdxj=temp1*drjkdxj
          dfcutjkdyj=temp1*drjkdyj
          dfcutjkdzj=temp1*drjkdzj
          dfcutjkdxk=-1.d0*dfcutjkdxj
          dfcutjkdyk=-1.d0*dfcutjkdyj
          dfcutjkdzk=-1.d0*dfcutjkdzj
          dfcutjkdxi=0.0d0
          dfcutjkdyi=0.0d0
          dfcutjkdzi=0.0d0

! costheta=(rjk**2 - rij**2 -rik**2)/(-2.d0*rij*rik)
          f=rjk**2 - rij**2 -rik**2
          g=-2.d0*rij*rik
          costheta=f/g
          costheta=1.d0-costheta ! avoid negative values

! calculate the derivatives of costheta
! (f/g)' = (f'g - fg')/g^2

            if(lforce)then
               dfdxi=-2.d0*rij*drijdxi - 2.d0*rik*drikdxi
               dfdyi=-2.d0*rij*drijdyi - 2.d0*rik*drikdyi
               dfdzi=-2.d0*rij*drijdzi - 2.d0*rik*drikdzi

               dfdxj=2.d0*rjk*drjkdxj - 2.d0*rij*drijdxj
               dfdyj=2.d0*rjk*drjkdyj - 2.d0*rij*drijdyj
               dfdzj=2.d0*rjk*drjkdzj - 2.d0*rij*drijdzj

               dfdxk=2.d0*rjk*drjkdxk - 2.d0*rik*drikdxk
               dfdyk=2.d0*rjk*drjkdyk - 2.d0*rik*drikdyk
               dfdzk=2.d0*rjk*drjkdzk - 2.d0*rik*drikdzk

               dgdxi=-2.d0*(drijdxi*rik + rij*drikdxi)
               dgdyi=-2.d0*(drijdyi*rik + rij*drikdyi)
               dgdzi=-2.d0*(drijdzi*rik + rij*drikdzi)

               dgdxj=-2.d0*drijdxj*rik
               dgdyj=-2.d0*drijdyj*rik
               dgdzj=-2.d0*drijdzj*rik

               dgdxk=-2.d0*rij*drikdxk
               dgdyk=-2.d0*rij*drikdyk
               dgdzk=-2.d0*rij*drikdzk

               dcosthetadxi=-(dfdxi*g - f*dgdxi)/g**2
               dcosthetadyi=-(dfdyi*g - f*dgdyi)/g**2
               dcosthetadzi=-(dfdzi*g - f*dgdzi)/g**2
               dcosthetadxj=-(dfdxj*g - f*dgdxj)/g**2
               dcosthetadyj=-(dfdyj*g - f*dgdyj)/g**2
               dcosthetadzj=-(dfdzj*g - f*dgdzj)/g**2
               dcosthetadxk=-(dfdxk*g - f*dgdxk)/g**2
               dcosthetadyk=-(dfdyk*g - f*dgdyk)/g**2
               dcosthetadzk=-(dfdzk*g - f*dgdzk)/g**2
            endif

! calculation of exponential term
      expxyz=dexp(-alpha*(rij**2+rik**2+rjk**2))
      temp1=-alpha*2.0d0*expxyz
      dexpxyzdxi=(rij*drijdxi+rik*drikdxi+rjk*drjkdxi)*temp1
      dexpxyzdyi=(rij*drijdyi+rik*drikdyi+rjk*drjkdyi)*temp1
      dexpxyzdzi=(rij*drijdzi+rik*drikdzi+rjk*drjkdzi)*temp1
      dexpxyzdxj=(rij*drijdxj+rik*drikdxj+rjk*drjkdxj)*temp1
      dexpxyzdyj=(rij*drijdyj+rik*drikdyj+rjk*drjkdyj)*temp1
      dexpxyzdzj=(rij*drijdzj+rik*drikdzj+rjk*drjkdzj)*temp1
      dexpxyzdxk=(rij*drijdxk+rik*drikdxk+rjk*drjkdxk)*temp1
      dexpxyzdyk=(rij*drijdyk+rik*drikdyk+rjk*drjkdyk)*temp1
      dexpxyzdzk=(rij*drijdzk+rik*drikdzk+rjk*drjkdzk)*temp1

      xx(ii,1,i)=xx(ii,1,i)+costheta*expxyz*fcutij*fcutik*fcutjk

! (fgh)' = f'gh + fg'h + fgh'
! (fghi)'= f'ghi + fg'hi + fgh'i + fghi'

! for the derivatives fcutjk is a constant

      if(lforce)then
! dxxii/dx_i
        temp1=(+dcosthetadxi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxi)
        temp2=(+dcosthetadxj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxj)
        temp3=(+dcosthetadxk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxk)
        dxdy(ii,i,i,1)=dxdy(ii,i,i,1)+temp1
        dxdy(ii,i,n,1)=dxdy(ii,i,n,1)+temp2
        dxdy(ii,i,m,1)=dxdy(ii,i,m,1)+temp3
        strs(1,1,ii,i)=strs(1,1,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,1,ii,i)=strs(2,1,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,1,ii,i)=strs(3,1,ii,i)+deltazj*temp2&
              +deltazk*temp3
! dxxii/dy_i
        temp1=(+dcosthetadyi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyi)
        temp2=(+dcosthetadyj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyj)
        temp3=(+dcosthetadyk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyk)
        dxdy(ii,i,i,2)=dxdy(ii,i,i,2)+temp1
        dxdy(ii,i,n,2)=dxdy(ii,i,n,2)+temp2
        dxdy(ii,i,m,2)=dxdy(ii,i,m,2)+temp3
        strs(1,2,ii,i)=strs(1,2,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,2,ii,i)=strs(2,2,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,2,ii,i)=strs(3,2,ii,i)+deltazj*temp2&
              +deltazk*temp3
! dxxii/dz_i
        temp1=(+dcosthetadzi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzi)
        temp2=(+dcosthetadzj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzj)
        temp3=(+dcosthetadzk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzk)
        dxdy(ii,i,i,3)=dxdy(ii,i,i,3)+temp1
        dxdy(ii,i,n,3)=dxdy(ii,i,n,3)+temp2
        dxdy(ii,i,m,3)=dxdy(ii,i,m,3)+temp3
        strs(1,3,ii,i)=strs(1,3,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,3,ii,i)=strs(2,3,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,3,ii,i)=strs(3,3,ii,i)+deltazj*temp2&
              +deltazk*temp3
      endif

            endif ! rjk .le. cutoff
          endif ! rik .le. cutoff
          endif ! j .ne. k
        enddo ! k
        endif ! rij .le. cutoff
      enddo ! j
      enddo ! i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THREE BODY TERM WITH GLOBAL CUTOFF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      elseif(ntype.eq.5) then ! three-body term with global cutoff
      do i=1,nat
        if(ltest) write(*,*)'calculating three body term for atom ',i

      do j=lsta(1,i),lsta(2,i) ! loop over all neighbors of this atom
        n=lstc(j)
        rij=lstb(j,4)
        if(rij.le.0.5d0) then
          write(*,*)'Error rij .le. 0.5',ntype
!          stop
        endif

        if(rij.le.cutoff)then

        deltaxj=-1.d0*(xyz(i,1)-lstb(j,1))
        deltayj=-1.d0*(xyz(i,2)-lstb(j,2))
        deltazj=-1.d0*(xyz(i,3)-lstb(j,3))
        drijdxi=-deltaxj/rij
        drijdyi=-deltayj/rij
        drijdzi=-deltazj/rij
        drijdxj=-1.d0*drijdxi
        drijdyj=-1.d0*drijdyi
        drijdzj=-1.d0*drijdzi
        drijdxk=0.0d0
        drijdyk=0.0d0
        drijdzk=0.0d0

        fcutij=0.5d0*(dcos(pi*rij/cutoff)+1.d0)
        temp1=0.5d0*(-dsin(pi*rij/cutoff))*(pi/cutoff)
        dfcutijdxi=temp1*drijdxi
        dfcutijdyi=temp1*drijdyi
        dfcutijdzi=temp1*drijdzi
        dfcutijdxj=-1.d0*dfcutijdxi
        dfcutijdyj=-1.d0*dfcutijdyi
        dfcutijdzj=-1.d0*dfcutijdzi
        dfcutijdxk=0.0d0
        dfcutijdyk=0.0d0
        dfcutijdzk=0.0d0

        do k=lsta(1,i),lsta(2,i) ! loop over all neighbors of this atom
          m=lstc(k)

!          if(j.ne.k) then ! neighbors j and k must be different
          if(k.gt.j) then ! neighbors j and k must be different

          rik=lstb(k,4)

          if(rik.le.cutoff)then
          deltaxk=-1.d0*(xyz(i,1)-lstb(k,1))
          deltayk=-1.d0*(xyz(i,2)-lstb(k,2))
          deltazk=-1.d0*(xyz(i,3)-lstb(k,3))
          drikdxi=(xyz(i,1)-lstb(k,1))/rik
          drikdyi=(xyz(i,2)-lstb(k,2))/rik
          drikdzi=(xyz(i,3)-lstb(k,3))/rik
          drikdxk=-1.d0*drikdxi
          drikdyk=-1.d0*drikdyi
          drikdzk=-1.d0*drikdzi
          drikdxj=0.0d0
          drikdyj=0.0d0
          drikdzj=0.0d0

          fcutik=0.5d0*(dcos(pi*rik/cutoff)+1.d0)
          temp1=0.5d0*(-dsin(pi*rik/cutoff))*(pi/cutoff)
          dfcutikdxi=temp1*drikdxi
          dfcutikdyi=temp1*drikdyi
          dfcutikdzi=temp1*drikdzi
          dfcutikdxj=0.0d0
          dfcutikdyj=0.0d0
          dfcutikdzj=0.0d0
          dfcutikdxk=-1.d0*dfcutikdxi
          dfcutikdyk=-1.d0*dfcutikdyi
          dfcutikdzk=-1.d0*dfcutikdzi

          rjk=(lstb(j,1)-lstb(k,1))**2 +&
              (lstb(j,2)-lstb(k,2))**2 +&
              (lstb(j,3)-lstb(k,3))**2
          rjk=dsqrt(rjk)
          if(rjk.le.0.5d0) then
            write(*,*)'Error rjk .le. 0.5',ntype
!            stop
          endif

          if(rjk.le.cutoff)then

          drjkdxj=(lstb(j,1)-lstb(k,1))/rjk
          drjkdyj=(lstb(j,2)-lstb(k,2))/rjk
          drjkdzj=(lstb(j,3)-lstb(k,3))/rjk
          drjkdxk=-1.d0*drjkdxj
          drjkdyk=-1.d0*drjkdyj
          drjkdzk=-1.d0*drjkdzj
          drjkdxi=0.0d0
          drjkdyi=0.0d0
          drjkdzi=0.0d0

          fcutjk=0.5d0*(dcos(pi*rjk/cutoff)+1.d0)
          temp1=0.5d0*(-dsin(pi*rjk/cutoff))*(pi/cutoff)
          dfcutjkdxj=temp1*drjkdxj
          dfcutjkdyj=temp1*drjkdyj
          dfcutjkdzj=temp1*drjkdzj
          dfcutjkdxk=-1.d0*dfcutjkdxj
          dfcutjkdyk=-1.d0*dfcutjkdyj
          dfcutjkdzk=-1.d0*dfcutjkdzj
          dfcutjkdxi=0.0d0
          dfcutjkdyi=0.0d0
          dfcutjkdzi=0.0d0

! costheta=(rjk**2 - rij**2 -rik**2)/(-2.d0*rij*rik)
          f=rjk**2 - rij**2 -rik**2
          g=-2.d0*rij*rik
          costheta=f/g
          costheta=costheta+1.d0 ! avoid negative values
          costheta=0.5d0*(costheta**2)
! summary: costheta=0.5*(1+f/g)**2

! calculate the derivatives of costheta
! (f/g)' = (f'g - fg')/g^2

            if(lforce)then
               dfdxi=-2.d0*rij*drijdxi - 2.d0*rik*drikdxi
               dfdyi=-2.d0*rij*drijdyi - 2.d0*rik*drikdyi
               dfdzi=-2.d0*rij*drijdzi - 2.d0*rik*drikdzi

               dfdxj=2.d0*rjk*drjkdxj - 2.d0*rij*drijdxj
               dfdyj=2.d0*rjk*drjkdyj - 2.d0*rij*drijdyj
               dfdzj=2.d0*rjk*drjkdzj - 2.d0*rij*drijdzj

               dfdxk=2.d0*rjk*drjkdxk - 2.d0*rik*drikdxk
               dfdyk=2.d0*rjk*drjkdyk - 2.d0*rik*drikdyk
               dfdzk=2.d0*rjk*drjkdzk - 2.d0*rik*drikdzk

               dgdxi=-2.d0*(drijdxi*rik + rij*drikdxi)
               dgdyi=-2.d0*(drijdyi*rik + rij*drikdyi)
               dgdzi=-2.d0*(drijdzi*rik + rij*drikdzi)

               dgdxj=-2.d0*drijdxj*rik
               dgdyj=-2.d0*drijdyj*rik
               dgdzj=-2.d0*drijdzj*rik

               dgdxk=-2.d0*rij*drikdxk
               dgdyk=-2.d0*rij*drikdyk
               dgdzk=-2.d0*rij*drikdzk

               dcosthetadxi=(1.d0+f/g)*(dfdxi*g - f*dgdxi)/g**2
               dcosthetadyi=(1.d0+f/g)*(dfdyi*g - f*dgdyi)/g**2
               dcosthetadzi=(1.d0+f/g)*(dfdzi*g - f*dgdzi)/g**2
               dcosthetadxj=(1.d0+f/g)*(dfdxj*g - f*dgdxj)/g**2
               dcosthetadyj=(1.d0+f/g)*(dfdyj*g - f*dgdyj)/g**2
               dcosthetadzj=(1.d0+f/g)*(dfdzj*g - f*dgdzj)/g**2
               dcosthetadxk=(1.d0+f/g)*(dfdxk*g - f*dgdxk)/g**2
               dcosthetadyk=(1.d0+f/g)*(dfdyk*g - f*dgdyk)/g**2
               dcosthetadzk=(1.d0+f/g)*(dfdzk*g - f*dgdzk)/g**2
            endif

! calculation of exponential term
      expxyz=dexp(-alpha*(rij**2+rik**2+rjk**2))
      temp1=-alpha*2.0d0*expxyz
      dexpxyzdxi=(rij*drijdxi+rik*drikdxi+rjk*drjkdxi)*temp1
      dexpxyzdyi=(rij*drijdyi+rik*drikdyi+rjk*drjkdyi)*temp1
      dexpxyzdzi=(rij*drijdzi+rik*drikdzi+rjk*drjkdzi)*temp1
      dexpxyzdxj=(rij*drijdxj+rik*drikdxj+rjk*drjkdxj)*temp1
      dexpxyzdyj=(rij*drijdyj+rik*drikdyj+rjk*drjkdyj)*temp1
      dexpxyzdzj=(rij*drijdzj+rik*drikdzj+rjk*drjkdzj)*temp1
      dexpxyzdxk=(rij*drijdxk+rik*drikdxk+rjk*drjkdxk)*temp1
      dexpxyzdyk=(rij*drijdyk+rik*drikdyk+rjk*drjkdyk)*temp1
      dexpxyzdzk=(rij*drijdzk+rik*drikdzk+rjk*drjkdzk)*temp1

      xx(ii,1,i)=xx(ii,1,i)+costheta*expxyz*fcutij*fcutik*fcutjk

! (fgh)' = f'gh + fg'h + fgh'
! (fghi)'= f'ghi + fg'hi + fgh'i + fghi'

      if(lforce)then
! for the derivatives fcutjk is a constant
        temp1=(+dcosthetadxi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxi)
        temp2=(+dcosthetadxj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxj)
        temp3=(+dcosthetadxk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxk)
        dxdy(ii,i,i,1)=dxdy(ii,i,i,1)+temp1
        dxdy(ii,i,n,1)=dxdy(ii,i,n,1)+temp2
        dxdy(ii,i,m,1)=dxdy(ii,i,m,1)+temp3
        strs(1,1,ii,i)=strs(1,1,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,1,ii,i)=strs(2,1,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,1,ii,i)=strs(3,1,ii,i)+deltazj*temp2&
              +deltazk*temp3

        temp1=(+dcosthetadyi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyi)
        temp2=(+dcosthetadyj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyj)
        temp3=(+dcosthetadyk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyk)
        dxdy(ii,i,i,2)=dxdy(ii,i,i,2)+temp1
        dxdy(ii,i,n,2)=dxdy(ii,i,n,2)+temp2
        dxdy(ii,i,m,2)=dxdy(ii,i,m,2)+temp3
        strs(1,2,ii,i)=strs(1,2,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,2,ii,i)=strs(2,2,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,2,ii,i)=strs(3,2,ii,i)+deltazj*temp2&
              +deltazk*temp3

        temp1=(+dcosthetadzi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzi)
        temp2=(+dcosthetadzj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzj)
        temp3=(+dcosthetadzk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzk)
        dxdy(ii,i,i,3)=dxdy(ii,i,i,3)+temp1
        dxdy(ii,i,n,3)=dxdy(ii,i,n,3)+temp2
        dxdy(ii,i,m,3)=dxdy(ii,i,m,3)+temp3
        strs(1,3,ii,i)=strs(1,3,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,3,ii,i)=strs(2,3,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,3,ii,i)=strs(3,3,ii,i)+deltazj*temp2&
              +deltazk*temp3
      endif !lforce

            endif ! rjk .le. cutoff
          endif ! rik .le. cutoff
          endif ! j .ne. k
        enddo ! k
        endif ! rij .le. cutoff
      enddo ! j
      enddo ! i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THREE BODY TERM WITH GLOBAL CUTOFF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      elseif(ntype.eq.6) then ! three-body term with global cutoff
      do i=1,nat
        if(ltest) write(*,*)'calculating three body term for atom ',i

      do j=lsta(1,i),lsta(2,i) ! loop over all neighbors of this atom
        n=lstc(j)
        rij=lstb(j,4)

        if(rij.le.cutoff)then

        deltaxj=-1.d0*(xyz(i,1)-lstb(j,1))
        deltayj=-1.d0*(xyz(i,2)-lstb(j,2))
        deltazj=-1.d0*(xyz(i,3)-lstb(j,3))
        drijdxi=-deltaxj/rij
        drijdyi=-deltayj/rij
        drijdzi=-deltazj/rij
        drijdxj=-1.d0*drijdxi
        drijdyj=-1.d0*drijdyi
        drijdzj=-1.d0*drijdzi
        drijdxk=0.0d0
        drijdyk=0.0d0
        drijdzk=0.0d0

        fcutij=0.5d0*(dcos(pi*rij/cutoff)+1.d0)
        temp1=0.5d0*(-dsin(pi*rij/cutoff))*(pi/cutoff)
        dfcutijdxi=temp1*drijdxi
        dfcutijdyi=temp1*drijdyi
        dfcutijdzi=temp1*drijdzi
        dfcutijdxj=-1.d0*dfcutijdxi
        dfcutijdyj=-1.d0*dfcutijdyi
        dfcutijdzj=-1.d0*dfcutijdzi
        dfcutijdxk=0.0d0
        dfcutijdyk=0.0d0
        dfcutijdzk=0.0d0

        do k=lsta(1,i),lsta(2,i) ! loop over all neighbors of this atom
          m=lstc(k)

!          if(j.ne.k) then ! neighbors j and k must be different
          if(k.gt.j) then ! neighbors j and k must be different

          rik=lstb(k,4)

          if(rik.le.cutoff)then

          deltaxk=-1.d0*(xyz(i,1)-lstb(k,1))
          deltayk=-1.d0*(xyz(i,2)-lstb(k,2))
          deltazk=-1.d0*(xyz(i,3)-lstb(k,3))
          drikdxi=(xyz(i,1)-lstb(k,1))/rik
          drikdyi=(xyz(i,2)-lstb(k,2))/rik
          drikdzi=(xyz(i,3)-lstb(k,3))/rik
          drikdxk=-1.d0*drikdxi
          drikdyk=-1.d0*drikdyi
          drikdzk=-1.d0*drikdzi
          drikdxj=0.0d0
          drikdyj=0.0d0
          drikdzj=0.0d0

          fcutik=0.5d0*(dcos(pi*rik/cutoff)+1.d0)
          temp1=0.5d0*(-dsin(pi*rik/cutoff))*(pi/cutoff)
          dfcutikdxi=temp1*drikdxi
          dfcutikdyi=temp1*drikdyi
          dfcutikdzi=temp1*drikdzi
          dfcutikdxj=0.0d0
          dfcutikdyj=0.0d0
          dfcutikdzj=0.0d0
          dfcutikdxk=-1.d0*dfcutikdxi
          dfcutikdyk=-1.d0*dfcutikdyi
          dfcutikdzk=-1.d0*dfcutikdzi

          rjk=(lstb(j,1)-lstb(k,1))**2 +&
              (lstb(j,2)-lstb(k,2))**2 +&
              (lstb(j,3)-lstb(k,3))**2
          rjk=dsqrt(rjk)

          if(rjk.le.cutoff)then

          drjkdxj=(lstb(j,1)-lstb(k,1))/rjk
          drjkdyj=(lstb(j,2)-lstb(k,2))/rjk
          drjkdzj=(lstb(j,3)-lstb(k,3))/rjk
          drjkdxk=-1.d0*drjkdxj
          drjkdyk=-1.d0*drjkdyj
          drjkdzk=-1.d0*drjkdzj
          drjkdxi=0.0d0
          drjkdyi=0.0d0
          drjkdzi=0.0d0

          fcutjk=0.5d0*(dcos(pi*rjk/cutoff)+1.d0)
          temp1=0.5d0*(-dsin(pi*rjk/cutoff))*(pi/cutoff)
          dfcutjkdxj=temp1*drjkdxj
          dfcutjkdyj=temp1*drjkdyj
          dfcutjkdzj=temp1*drjkdzj
          dfcutjkdxk=-1.d0*dfcutjkdxj
          dfcutjkdyk=-1.d0*dfcutjkdyj
          dfcutjkdzk=-1.d0*dfcutjkdzj
          dfcutjkdxi=0.0d0
          dfcutjkdyi=0.0d0
          dfcutjkdzi=0.0d0

! costheta=(rjk**2 - rij**2 -rik**2)/(-2.d0*rij*rik)
          f=rjk**2 - rij**2 -rik**2
          g=-2.d0*rij*rik
          costheta=f/g
          costheta=1.d0-costheta ! avoid negative values
          costheta=0.5d0*(costheta**2)
! summary: costheta=0.5*(1-f/g)**2

! calculate the derivatives of costheta
! (f/g)' = (f'g - fg')/g^2

            if(lforce)then
               dfdxi=-2.d0*rij*drijdxi - 2.d0*rik*drikdxi
               dfdyi=-2.d0*rij*drijdyi - 2.d0*rik*drikdyi
               dfdzi=-2.d0*rij*drijdzi - 2.d0*rik*drikdzi

               dfdxj=2.d0*rjk*drjkdxj - 2.d0*rij*drijdxj
               dfdyj=2.d0*rjk*drjkdyj - 2.d0*rij*drijdyj
               dfdzj=2.d0*rjk*drjkdzj - 2.d0*rij*drijdzj

               dfdxk=2.d0*rjk*drjkdxk - 2.d0*rik*drikdxk
               dfdyk=2.d0*rjk*drjkdyk - 2.d0*rik*drikdyk
               dfdzk=2.d0*rjk*drjkdzk - 2.d0*rik*drikdzk

               dgdxi=-2.d0*(drijdxi*rik + rij*drikdxi)
               dgdyi=-2.d0*(drijdyi*rik + rij*drikdyi)
               dgdzi=-2.d0*(drijdzi*rik + rij*drikdzi)

               dgdxj=-2.d0*drijdxj*rik
               dgdyj=-2.d0*drijdyj*rik
               dgdzj=-2.d0*drijdzj*rik

               dgdxk=-2.d0*rij*drikdxk
               dgdyk=-2.d0*rij*drikdyk
               dgdzk=-2.d0*rij*drikdzk

               dcosthetadxi=-(1.d0-f/g)*(dfdxi*g - f*dgdxi)/g**2
               dcosthetadyi=-(1.d0-f/g)*(dfdyi*g - f*dgdyi)/g**2
               dcosthetadzi=-(1.d0-f/g)*(dfdzi*g - f*dgdzi)/g**2
               dcosthetadxj=-(1.d0-f/g)*(dfdxj*g - f*dgdxj)/g**2
               dcosthetadyj=-(1.d0-f/g)*(dfdyj*g - f*dgdyj)/g**2
               dcosthetadzj=-(1.d0-f/g)*(dfdzj*g - f*dgdzj)/g**2
               dcosthetadxk=-(1.d0-f/g)*(dfdxk*g - f*dgdxk)/g**2
               dcosthetadyk=-(1.d0-f/g)*(dfdyk*g - f*dgdyk)/g**2
               dcosthetadzk=-(1.d0-f/g)*(dfdzk*g - f*dgdzk)/g**2
            endif


! calculation of exponential term
      expxyz=dexp(-alpha*(rij**2+rik**2+rjk**2))
      temp1=-alpha*2.0d0*expxyz
      dexpxyzdxi=(rij*drijdxi+rik*drikdxi+rjk*drjkdxi)*temp1
      dexpxyzdyi=(rij*drijdyi+rik*drikdyi+rjk*drjkdyi)*temp1
      dexpxyzdzi=(rij*drijdzi+rik*drikdzi+rjk*drjkdzi)*temp1
      dexpxyzdxj=(rij*drijdxj+rik*drikdxj+rjk*drjkdxj)*temp1
      dexpxyzdyj=(rij*drijdyj+rik*drikdyj+rjk*drjkdyj)*temp1
      dexpxyzdzj=(rij*drijdzj+rik*drikdzj+rjk*drjkdzj)*temp1
      dexpxyzdxk=(rij*drijdxk+rik*drikdxk+rjk*drjkdxk)*temp1
      dexpxyzdyk=(rij*drijdyk+rik*drikdyk+rjk*drjkdyk)*temp1
      dexpxyzdzk=(rij*drijdzk+rik*drikdzk+rjk*drjkdzk)*temp1

      xx(ii,1,i)=xx(ii,1,i)+costheta*expxyz*fcutij*fcutik*fcutjk

! (fgh)' = f'gh + fg'h + fgh'
! (fghi)'= f'ghi + fg'hi + fgh'i + fghi'

      if(lforce)then
! for the derivatives fcutjk is a constant
        temp1=(+dcosthetadxi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxi)
        temp2=(+dcosthetadxj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxj)
        temp3=(+dcosthetadxk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxk)
        dxdy(ii,i,i,1)=dxdy(ii,i,i,1)+temp1
        dxdy(ii,i,n,1)=dxdy(ii,i,n,1)+temp2
        dxdy(ii,i,m,1)=dxdy(ii,i,m,1)+temp3
        strs(1,1,ii,i)=strs(1,1,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,1,ii,i)=strs(2,1,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,1,ii,i)=strs(3,1,ii,i)+deltazj*temp2&
              +deltazk*temp3

        temp1=(+dcosthetadyi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyi)
        temp2=(+dcosthetadyj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyj)
        temp3=(+dcosthetadyk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyk)
        dxdy(ii,i,i,2)=dxdy(ii,i,i,2)+temp1
        dxdy(ii,i,n,2)=dxdy(ii,i,n,2)+temp2
        dxdy(ii,i,m,2)=dxdy(ii,i,m,2)+temp3
        strs(1,2,ii,i)=strs(1,2,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,2,ii,i)=strs(2,2,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,2,ii,i)=strs(3,2,ii,i)+deltazj*temp2&
              +deltazk*temp3

        temp1=(+dcosthetadzi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzi)
        temp2=(+dcosthetadzj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzj)
        temp3=(+dcosthetadzk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzk)
        dxdy(ii,i,i,3)=dxdy(ii,i,i,3)+temp1
        dxdy(ii,i,n,3)=dxdy(ii,i,n,3)+temp2
        dxdy(ii,i,m,3)=dxdy(ii,i,m,3)+temp3
        strs(1,3,ii,i)=strs(1,3,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,3,ii,i)=strs(2,3,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,3,ii,i)=strs(3,3,ii,i)+deltazj*temp2&
              +deltazk*temp3
      endif ! lforce

            endif ! rjk .le. cutoff
          endif ! rik .le. cutoff
          endif ! j .ne. k
        enddo ! k
        endif ! rij .le. cutoff
      enddo ! j
      enddo ! i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THREE BODY TERM WITH GLOBAL CUTOFF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      elseif(ntype.eq.7) then ! three-body term with global cutoff
      do i=1,nat
        if(ltest) write(*,*)'calculating three body term for atom ',i

      do j=lsta(1,i),lsta(2,i) ! loop over all neighbors of this atom
        n=lstc(j)
        rij=lstb(j,4)

        if(rij.le.cutoff)then

        deltaxj=-1.d0*(xyz(i,1)-lstb(j,1))
        deltayj=-1.d0*(xyz(i,2)-lstb(j,2))
        deltazj=-1.d0*(xyz(i,3)-lstb(j,3))
        drijdxi=-deltaxj/rij
        drijdyi=-deltayj/rij
        drijdzi=-deltazj/rij
        drijdxj=-1.d0*drijdxi
        drijdyj=-1.d0*drijdyi
        drijdzj=-1.d0*drijdzi
        drijdxk=0.0d0
        drijdyk=0.0d0
        drijdzk=0.0d0

        fcutij=0.5d0*(dcos(pi*rij/cutoff)+1.d0)
        temp1=0.5d0*(-dsin(pi*rij/cutoff))*(pi/cutoff)
        dfcutijdxi=temp1*drijdxi
        dfcutijdyi=temp1*drijdyi
        dfcutijdzi=temp1*drijdzi
        dfcutijdxj=-1.d0*dfcutijdxi
        dfcutijdyj=-1.d0*dfcutijdyi
        dfcutijdzj=-1.d0*dfcutijdzi
        dfcutijdxk=0.0d0
        dfcutijdyk=0.0d0
        dfcutijdzk=0.0d0

        do k=lsta(1,i),lsta(2,i) ! loop over all neighbors of this atom
          m=lstc(k)

!          if(j.ne.k) then ! neighbors j and k must be different
          if(k.gt.j) then ! neighbors j and k must be different

          rik=lstb(k,4)

          if(rik.le.cutoff)then

          deltaxk=-1.d0*(xyz(i,1)-lstb(k,1))
          deltayk=-1.d0*(xyz(i,2)-lstb(k,2))
          deltazk=-1.d0*(xyz(i,3)-lstb(k,3))
          drikdxi=(xyz(i,1)-lstb(k,1))/rik
          drikdyi=(xyz(i,2)-lstb(k,2))/rik
          drikdzi=(xyz(i,3)-lstb(k,3))/rik
          drikdxk=-1.d0*drikdxi
          drikdyk=-1.d0*drikdyi
          drikdzk=-1.d0*drikdzi
          drikdxj=0.0d0
          drikdyj=0.0d0
          drikdzj=0.0d0

          fcutik=0.5d0*(dcos(pi*rik/cutoff)+1.d0)
          temp1=0.5d0*(-dsin(pi*rik/cutoff))*(pi/cutoff)
          dfcutikdxi=temp1*drikdxi
          dfcutikdyi=temp1*drikdyi
          dfcutikdzi=temp1*drikdzi
          dfcutikdxj=0.0d0
          dfcutikdyj=0.0d0
          dfcutikdzj=0.0d0
          dfcutikdxk=-1.d0*dfcutikdxi
          dfcutikdyk=-1.d0*dfcutikdyi
          dfcutikdzk=-1.d0*dfcutikdzi

          rjk=(lstb(j,1)-lstb(k,1))**2 +&
              (lstb(j,2)-lstb(k,2))**2 +&
              (lstb(j,3)-lstb(k,3))**2
          rjk=dsqrt(rjk)
          if(rjk.le.0.5d0) then
            write(*,*)'Error rjk .le. 0.5',ntype
!            stop
          endif

          if(rjk.le.cutoff)then

          drjkdxj=(lstb(j,1)-lstb(k,1))/rjk
          drjkdyj=(lstb(j,2)-lstb(k,2))/rjk
          drjkdzj=(lstb(j,3)-lstb(k,3))/rjk
          drjkdxk=-1.d0*drjkdxj
          drjkdyk=-1.d0*drjkdyj
          drjkdzk=-1.d0*drjkdzj
          drjkdxi=0.0d0
          drjkdyi=0.0d0
          drjkdzi=0.0d0

          fcutjk=0.5d0*(dcos(pi*rjk/cutoff)+1.d0)
          temp1=0.5d0*(-dsin(pi*rjk/cutoff))*(pi/cutoff)
          dfcutjkdxj=temp1*drjkdxj
          dfcutjkdyj=temp1*drjkdyj
          dfcutjkdzj=temp1*drjkdzj
          dfcutjkdxk=-1.d0*dfcutjkdxj
          dfcutjkdyk=-1.d0*dfcutjkdyj
          dfcutjkdzk=-1.d0*dfcutjkdzj
          dfcutjkdxi=0.0d0
          dfcutjkdyi=0.0d0
          dfcutjkdzi=0.0d0

! costheta=(rjk**2 - rij**2 -rik**2)/(-2.d0*rij*rik)
            f=rjk**2 - rij**2 -rik**2
            g=-2.d0*rij*rik
            costheta=f/g
            costheta=costheta+1.d0 ! avoid negative values
            costheta=0.125d0*(costheta**4)
! summary: costheta=1/8*(1+f/g)**4

! calculate the derivatives of costheta
! (f/g)' = (f'g - fg')/g^2

            if(lforce)then
               dfdxi=-2.d0*rij*drijdxi - 2.d0*rik*drikdxi
               dfdyi=-2.d0*rij*drijdyi - 2.d0*rik*drikdyi
               dfdzi=-2.d0*rij*drijdzi - 2.d0*rik*drikdzi

               dfdxj=2.d0*rjk*drjkdxj - 2.d0*rij*drijdxj
               dfdyj=2.d0*rjk*drjkdyj - 2.d0*rij*drijdyj
               dfdzj=2.d0*rjk*drjkdzj - 2.d0*rij*drijdzj

               dfdxk=2.d0*rjk*drjkdxk - 2.d0*rik*drikdxk
               dfdyk=2.d0*rjk*drjkdyk - 2.d0*rik*drikdyk
               dfdzk=2.d0*rjk*drjkdzk - 2.d0*rik*drikdzk

               dgdxi=-2.d0*(drijdxi*rik + rij*drikdxi)
               dgdyi=-2.d0*(drijdyi*rik + rij*drikdyi)
               dgdzi=-2.d0*(drijdzi*rik + rij*drikdzi)

               dgdxj=-2.d0*drijdxj*rik
               dgdyj=-2.d0*drijdyj*rik
               dgdzj=-2.d0*drijdzj*rik

               dgdxk=-2.d0*rij*drikdxk
               dgdyk=-2.d0*rij*drikdyk
               dgdzk=-2.d0*rij*drikdzk

               dcosthetadxi=0.5d0*(1.d0+f/g)**3*(dfdxi*g-f*dgdxi)/g**2
               dcosthetadyi=0.5d0*(1.d0+f/g)**3*(dfdyi*g-f*dgdyi)/g**2
               dcosthetadzi=0.5d0*(1.d0+f/g)**3*(dfdzi*g-f*dgdzi)/g**2
               dcosthetadxj=0.5d0*(1.d0+f/g)**3*(dfdxj*g-f*dgdxj)/g**2
               dcosthetadyj=0.5d0*(1.d0+f/g)**3*(dfdyj*g-f*dgdyj)/g**2
               dcosthetadzj=0.5d0*(1.d0+f/g)**3*(dfdzj*g-f*dgdzj)/g**2
               dcosthetadxk=0.5d0*(1.d0+f/g)**3*(dfdxk*g-f*dgdxk)/g**2
               dcosthetadyk=0.5d0*(1.d0+f/g)**3*(dfdyk*g-f*dgdyk)/g**2
               dcosthetadzk=0.5d0*(1.d0+f/g)**3*(dfdzk*g-f*dgdzk)/g**2
            endif

! calculation of exponential term
      expxyz=dexp(-alpha*(rij**2+rik**2+rjk**2))
      temp1=-alpha*2.0d0*expxyz
      dexpxyzdxi=(rij*drijdxi+rik*drikdxi+rjk*drjkdxi)*temp1
      dexpxyzdyi=(rij*drijdyi+rik*drikdyi+rjk*drjkdyi)*temp1
      dexpxyzdzi=(rij*drijdzi+rik*drikdzi+rjk*drjkdzi)*temp1
      dexpxyzdxj=(rij*drijdxj+rik*drikdxj+rjk*drjkdxj)*temp1
      dexpxyzdyj=(rij*drijdyj+rik*drikdyj+rjk*drjkdyj)*temp1
      dexpxyzdzj=(rij*drijdzj+rik*drikdzj+rjk*drjkdzj)*temp1
      dexpxyzdxk=(rij*drijdxk+rik*drikdxk+rjk*drjkdxk)*temp1
      dexpxyzdyk=(rij*drijdyk+rik*drikdyk+rjk*drjkdyk)*temp1
      dexpxyzdzk=(rij*drijdzk+rik*drikdzk+rjk*drjkdzk)*temp1

      xx(ii,1,i)=xx(ii,1,i)+costheta*expxyz*fcutij*fcutik*fcutjk

! (fgh)' = f'gh + fg'h + fgh'
! (fghi)'= f'ghi + fg'hi + fgh'i + fghi'

      if(lforce)then
! for the derivatives fcutjk is a constant
        temp1=(+dcosthetadxi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxi)
        temp2=(+dcosthetadxj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxj)
        temp3=(+dcosthetadxk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxk)
        dxdy(ii,i,i,1)=dxdy(ii,i,i,1)+temp1
        dxdy(ii,i,n,1)=dxdy(ii,i,n,1)+temp2
        dxdy(ii,i,m,1)=dxdy(ii,i,m,1)+temp3
        strs(1,1,ii,i)=strs(1,1,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,1,ii,i)=strs(2,1,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,1,ii,i)=strs(3,1,ii,i)+deltazj*temp2&
              +deltazk*temp3

        temp1=(+dcosthetadyi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyi)
        temp2=(+dcosthetadyj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyj)
        temp3=(+dcosthetadyk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyk)
        dxdy(ii,i,i,2)=dxdy(ii,i,i,2)+temp1
        dxdy(ii,i,n,2)=dxdy(ii,i,n,2)+temp2
        dxdy(ii,i,m,2)=dxdy(ii,i,m,2)+temp3
        strs(1,2,ii,i)=strs(1,2,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,2,ii,i)=strs(2,2,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,2,ii,i)=strs(3,2,ii,i)+deltazj*temp2&
              +deltazk*temp3

        temp1=(+dcosthetadzi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzi)
        temp2=(+dcosthetadzj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzj)
        temp3=(+dcosthetadzk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzk)
        dxdy(ii,i,i,3)=dxdy(ii,i,i,3)+temp1
        dxdy(ii,i,n,3)=dxdy(ii,i,n,3)+temp2
        dxdy(ii,i,m,3)=dxdy(ii,i,m,3)+temp3
        strs(1,3,ii,i)=strs(1,3,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,3,ii,i)=strs(2,3,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,3,ii,i)=strs(3,3,ii,i)+deltazj*temp2&
              +deltazk*temp3
      endif ! lforce

            endif ! rjk .le. cutoff
          endif ! rik .le. cutoff
          endif ! j .ne. k
        enddo ! k
        endif ! rij .le. cutoff
      enddo ! j
      enddo ! i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THREE BODY TERM WITH GLOBAL CUTOFF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      elseif(ntype.eq.8) then ! three-body term with global cutoff
      do i=1,nat
        if(ltest) write(*,*)'calculating three body term for atom ',i

      do j=lsta(1,i),lsta(2,i) ! loop over all neighbors of this atom
        n=lstc(j)
        rij=lstb(j,4)

        if(rij.le.cutoff)then

        deltaxj=-1.d0*(xyz(i,1)-lstb(j,1))
        deltayj=-1.d0*(xyz(i,2)-lstb(j,2))
        deltazj=-1.d0*(xyz(i,3)-lstb(j,3))
        drijdxi=-deltaxj/rij
        drijdyi=-deltayj/rij
        drijdzi=-deltazj/rij
        drijdxj=-1.d0*drijdxi
        drijdyj=-1.d0*drijdyi
        drijdzj=-1.d0*drijdzi
        drijdxk=0.0d0
        drijdyk=0.0d0
        drijdzk=0.0d0

        fcutij=0.5d0*(dcos(pi*rij/cutoff)+1.d0)
        temp1=0.5d0*(-dsin(pi*rij/cutoff))*(pi/cutoff)
        dfcutijdxi=temp1*drijdxi
        dfcutijdyi=temp1*drijdyi
        dfcutijdzi=temp1*drijdzi
        dfcutijdxj=-1.d0*dfcutijdxi
        dfcutijdyj=-1.d0*dfcutijdyi
        dfcutijdzj=-1.d0*dfcutijdzi
        dfcutijdxk=0.0d0
        dfcutijdyk=0.0d0
        dfcutijdzk=0.0d0

        do k=lsta(1,i),lsta(2,i) ! loop over all neighbors of this atom
          m=lstc(k)

!          if(j.ne.k) then ! neighbors j and k must be different
          if(k.gt.j) then ! neighbors j and k must be different

          rik=lstb(k,4)

          if(rik.le.cutoff)then

          deltaxk=-1.d0*(xyz(i,1)-lstb(k,1))
          deltayk=-1.d0*(xyz(i,2)-lstb(k,2))
          deltazk=-1.d0*(xyz(i,3)-lstb(k,3))
          drikdxi=(xyz(i,1)-lstb(k,1))/rik
          drikdyi=(xyz(i,2)-lstb(k,2))/rik
          drikdzi=(xyz(i,3)-lstb(k,3))/rik
          drikdxk=-1.d0*drikdxi
          drikdyk=-1.d0*drikdyi
          drikdzk=-1.d0*drikdzi
          drikdxj=0.0d0
          drikdyj=0.0d0
          drikdzj=0.0d0

          fcutik=0.5d0*(dcos(pi*rik/cutoff)+1.d0)
          temp1=0.5d0*(-dsin(pi*rik/cutoff))*(pi/cutoff)
          dfcutikdxi=temp1*drikdxi
          dfcutikdyi=temp1*drikdyi
          dfcutikdzi=temp1*drikdzi
          dfcutikdxj=0.0d0
          dfcutikdyj=0.0d0
          dfcutikdzj=0.0d0
          dfcutikdxk=-1.d0*dfcutikdxi
          dfcutikdyk=-1.d0*dfcutikdyi
          dfcutikdzk=-1.d0*dfcutikdzi

          rjk=(lstb(j,1)-lstb(k,1))**2 +&
              (lstb(j,2)-lstb(k,2))**2 +&
              (lstb(j,3)-lstb(k,3))**2
          rjk=dsqrt(rjk)

          if(rjk.le.cutoff)then

          drjkdxj=(lstb(j,1)-lstb(k,1))/rjk
          drjkdyj=(lstb(j,2)-lstb(k,2))/rjk
          drjkdzj=(lstb(j,3)-lstb(k,3))/rjk
          drjkdxk=-1.d0*drjkdxj
          drjkdyk=-1.d0*drjkdyj
          drjkdzk=-1.d0*drjkdzj
          drjkdxi=0.0d0
          drjkdyi=0.0d0
          drjkdzi=0.0d0

          fcutjk=0.5d0*(dcos(pi*rjk/cutoff)+1.d0)
          temp1=0.5d0*(-dsin(pi*rjk/cutoff))*(pi/cutoff)
          dfcutjkdxj=temp1*drjkdxj
          dfcutjkdyj=temp1*drjkdyj
          dfcutjkdzj=temp1*drjkdzj
          dfcutjkdxk=-1.d0*dfcutjkdxj
          dfcutjkdyk=-1.d0*dfcutjkdyj
          dfcutjkdzk=-1.d0*dfcutjkdzj
          dfcutjkdxi=0.0d0
          dfcutjkdyi=0.0d0
          dfcutjkdzi=0.0d0

! costheta=(rjk**2 - rij**2 -rik**2)/(-2.d0*rij*rik)
          f=rjk**2 - rij**2 -rik**2
          g=-2.d0*rij*rik
          costheta=f/g
          costheta=1.d0-costheta ! avoid negative values
          costheta=0.125d0*(costheta**4)
! summary: costheta=1/8*(1-f/g)**4

! calculate the derivatives of costheta
! (f/g)' = (f'g - fg')/g^2

            if(lforce) then
               dfdxi=-2.d0*rij*drijdxi - 2.d0*rik*drikdxi
               dfdyi=-2.d0*rij*drijdyi - 2.d0*rik*drikdyi
               dfdzi=-2.d0*rij*drijdzi - 2.d0*rik*drikdzi

               dfdxj=2.d0*rjk*drjkdxj - 2.d0*rij*drijdxj
               dfdyj=2.d0*rjk*drjkdyj - 2.d0*rij*drijdyj
               dfdzj=2.d0*rjk*drjkdzj - 2.d0*rij*drijdzj

               dfdxk=2.d0*rjk*drjkdxk - 2.d0*rik*drikdxk
               dfdyk=2.d0*rjk*drjkdyk - 2.d0*rik*drikdyk
               dfdzk=2.d0*rjk*drjkdzk - 2.d0*rik*drikdzk

               dgdxi=-2.d0*(drijdxi*rik + rij*drikdxi)
               dgdyi=-2.d0*(drijdyi*rik + rij*drikdyi)
               dgdzi=-2.d0*(drijdzi*rik + rij*drikdzi)

               dgdxj=-2.d0*drijdxj*rik
               dgdyj=-2.d0*drijdyj*rik
               dgdzj=-2.d0*drijdzj*rik

               dgdxk=-2.d0*rij*drikdxk
               dgdyk=-2.d0*rij*drikdyk
               dgdzk=-2.d0*rij*drikdzk

               dcosthetadxi=-0.5d0*(1.d0-f/g)**3*(dfdxi*g-f*dgdxi)/g**2
               dcosthetadyi=-0.5d0*(1.d0-f/g)**3*(dfdyi*g-f*dgdyi)/g**2
               dcosthetadzi=-0.5d0*(1.d0-f/g)**3*(dfdzi*g-f*dgdzi)/g**2
               dcosthetadxj=-0.5d0*(1.d0-f/g)**3*(dfdxj*g-f*dgdxj)/g**2
               dcosthetadyj=-0.5d0*(1.d0-f/g)**3*(dfdyj*g-f*dgdyj)/g**2
               dcosthetadzj=-0.5d0*(1.d0-f/g)**3*(dfdzj*g-f*dgdzj)/g**2
               dcosthetadxk=-0.5d0*(1.d0-f/g)**3*(dfdxk*g-f*dgdxk)/g**2
               dcosthetadyk=-0.5d0*(1.d0-f/g)**3*(dfdyk*g-f*dgdyk)/g**2
               dcosthetadzk=-0.5d0*(1.d0-f/g)**3*(dfdzk*g-f*dgdzk)/g**2
            endif


! calculation of exponential term
! CAUTION: Is there a square missing here?
      expxyz=dexp(-alpha*(rij**2+rik**2+rjk**2))
      temp1=-alpha*2.0d0*expxyz
      dexpxyzdxi=(rij*drijdxi+rik*drikdxi+rjk*drjkdxi)*temp1
      dexpxyzdyi=(rij*drijdyi+rik*drikdyi+rjk*drjkdyi)*temp1
      dexpxyzdzi=(rij*drijdzi+rik*drikdzi+rjk*drjkdzi)*temp1
      dexpxyzdxj=(rij*drijdxj+rik*drikdxj+rjk*drjkdxj)*temp1
      dexpxyzdyj=(rij*drijdyj+rik*drikdyj+rjk*drjkdyj)*temp1
      dexpxyzdzj=(rij*drijdzj+rik*drikdzj+rjk*drjkdzj)*temp1
      dexpxyzdxk=(rij*drijdxk+rik*drikdxk+rjk*drjkdxk)*temp1
      dexpxyzdyk=(rij*drijdyk+rik*drikdyk+rjk*drjkdyk)*temp1
      dexpxyzdzk=(rij*drijdzk+rik*drikdzk+rjk*drjkdzk)*temp1

      xx(ii,1,i)=xx(ii,1,i)+costheta*expxyz*fcutij*fcutik*fcutjk

! (fgh)' = f'gh + fg'h + fgh'
! (fghi)'= f'ghi + fg'hi + fgh'i + fghi'

      if(lforce)then
! for the derivatives fcutjk is a constant
        temp1=(+dcosthetadxi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxi)
        temp2=(+dcosthetadxj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxj)
        temp3=(+dcosthetadxk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxk)
        dxdy(ii,i,i,1)=dxdy(ii,i,i,1)+temp1
        dxdy(ii,i,n,1)=dxdy(ii,i,n,1)+temp2
        dxdy(ii,i,m,1)=dxdy(ii,i,m,1)+temp3
        strs(1,1,ii,i)=strs(1,1,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,1,ii,i)=strs(2,1,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,1,ii,i)=strs(3,1,ii,i)+deltazj*temp2&
              +deltazk*temp3

        temp1=(+dcosthetadyi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyi)
        temp2=(+dcosthetadyj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyj)
        temp3=(+dcosthetadyk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyk)
        dxdy(ii,i,i,2)=dxdy(ii,i,i,2)+temp1
        dxdy(ii,i,n,2)=dxdy(ii,i,n,2)+temp2
        dxdy(ii,i,m,2)=dxdy(ii,i,m,2)+temp3
        strs(1,2,ii,i)=strs(1,2,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,2,ii,i)=strs(2,2,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,2,ii,i)=strs(3,2,ii,i)+deltazj*temp2&
              +deltazk*temp3

        temp1=(+dcosthetadzi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzi)
        temp2=(+dcosthetadzj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzj)
        temp3=(+dcosthetadzk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzk)
        dxdy(ii,i,i,3)=dxdy(ii,i,i,3)+temp1
        dxdy(ii,i,n,3)=dxdy(ii,i,n,3)+temp2
        dxdy(ii,i,m,3)=dxdy(ii,i,m,3)+temp3
        strs(1,3,ii,i)=strs(1,3,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,3,ii,i)=strs(2,3,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,3,ii,i)=strs(3,3,ii,i)+deltazj*temp2&
              +deltazk*temp3
      endif ! lforce

            endif ! rjk .le. cutoff
          endif ! rik .le. cutoff
          endif ! j .ne. k
        enddo ! k
        endif ! rij .le. cutoff
      enddo ! j
      enddo ! i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THREE BODY TERM WITH GLOBAL CUTOFF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      elseif(ntype.eq.9) then ! three-body term with global cutoff
      do i=1,nat
        if(ltest) write(*,*)'calculating three body term for atom ',i

      do j=lsta(1,i),lsta(2,i) ! loop over all neighbors of this atom
        n=lstc(j)
        rij=lstb(j,4)

        if(rij.le.cutoff)then

        deltaxj=-1.d0*(xyz(i,1)-lstb(j,1))
        deltayj=-1.d0*(xyz(i,2)-lstb(j,2))
        deltazj=-1.d0*(xyz(i,3)-lstb(j,3))
        drijdxi=-deltaxj/rij
        drijdyi=-deltayj/rij
        drijdzi=-deltazj/rij
        drijdxj=-1.d0*drijdxi
        drijdyj=-1.d0*drijdyi
        drijdzj=-1.d0*drijdzi
        drijdxk=0.0d0
        drijdyk=0.0d0
        drijdzk=0.0d0

        fcutij=0.5d0*(dcos(pi*rij/cutoff)+1.d0)
        temp1=0.5d0*(-dsin(pi*rij/cutoff))*(pi/cutoff)
        dfcutijdxi=temp1*drijdxi
        dfcutijdyi=temp1*drijdyi
        dfcutijdzi=temp1*drijdzi
        dfcutijdxj=-1.d0*dfcutijdxi
        dfcutijdyj=-1.d0*dfcutijdyi
        dfcutijdzj=-1.d0*dfcutijdzi
        dfcutijdxk=0.0d0
        dfcutijdyk=0.0d0
        dfcutijdzk=0.0d0

        do k=lsta(1,i),lsta(2,i) ! loop over all neighbors of this atom
          m=lstc(k)

!          if(j.ne.k) then ! neighbors j and k must be different
          if(k.gt.j) then ! neighbors j and k must be different

          rik=lstb(k,4)

          if(rik.le.cutoff)then
          deltaxk=-1.d0*(xyz(i,1)-lstb(k,1))
          deltayk=-1.d0*(xyz(i,2)-lstb(k,2))
          deltazk=-1.d0*(xyz(i,3)-lstb(k,3))
          drikdxi=(xyz(i,1)-lstb(k,1))/rik
          drikdyi=(xyz(i,2)-lstb(k,2))/rik
          drikdzi=(xyz(i,3)-lstb(k,3))/rik
          drikdxk=-1.d0*drikdxi
          drikdyk=-1.d0*drikdyi
          drikdzk=-1.d0*drikdzi
          drikdxj=0.0d0
          drikdyj=0.0d0
          drikdzj=0.0d0

          fcutik=0.5d0*(dcos(pi*rik/cutoff)+1.d0)
          temp1=0.5d0*(-dsin(pi*rik/cutoff))*(pi/cutoff)
          dfcutikdxi=temp1*drikdxi
          dfcutikdyi=temp1*drikdyi
          dfcutikdzi=temp1*drikdzi
          dfcutikdxj=0.0d0
          dfcutikdyj=0.0d0
          dfcutikdzj=0.0d0
          dfcutikdxk=-1.d0*dfcutikdxi
          dfcutikdyk=-1.d0*dfcutikdyi
          dfcutikdzk=-1.d0*dfcutikdzi

          rjk=(lstb(j,1)-lstb(k,1))**2 +&
              (lstb(j,2)-lstb(k,2))**2 +&
              (lstb(j,3)-lstb(k,3))**2
          rjk=dsqrt(rjk)

          if(rjk.le.cutoff)then

          drjkdxj=(lstb(j,1)-lstb(k,1))/rjk
          drjkdyj=(lstb(j,2)-lstb(k,2))/rjk
          drjkdzj=(lstb(j,3)-lstb(k,3))/rjk
          drjkdxk=-1.d0*drjkdxj
          drjkdyk=-1.d0*drjkdyj
          drjkdzk=-1.d0*drjkdzj
          drjkdxi=0.0d0
          drjkdyi=0.0d0
          drjkdzi=0.0d0

          fcutjk=0.5d0*(dcos(pi*rjk/cutoff)+1.d0)
          temp1=0.5d0*(-dsin(pi*rjk/cutoff))*(pi/cutoff)
          dfcutjkdxj=temp1*drjkdxj
          dfcutjkdyj=temp1*drjkdyj
          dfcutjkdzj=temp1*drjkdzj
          dfcutjkdxk=-1.d0*dfcutjkdxj
          dfcutjkdyk=-1.d0*dfcutjkdyj
          dfcutjkdzk=-1.d0*dfcutjkdzj
          dfcutjkdxi=0.0d0
          dfcutjkdyi=0.0d0
          dfcutjkdzi=0.0d0

! costheta=(rjk**2 - rij**2 -rik**2)/(-2.d0*rij*rik)
            f=rjk**2 - rij**2 -rik**2
            g=-2.d0*rij*rik
            costheta=f/g
            costheta=costheta+1.d0 ! avoid negative values
            costheta=(1.d0/32768.d0)*(costheta**16)
! summary: costheta=1/32768*(1+f/g)**16

! calculate the derivatives of costheta
! (f/g)' = (f'g - fg')/g^2

            if(lforce) then
               dfdxi=-2.d0*rij*drijdxi - 2.d0*rik*drikdxi
               dfdyi=-2.d0*rij*drijdyi - 2.d0*rik*drikdyi
               dfdzi=-2.d0*rij*drijdzi - 2.d0*rik*drikdzi

               dfdxj=2.d0*rjk*drjkdxj - 2.d0*rij*drijdxj
               dfdyj=2.d0*rjk*drjkdyj - 2.d0*rij*drijdyj
               dfdzj=2.d0*rjk*drjkdzj - 2.d0*rij*drijdzj

               dfdxk=2.d0*rjk*drjkdxk - 2.d0*rik*drikdxk
               dfdyk=2.d0*rjk*drjkdyk - 2.d0*rik*drikdyk
               dfdzk=2.d0*rjk*drjkdzk - 2.d0*rik*drikdzk

               dgdxi=-2.d0*(drijdxi*rik + rij*drikdxi)
               dgdyi=-2.d0*(drijdyi*rik + rij*drikdyi)
               dgdzi=-2.d0*(drijdzi*rik + rij*drikdzi)

               dgdxj=-2.d0*drijdxj*rik
               dgdyj=-2.d0*drijdyj*rik
               dgdzj=-2.d0*drijdzj*rik

               dgdxk=-2.d0*rij*drikdxk
               dgdyk=-2.d0*rij*drikdyk
               dgdzk=-2.d0*rij*drikdzk

       dcosthetadxi=(1.d0/2048.d0)*(1.d0+f/g)**15*(dfdxi*g-f*dgdxi)/g**2
       dcosthetadyi=(1.d0/2048.d0)*(1.d0+f/g)**15*(dfdyi*g-f*dgdyi)/g**2
       dcosthetadzi=(1.d0/2048.d0)*(1.d0+f/g)**15*(dfdzi*g-f*dgdzi)/g**2
       dcosthetadxj=(1.d0/2048.d0)*(1.d0+f/g)**15*(dfdxj*g-f*dgdxj)/g**2
       dcosthetadyj=(1.d0/2048.d0)*(1.d0+f/g)**15*(dfdyj*g-f*dgdyj)/g**2
       dcosthetadzj=(1.d0/2048.d0)*(1.d0+f/g)**15*(dfdzj*g-f*dgdzj)/g**2
       dcosthetadxk=(1.d0/2048.d0)*(1.d0+f/g)**15*(dfdxk*g-f*dgdxk)/g**2
       dcosthetadyk=(1.d0/2048.d0)*(1.d0+f/g)**15*(dfdyk*g-f*dgdyk)/g**2
       dcosthetadzk=(1.d0/2048.d0)*(1.d0+f/g)**15*(dfdzk*g-f*dgdzk)/g**2
           endif

! calculation of exponential term
      expxyz=dexp(-alpha*(rij**2+rik**2+rjk**2))
      temp1=-alpha*2.0d0*expxyz
      dexpxyzdxi=(rij*drijdxi+rik*drikdxi+rjk*drjkdxi)*temp1
      dexpxyzdyi=(rij*drijdyi+rik*drikdyi+rjk*drjkdyi)*temp1
      dexpxyzdzi=(rij*drijdzi+rik*drikdzi+rjk*drjkdzi)*temp1
      dexpxyzdxj=(rij*drijdxj+rik*drikdxj+rjk*drjkdxj)*temp1
      dexpxyzdyj=(rij*drijdyj+rik*drikdyj+rjk*drjkdyj)*temp1
      dexpxyzdzj=(rij*drijdzj+rik*drikdzj+rjk*drjkdzj)*temp1
      dexpxyzdxk=(rij*drijdxk+rik*drikdxk+rjk*drjkdxk)*temp1
      dexpxyzdyk=(rij*drijdyk+rik*drikdyk+rjk*drjkdyk)*temp1
      dexpxyzdzk=(rij*drijdzk+rik*drikdzk+rjk*drjkdzk)*temp1

      xx(ii,1,i)=xx(ii,1,i)+costheta*expxyz*fcutij*fcutik*fcutjk

! (fgh)' = f'gh + fg'h + fgh'
! (fghi)'= f'ghi + fg'hi + fgh'i + fghi'

      if(lforce)then
! for the derivatives fcutjk is a constant
        temp1=(+dcosthetadxi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxi)
        temp2=(+dcosthetadxj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxj)
        temp3=(+dcosthetadxk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxk)
        dxdy(ii,i,i,1)=dxdy(ii,i,i,1)+temp1
        dxdy(ii,i,n,1)=dxdy(ii,i,n,1)+temp2
        dxdy(ii,i,m,1)=dxdy(ii,i,m,1)+temp3
        strs(1,1,ii,i)=strs(1,1,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,1,ii,i)=strs(2,1,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,1,ii,i)=strs(3,1,ii,i)+deltazj*temp2&
              +deltazk*temp3

        temp1=(+dcosthetadyi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyi)
        temp2=(+dcosthetadyj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyj)
        temp3=(+dcosthetadyk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyk)
        dxdy(ii,i,i,2)=dxdy(ii,i,i,2)+temp1
        dxdy(ii,i,n,2)=dxdy(ii,i,n,2)+temp2
        dxdy(ii,i,m,2)=dxdy(ii,i,m,2)+temp3
        strs(1,2,ii,i)=strs(1,2,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,2,ii,i)=strs(2,2,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,2,ii,i)=strs(3,2,ii,i)+deltazj*temp2&
              +deltazk*temp3

        temp1=(+dcosthetadzi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzi)
        temp2=(+dcosthetadzj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzj)
        temp3=(+dcosthetadzk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzk)
        dxdy(ii,i,i,3)=dxdy(ii,i,i,3)+temp1
        dxdy(ii,i,n,3)=dxdy(ii,i,n,3)+temp2
        dxdy(ii,i,m,3)=dxdy(ii,i,m,3)+temp3
        strs(1,3,ii,i)=strs(1,3,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,3,ii,i)=strs(2,3,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,3,ii,i)=strs(3,3,ii,i)+deltazj*temp2&
              +deltazk*temp3
      endif !lforce

            endif ! rjk .le. cutoff
          endif ! rik .le. cutoff
          endif ! j .ne. k
        enddo ! k
        endif ! rij .le. cutoff
      enddo ! j
      enddo ! i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THREE BODY TERM WITH GLOBAL CUTOFF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      elseif(ntype.eq.10) then ! three-body term with global cutoff
      do i=1,nat
        if(ltest) write(*,*)'calculating three body term for atom ',i

      do j=lsta(1,i),lsta(2,i) ! loop over all neighbors of this atom
        n=lstc(j)
        rij=lstb(j,4)

        if(rij.le.cutoff)then

        deltaxj=-1.d0*(xyz(i,1)-lstb(j,1))
        deltayj=-1.d0*(xyz(i,2)-lstb(j,2))
        deltazj=-1.d0*(xyz(i,3)-lstb(j,3))
        drijdxi=-deltaxj/rij
        drijdyi=-deltayj/rij
        drijdzi=-deltazj/rij
        drijdxj=-1.d0*drijdxi
        drijdyj=-1.d0*drijdyi
        drijdzj=-1.d0*drijdzi
        drijdxk=0.0d0
        drijdyk=0.0d0
        drijdzk=0.0d0

        fcutij=0.5d0*(dcos(pi*rij/cutoff)+1.d0)
        temp1=0.5d0*(-dsin(pi*rij/cutoff))*(pi/cutoff)
        dfcutijdxi=temp1*drijdxi
        dfcutijdyi=temp1*drijdyi
        dfcutijdzi=temp1*drijdzi
        dfcutijdxj=-1.d0*dfcutijdxi
        dfcutijdyj=-1.d0*dfcutijdyi
        dfcutijdzj=-1.d0*dfcutijdzi
        dfcutijdxk=0.0d0
        dfcutijdyk=0.0d0
        dfcutijdzk=0.0d0

        do k=lsta(1,i),lsta(2,i) ! loop over all neighbors of this atom
          m=lstc(k)

!          if(j.ne.k) then ! neighbors j and k must be different
          if(k.gt.j) then ! neighbors j and k must be different

          rik=lstb(k,4)

          if(rik.le.cutoff)then
          deltaxk=-1.d0*(xyz(i,1)-lstb(k,1))
          deltayk=-1.d0*(xyz(i,2)-lstb(k,2))
          deltazk=-1.d0*(xyz(i,3)-lstb(k,3))
          drikdxi=(xyz(i,1)-lstb(k,1))/rik
          drikdyi=(xyz(i,2)-lstb(k,2))/rik
          drikdzi=(xyz(i,3)-lstb(k,3))/rik
          drikdxk=-1.d0*drikdxi
          drikdyk=-1.d0*drikdyi
          drikdzk=-1.d0*drikdzi
          drikdxj=0.0d0
          drikdyj=0.0d0
          drikdzj=0.0d0

          fcutik=0.5d0*(dcos(pi*rik/cutoff)+1.d0)
          temp1=0.5d0*(-dsin(pi*rik/cutoff))*(pi/cutoff)
          dfcutikdxi=temp1*drikdxi
          dfcutikdyi=temp1*drikdyi
          dfcutikdzi=temp1*drikdzi
          dfcutikdxj=0.0d0
          dfcutikdyj=0.0d0
          dfcutikdzj=0.0d0
          dfcutikdxk=-1.d0*dfcutikdxi
          dfcutikdyk=-1.d0*dfcutikdyi
          dfcutikdzk=-1.d0*dfcutikdzi

          rjk=(lstb(j,1)-lstb(k,1))**2 +&
              (lstb(j,2)-lstb(k,2))**2 +&
              (lstb(j,3)-lstb(k,3))**2
          rjk=dsqrt(rjk)

          if(rjk.le.cutoff)then

          drjkdxj=(lstb(j,1)-lstb(k,1))/rjk
          drjkdyj=(lstb(j,2)-lstb(k,2))/rjk
          drjkdzj=(lstb(j,3)-lstb(k,3))/rjk
          drjkdxk=-1.d0*drjkdxj
          drjkdyk=-1.d0*drjkdyj
          drjkdzk=-1.d0*drjkdzj
          drjkdxi=0.0d0
          drjkdyi=0.0d0
          drjkdzi=0.0d0

          fcutjk=0.5d0*(dcos(pi*rjk/cutoff)+1.d0)
          temp1=0.5d0*(-dsin(pi*rjk/cutoff))*(pi/cutoff)
          dfcutjkdxj=temp1*drjkdxj
          dfcutjkdyj=temp1*drjkdyj
          dfcutjkdzj=temp1*drjkdzj
          dfcutjkdxk=-1.d0*dfcutjkdxj
          dfcutjkdyk=-1.d0*dfcutjkdyj
          dfcutjkdzk=-1.d0*dfcutjkdzj
          dfcutjkdxi=0.0d0
          dfcutjkdyi=0.0d0
          dfcutjkdzi=0.0d0

! costheta=(rjk**2 - rij**2 -rik**2)/(-2.d0*rij*rik)
            f=rjk**2 - rij**2 -rik**2
            g=-2.d0*rij*rik
            costheta=f/g
            costheta=1.d0-costheta ! avoid negative values
            costheta=(1.d0/32768.d0)*(costheta**16)
! summary: costheta=1/32768*(1-f/g)**16

! calculate the derivatives of costheta
! (f/g)' = (f'g - fg')/g^2

            if(lforce)then
               dfdxi=-2.d0*rij*drijdxi - 2.d0*rik*drikdxi
               dfdyi=-2.d0*rij*drijdyi - 2.d0*rik*drikdyi
               dfdzi=-2.d0*rij*drijdzi - 2.d0*rik*drikdzi

               dfdxj=2.d0*rjk*drjkdxj - 2.d0*rij*drijdxj
               dfdyj=2.d0*rjk*drjkdyj - 2.d0*rij*drijdyj
               dfdzj=2.d0*rjk*drjkdzj - 2.d0*rij*drijdzj

               dfdxk=2.d0*rjk*drjkdxk - 2.d0*rik*drikdxk
               dfdyk=2.d0*rjk*drjkdyk - 2.d0*rik*drikdyk
               dfdzk=2.d0*rjk*drjkdzk - 2.d0*rik*drikdzk

               dgdxi=-2.d0*(drijdxi*rik + rij*drikdxi)
               dgdyi=-2.d0*(drijdyi*rik + rij*drikdyi)
               dgdzi=-2.d0*(drijdzi*rik + rij*drikdzi)

               dgdxj=-2.d0*drijdxj*rik
               dgdyj=-2.d0*drijdyj*rik
               dgdzj=-2.d0*drijdzj*rik

               dgdxk=-2.d0*rij*drikdxk
               dgdyk=-2.d0*rij*drikdyk
               dgdzk=-2.d0*rij*drikdzk

      dcosthetadxi=-(1.d0/2048.d0)*(1.d0-f/g)**15*(dfdxi*g-f*dgdxi)/g**2
      dcosthetadyi=-(1.d0/2048.d0)*(1.d0-f/g)**15*(dfdyi*g-f*dgdyi)/g**2
      dcosthetadzi=-(1.d0/2048.d0)*(1.d0-f/g)**15*(dfdzi*g-f*dgdzi)/g**2
      dcosthetadxj=-(1.d0/2048.d0)*(1.d0-f/g)**15*(dfdxj*g-f*dgdxj)/g**2
      dcosthetadyj=-(1.d0/2048.d0)*(1.d0-f/g)**15*(dfdyj*g-f*dgdyj)/g**2
      dcosthetadzj=-(1.d0/2048.d0)*(1.d0-f/g)**15*(dfdzj*g-f*dgdzj)/g**2
      dcosthetadxk=-(1.d0/2048.d0)*(1.d0-f/g)**15*(dfdxk*g-f*dgdxk)/g**2
      dcosthetadyk=-(1.d0/2048.d0)*(1.d0-f/g)**15*(dfdyk*g-f*dgdyk)/g**2
      dcosthetadzk=-(1.d0/2048.d0)*(1.d0-f/g)**15*(dfdzk*g-f*dgdzk)/g**2
           endif

! calculation of exponential term
      expxyz=dexp(-alpha*(rij**2+rik**2+rjk**2))
      temp1=-alpha*2.0d0*expxyz
      dexpxyzdxi=(rij*drijdxi+rik*drikdxi+rjk*drjkdxi)*temp1
      dexpxyzdyi=(rij*drijdyi+rik*drikdyi+rjk*drjkdyi)*temp1
      dexpxyzdzi=(rij*drijdzi+rik*drikdzi+rjk*drjkdzi)*temp1
      dexpxyzdxj=(rij*drijdxj+rik*drikdxj+rjk*drjkdxj)*temp1
      dexpxyzdyj=(rij*drijdyj+rik*drikdyj+rjk*drjkdyj)*temp1
      dexpxyzdzj=(rij*drijdzj+rik*drikdzj+rjk*drjkdzj)*temp1
      dexpxyzdxk=(rij*drijdxk+rik*drikdxk+rjk*drjkdxk)*temp1
      dexpxyzdyk=(rij*drijdyk+rik*drikdyk+rjk*drjkdyk)*temp1
      dexpxyzdzk=(rij*drijdzk+rik*drikdzk+rjk*drjkdzk)*temp1

      xx(ii,1,i)=xx(ii,1,i)+costheta*expxyz*fcutij*fcutik*fcutjk

! (fgh)' = f'gh + fg'h + fgh'
! (fghi)'= f'ghi + fg'hi + fgh'i + fghi'

      if(lforce)then
! for the derivatives fcutjk is a constant
        temp1=(+dcosthetadxi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxi)
        temp2=(+dcosthetadxj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxj)
        temp3=(+dcosthetadxk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxk)
        dxdy(ii,i,i,1)=dxdy(ii,i,i,1)+temp1
        dxdy(ii,i,n,1)=dxdy(ii,i,n,1)+temp2
        dxdy(ii,i,m,1)=dxdy(ii,i,m,1)+temp3
        strs(1,1,ii,i)=strs(1,1,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,1,ii,i)=strs(2,1,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,1,ii,i)=strs(3,1,ii,i)+deltazj*temp2&
              +deltazk*temp3

        temp1=(+dcosthetadyi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyi)
        temp2=(+dcosthetadyj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyj)
        temp3=(+dcosthetadyk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyk)
        dxdy(ii,i,i,2)=dxdy(ii,i,i,2)+temp1
        dxdy(ii,i,n,2)=dxdy(ii,i,n,2)+temp2
        dxdy(ii,i,m,2)=dxdy(ii,i,m,2)+temp3
        strs(1,2,ii,i)=strs(1,2,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,2,ii,i)=strs(2,2,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,2,ii,i)=strs(3,2,ii,i)+deltazj*temp2&
              +deltazk*temp3

        temp1=(+dcosthetadzi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzi)
        temp2=(+dcosthetadzj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzj)
        temp3=(+dcosthetadzk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzk)
        dxdy(ii,i,i,3)=dxdy(ii,i,i,3)+temp1
        dxdy(ii,i,n,3)=dxdy(ii,i,n,3)+temp2
        dxdy(ii,i,m,3)=dxdy(ii,i,m,3)+temp3
        strs(1,3,ii,i)=strs(1,3,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,3,ii,i)=strs(2,3,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,3,ii,i)=strs(3,3,ii,i)+deltazj*temp2&
              +deltazk*temp3
      endif

            endif ! rjk .le. cutoff
          endif ! rik .le. cutoff
          endif ! j .ne. k
        enddo ! k
        endif ! rij .le. cutoff
      enddo ! j
      enddo ! i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dangling bond vector
! cutoff is used from neural.in, alpha is not used
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      elseif(ntype.eq.11) then ! dangling bond vector
      do i=1,nat
        if(ltest) write(*,*)'calculating dangling bond for atom ',i
        dangvec(1)=0.0d0
        dangvec(2)=0.0d0
        dangvec(3)=0.0d0

      do j=lsta(1,i),lsta(2,i) ! loop over all neighbors of atom i
        n=lstc(j)
        rij=lstb(j,4)
        if(rij.le.cutoff) then ! rij < cutoff

        deltaxj=xyz(i,1)-lstb(j,1)
        deltayj=xyz(i,2)-lstb(j,2)
        deltazj=xyz(i,3)-lstb(j,3)

        fcutij=0.5d0*(dcos(pi*rij/cutoff)+1.d0)

        dangvec(1)=dangvec(1)+deltaxj*fcutij
        dangvec(2)=dangvec(2)+deltayj*fcutij
        dangvec(3)=dangvec(3)+deltazj*fcutij
        endif ! rij .le. cutoff
      enddo ! j

! calculate the symmetry function (length of the dangling bond vector):
      dangvector=0.0d0
      dangvector=dangvec(1)**2 + dangvec(2)**2 + dangvec(3)**2
      xx(ii,1,i)=dangvector

! calculate forces
      if(lforce) then
      do j=lsta(1,i),lsta(2,i) ! loop over all neighbors of atom i
        n=lstc(j)
        rij=lstb(j,4)

        if(rij.le.cutoff) then ! rij < cutoff

        deltaxj=(xyz(i,1)-lstb(j,1))
        deltayj=(xyz(i,2)-lstb(j,2))
        deltazj=(xyz(i,3)-lstb(j,3))
        drijdxi=deltaxj/rij
        drijdyi=deltayj/rij
        drijdzi=deltazj/rij
        drijdxj=-1.d0*drijdxi
        drijdyj=-1.d0*drijdyi
        drijdzj=-1.d0*drijdzi

        fcutij=0.5d0*(dcos(pi*rij/cutoff)+1.d0)
        temp1=0.5d0*(-dsin(pi*rij/cutoff))*(pi/cutoff)
        dfcutijdxi=temp1*drijdxi
        dfcutijdyi=temp1*drijdyi
        dfcutijdzi=temp1*drijdzi
        dfcutijdxj=-1.d0*dfcutijdxi
        dfcutijdyj=-1.d0*dfcutijdyi
        dfcutijdzj=-1.d0*dfcutijdzi

! d/dx_i !ok
        dxdy(ii,i,i,1)=dxdy(ii,i,i,1)&
        + 2.d0*dangvec(1)*fcutij&
        + 2.d0*dangvec(1)*deltaxj*dfcutijdxi&
        + 2.d0*dangvec(2)*deltayj*dfcutijdxi&
        + 2.d0*dangvec(3)*deltazj*dfcutijdxi

! d/dx_j
        dxdy(ii,i,n,1)=dxdy(ii,i,n,1)&
        - 2.d0*dangvec(1)*fcutij&
        + 2.d0*dangvec(1)*deltaxj*dfcutijdxj&
        + 2.d0*dangvec(2)*deltayj*dfcutijdxj&
        + 2.d0*dangvec(3)*deltazj*dfcutijdxj

! d/dy_i
        dxdy(ii,i,i,2)=dxdy(ii,i,i,2)&
        + 2.d0*dangvec(1)*deltaxj*dfcutijdyi&
        + 2.d0*dangvec(2)*fcutij&
        + 2.d0*dangvec(2)*deltayj*dfcutijdyi&
        + 2.d0*dangvec(3)*deltazj*dfcutijdyi

! d/dy_j
        dxdy(ii,i,n,2)=dxdy(ii,i,n,2)&
        + 2.d0*dangvec(1)*deltaxj*dfcutijdyj&
        - 2.d0*dangvec(2)*fcutij&
        + 2.d0*dangvec(2)*deltayj*dfcutijdyj&
        + 2.d0*dangvec(3)*deltazj*dfcutijdyj

! d/dz_i
        dxdy(ii,i,i,3)=dxdy(ii,i,i,3)&
        + 2.d0*dangvec(1)*deltaxj*dfcutijdzi&
        + 2.d0*dangvec(2)*deltayj*dfcutijdzi&
        + 2.d0*dangvec(3)*fcutij&
        + 2.d0*dangvec(3)*deltazj*dfcutijdzi

! d/dz_j
        !dxdy(ii,i,n,3)=dxdy(ii,i,n,3)+&
  ! --------------------fyuhaoy--------------------
	dxdy(ii,i,n,3)=dxdy(ii,i,n,3)&
  ! --------------------fyuhaoy--------------------
        + 2.d0*dangvec(1)*deltaxj*dfcutijdzj&
        + 2.d0*dangvec(2)*deltayj*dfcutijdzj&
        - 2.d0*dangvec(3)*fcutij&
        + 2.d0*dangvec(3)*deltazj*dfcutijdzj

!! stress
!! dxx/dx
        strs(1,1,ii,i)=strs(1,1,ii,i)&
        +deltaxj*2.0d0*(&
        - dangvec(1)*fcutij&
        + dangvec(1)*deltaxj*dfcutijdxj&
        + dangvec(2)*deltayj*dfcutijdxj&
        + dangvec(3)*deltazj*dfcutijdxj)

        strs(2,1,ii,i)=strs(2,1,ii,i)&
        +deltayj*2.0d0*(&
        - dangvec(1)*fcutij&
        + dangvec(1)*deltaxj*dfcutijdxj&
        + dangvec(2)*deltayj*dfcutijdxj&
        + dangvec(3)*deltazj*dfcutijdxj)

        strs(3,1,ii,i)=strs(3,1,ii,i)&
        +deltazj*2.0d0*(&
        - dangvec(1)*fcutij&
        + dangvec(1)*deltaxj*dfcutijdxj&
        + dangvec(2)*deltayj*dfcutijdxj&
        + dangvec(3)*deltazj*dfcutijdxj)
!! dxx/dy
        strs(1,2,ii,i)=strs(1,2,ii,i)&
        +deltaxj*2.0d0*(&
        + dangvec(1)*deltaxj*dfcutijdyj&
        - dangvec(2)*fcutij&
        + dangvec(2)*deltayj*dfcutijdyj&
        + dangvec(3)*deltazj*dfcutijdyj)

        strs(2,2,ii,i)=strs(2,2,ii,i)&
        +deltayj*2.0d0*(&
        + dangvec(1)*deltaxj*dfcutijdyj&
        - dangvec(2)*fcutij&
        + dangvec(2)*deltayj*dfcutijdyj&
        + dangvec(3)*deltazj*dfcutijdyj)

        strs(3,2,ii,i)=strs(3,2,ii,i)&
        +deltazj*2.0d0*(&
        + dangvec(1)*deltaxj*dfcutijdyj&
        - dangvec(2)*fcutij&
        + dangvec(2)*deltayj*dfcutijdyj&
        + dangvec(3)*deltazj*dfcutijdyj)

!! dxx/dz
        strs(1,3,ii,i)=strs(1,3,ii,i)&
        +deltaxj*2.0d0*(&
        + dangvec(1)*deltaxj*dfcutijdzj&
        + dangvec(2)*deltayj*dfcutijdzj&
        - dangvec(3)*fcutij&
        + dangvec(3)*deltazj*dfcutijdzj)

        strs(2,3,ii,i)=strs(2,3,ii,i)&
        +deltayj*2.0d0*(&
        + dangvec(1)*deltaxj*dfcutijdzj&
        + dangvec(2)*deltayj*dfcutijdzj&
        - dangvec(3)*fcutij&
        + dangvec(3)*deltazj*dfcutijdzj)

        strs(3,3,ii,i)=strs(3,3,ii,i)&
        +deltazj*2.0d0*(&
        + dangvec(1)*deltaxj*dfcutijdzj&
        + dangvec(2)*deltayj*dfcutijdzj&
        - dangvec(3)*fcutij&
        + dangvec(3)*deltazj*dfcutijdzj)

        endif ! rij .le. cutoff
      enddo ! j
      endif ! lforce

      enddo ! i

! end dangling bond vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      elseif(ntype.eq.12) then ! three-body term with global cutoff
      do i=1,nat
        if(ltest) write(*,*)'calculating three body term for atom ',i

      do j=lsta(1,i),lsta(2,i) ! loop over all neighbors of this atom
        n=lstc(j)
        rij=lstb(j,4)

        if(rij.le.cutoff)then

        deltaxj=-1.d0*(xyz(i,1)-lstb(j,1))
        deltayj=-1.d0*(xyz(i,2)-lstb(j,2))
        deltazj=-1.d0*(xyz(i,3)-lstb(j,3))
        drijdxi=-deltaxj/rij
        drijdyi=-deltayj/rij
        drijdzi=-deltazj/rij
        drijdxj=-1.d0*drijdxi
        drijdyj=-1.d0*drijdyi
        drijdzj=-1.d0*drijdzi
        drijdxk=0.0d0
        drijdyk=0.0d0
        drijdzk=0.0d0

        fcutij=0.5d0*(dcos(pi*rij/cutoff)+1.d0)
        temp1=0.5d0*(-dsin(pi*rij/cutoff))*(pi/cutoff)
        dfcutijdxi=temp1*drijdxi
        dfcutijdyi=temp1*drijdyi
        dfcutijdzi=temp1*drijdzi
        dfcutijdxj=-1.d0*dfcutijdxi
        dfcutijdyj=-1.d0*dfcutijdyi
        dfcutijdzj=-1.d0*dfcutijdzi
        dfcutijdxk=0.0d0
        dfcutijdyk=0.0d0
        dfcutijdzk=0.0d0

        do k=lsta(1,i),lsta(2,i) ! loop over all neighbors of this atom
          m=lstc(k)

!          if(j.ne.k) then ! neighbors j and k must be different
          if(k.gt.j) then ! neighbors j and k must be different

          rik=lstb(k,4)

          if(rik.le.cutoff)then

          deltaxk=-1.d0*(xyz(i,1)-lstb(k,1))
          deltayk=-1.d0*(xyz(i,2)-lstb(k,2))
          deltazk=-1.d0*(xyz(i,3)-lstb(k,3))
          drikdxi=-deltaxk/rik
          drikdyi=-deltayk/rik
          drikdzi=-deltazk/rik
          drikdxk=-1.d0*drikdxi
          drikdyk=-1.d0*drikdyi
          drikdzk=-1.d0*drikdzi
          drikdxj=0.0d0
          drikdyj=0.0d0
          drikdzj=0.0d0

          fcutik=0.5d0*(dcos(pi*rik/cutoff)+1.d0)
          temp1=0.5d0*(-dsin(pi*rik/cutoff))*(pi/cutoff)
          dfcutikdxi=temp1*drikdxi
          dfcutikdyi=temp1*drikdyi
          dfcutikdzi=temp1*drikdzi
          dfcutikdxj=0.0d0
          dfcutikdyj=0.0d0
          dfcutikdzj=0.0d0
          dfcutikdxk=-1.d0*dfcutikdxi
          dfcutikdyk=-1.d0*dfcutikdyi
          dfcutikdzk=-1.d0*dfcutikdzi

          rjk=(lstb(j,1)-lstb(k,1))**2 +&
              (lstb(j,2)-lstb(k,2))**2 +&
              (lstb(j,3)-lstb(k,3))**2
          rjk=dsqrt(rjk)

          if(rjk.le.0.5d0) then
            write(*,*)'Error rjk .le. 0.5',ntype
!            stop
          endif

          if(rjk.le.cutoff) then
          drjkdxj=(lstb(j,1)-lstb(k,1))/rjk
          drjkdyj=(lstb(j,2)-lstb(k,2))/rjk
          drjkdzj=(lstb(j,3)-lstb(k,3))/rjk
          drjkdxk=-1.d0*drjkdxj
          drjkdyk=-1.d0*drjkdyj
          drjkdzk=-1.d0*drjkdzj
          drjkdxi=0.0d0
          drjkdyi=0.0d0
          drjkdzi=0.0d0

          fcutjk=0.5d0*(dcos(pi*rjk/cutoff)+1.d0)
          temp1=0.5d0*(-dsin(pi*rjk/cutoff))*(pi/cutoff)
          dfcutjkdxj=temp1*drjkdxj
          dfcutjkdyj=temp1*drjkdyj
          dfcutjkdzj=temp1*drjkdzj
          dfcutjkdxk=-1.d0*dfcutjkdxj
          dfcutjkdyk=-1.d0*dfcutjkdyj
          dfcutjkdzk=-1.d0*dfcutjkdzj
          dfcutjkdxi=0.0d0
          dfcutjkdyi=0.0d0
          dfcutjkdzi=0.0d0

! costheta=(rjk**2 - rij**2 -rik**2)/(-2.d0*rij*rik)
            f=rjk**2 - rij**2 -rik**2
            g=-2.d0*rij*rik
            if(lforce) then
               dfdxi=-2.d0*rij*drijdxi - 2.d0*rik*drikdxi
               dfdyi=-2.d0*rij*drijdyi - 2.d0*rik*drikdyi
               dfdzi=-2.d0*rij*drijdzi - 2.d0*rik*drikdzi

               dfdxj=2.d0*rjk*drjkdxj - 2.d0*rij*drijdxj
               dfdyj=2.d0*rjk*drjkdyj - 2.d0*rij*drijdyj
               dfdzj=2.d0*rjk*drjkdzj - 2.d0*rij*drijdzj

               dfdxk=2.d0*rjk*drjkdxk - 2.d0*rik*drikdxk
               dfdyk=2.d0*rjk*drjkdyk - 2.d0*rik*drikdyk
               dfdzk=2.d0*rjk*drjkdzk - 2.d0*rik*drikdzk

               dgdxi=-2.d0*(drijdxi*rik + rij*drikdxi)
               dgdyi=-2.d0*(drijdyi*rik + rij*drikdyi)
               dgdzi=-2.d0*(drijdzi*rik + rij*drikdzi)

               dgdxj=-2.d0*drijdxj*rik
               dgdyj=-2.d0*drijdyj*rik
               dgdzj=-2.d0*drijdzj*rik

               dgdxk=-2.d0*rij*drikdxk
               dgdyk=-2.d0*rij*drikdyk
               dgdzk=-2.d0*rij*drikdzk
            endif

            costheta=f/g
            if(costheta.lt.-1.0d0)then
              costheta=costheta+0.000000000001d0
            endif
            if(costheta.lt.-1.0d0)then
              write(*,'(a,f20.16)')'ERROR costheta ',costheta
              stop
            endif
            if(costheta.gt.1.0d0)then
              costheta=costheta-0.000000000001d0
            endif
            if(costheta.gt.1.0d0)then
              write(*,'(a,f20.16)')'ERROR costheta ',costheta
              stop
            endif

            if(lforce) then
               dcosthetadxi=(dfdxi*g - f*dgdxi)/g**2
               dcosthetadyi=(dfdyi*g - f*dgdyi)/g**2
               dcosthetadzi=(dfdzi*g - f*dgdzi)/g**2
               dcosthetadxj=(dfdxj*g - f*dgdxj)/g**2
               dcosthetadyj=(dfdyj*g - f*dgdyj)/g**2
               dcosthetadzj=(dfdzj*g - f*dgdzj)/g**2
               dcosthetadxk=(dfdxk*g - f*dgdxk)/g**2
               dcosthetadyk=(dfdyk*g - f*dgdyk)/g**2
               dcosthetadzk=(dfdzk*g - f*dgdzk)/g**2
            endif

            theta=acos(costheta)
            theta=2.d0*theta

            if(lforce) then
              temp1=-2.d0/sqrt(1-costheta**2)
              dthetadxi=temp1*dcosthetadxi
              dthetadyi=temp1*dcosthetadyi
              dthetadzi=temp1*dcosthetadzi
              dthetadxj=temp1*dcosthetadxj
              dthetadyj=temp1*dcosthetadyj
              dthetadzj=temp1*dcosthetadzj
              dthetadxk=temp1*dcosthetadxk
              dthetadyk=temp1*dcosthetadyk
              dthetadzk=temp1*dcosthetadzk
            endif

            costheta=cos(theta)
            costheta=costheta+1.d0 ! avoid negative values
            if(lforce) then
               temp1=-1.d0*sin(theta)
               dcosthetadxi=temp1*dthetadxi
               dcosthetadyi=temp1*dthetadyi
               dcosthetadzi=temp1*dthetadzi
               dcosthetadxj=temp1*dthetadxj
               dcosthetadyj=temp1*dthetadyj
               dcosthetadzj=temp1*dthetadzj
               dcosthetadxk=temp1*dthetadxk
               dcosthetadyk=temp1*dthetadyk
               dcosthetadzk=temp1*dthetadzk
            endif

! calculation of exponential term
      expxyz=dexp(-alpha*(rij**2+rik**2+rjk**2))
      temp1=-alpha*2.0d0*expxyz
      dexpxyzdxi=(rij*drijdxi+rik*drikdxi+rjk*drjkdxi)*temp1
      dexpxyzdyi=(rij*drijdyi+rik*drikdyi+rjk*drjkdyi)*temp1
      dexpxyzdzi=(rij*drijdzi+rik*drikdzi+rjk*drjkdzi)*temp1
      dexpxyzdxj=(rij*drijdxj+rik*drikdxj+rjk*drjkdxj)*temp1
      dexpxyzdyj=(rij*drijdyj+rik*drikdyj+rjk*drjkdyj)*temp1
      dexpxyzdzj=(rij*drijdzj+rik*drikdzj+rjk*drjkdzj)*temp1
      dexpxyzdxk=(rij*drijdxk+rik*drikdxk+rjk*drjkdxk)*temp1
      dexpxyzdyk=(rij*drijdyk+rik*drikdyk+rjk*drjkdyk)*temp1
      dexpxyzdzk=(rij*drijdzk+rik*drikdzk+rjk*drjkdzk)*temp1

      xx(ii,1,i)=xx(ii,1,i)+costheta*expxyz*fcutij*fcutik*fcutjk

      if(lforce) then
! dxxii/dx_i
        temp1=(+dcosthetadxi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxi)
        temp2=(+dcosthetadxj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxj)
        temp3=(+dcosthetadxk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxk)
        dxdy(ii,i,i,1)=dxdy(ii,i,i,1)+temp1
        dxdy(ii,i,n,1)=dxdy(ii,i,n,1)+temp2
        dxdy(ii,i,m,1)=dxdy(ii,i,m,1)+temp3
        strs(1,1,ii,i)=strs(1,1,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,1,ii,i)=strs(2,1,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,1,ii,i)=strs(3,1,ii,i)+deltazj*temp2&
              +deltazk*temp3
! dxxii/dy_i
        temp1=(+dcosthetadyi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyi)
        temp2=(+dcosthetadyj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyj)
        temp3=(+dcosthetadyk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyk)
        dxdy(ii,i,i,2)=dxdy(ii,i,i,2)+temp1
        dxdy(ii,i,n,2)=dxdy(ii,i,n,2)+temp2
        dxdy(ii,i,m,2)=dxdy(ii,i,m,2)+temp3
        strs(1,2,ii,i)=strs(1,2,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,2,ii,i)=strs(2,2,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,2,ii,i)=strs(3,2,ii,i)+deltazj*temp2&
              +deltazk*temp3
!! dxxii/dz_i
        temp1=(+dcosthetadzi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzi)
        temp2=(+dcosthetadzj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzj)
        temp3=(+dcosthetadzk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzk)
        dxdy(ii,i,i,3)=dxdy(ii,i,i,3)+temp1
        dxdy(ii,i,n,3)=dxdy(ii,i,n,3)+temp2
        dxdy(ii,i,m,3)=dxdy(ii,i,m,3)+temp3
        strs(1,3,ii,i)=strs(1,3,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,3,ii,i)=strs(2,3,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,3,ii,i)=strs(3,3,ii,i)+deltazj*temp2&
              +deltazk*temp3
      endif

            endif ! rjk .le. cutoff
            endif ! rik .le. cutoff
          endif ! j .ne. k
        enddo ! k
        endif ! rij .le. cutoff
      enddo ! j
      enddo ! i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      elseif(ntype.eq.13) then ! three-body term with global cutoff
      do i=1,nat
        if(ltest) write(*,*)'calculating three body term for atom ',i

      do j=lsta(1,i),lsta(2,i) ! loop over all neighbors of this atom
        n=lstc(j)
        rij=lstb(j,4)

        if(rij.le.cutoff)then

        deltaxj=-1.d0*(xyz(i,1)-lstb(j,1))
        deltayj=-1.d0*(xyz(i,2)-lstb(j,2))
        deltazj=-1.d0*(xyz(i,3)-lstb(j,3))
        drijdxi=-deltaxj/rij
        drijdyi=-deltayj/rij
        drijdzi=-deltazj/rij
        drijdxj=-1.d0*drijdxi
        drijdyj=-1.d0*drijdyi
        drijdzj=-1.d0*drijdzi
        drijdxk=0.0d0
        drijdyk=0.0d0
        drijdzk=0.0d0

        fcutij=0.5d0*(dcos(pi*rij/cutoff)+1.d0)
        temp1=0.5d0*(-dsin(pi*rij/cutoff))*(pi/cutoff)
        dfcutijdxi=temp1*drijdxi
        dfcutijdyi=temp1*drijdyi
        dfcutijdzi=temp1*drijdzi
        dfcutijdxj=-1.d0*dfcutijdxi
        dfcutijdyj=-1.d0*dfcutijdyi
        dfcutijdzj=-1.d0*dfcutijdzi
        dfcutijdxk=0.0d0
        dfcutijdyk=0.0d0
        dfcutijdzk=0.0d0

        do k=lsta(1,i),lsta(2,i) ! loop over all neighbors of this atom
          m=lstc(k)

!          if(j.ne.k) then ! neighbors j and k must be different
          if(k.gt.j) then ! neighbors j and k must be different

          rik=lstb(k,4)

          if(rik.le.cutoff)then

          deltaxk=-1.d0*(xyz(i,1)-lstb(k,1))
          deltayk=-1.d0*(xyz(i,2)-lstb(k,2))
          deltazk=-1.d0*(xyz(i,3)-lstb(k,3))
          drikdxi=-deltaxk/rik
          drikdyi=-deltayk/rik
          drikdzi=-deltazk/rik
          drikdxk=-1.d0*drikdxi
          drikdyk=-1.d0*drikdyi
          drikdzk=-1.d0*drikdzi
          drikdxj=0.0d0
          drikdyj=0.0d0
          drikdzj=0.0d0

          fcutik=0.5d0*(dcos(pi*rik/cutoff)+1.d0)
          temp1=0.5d0*(-dsin(pi*rik/cutoff))*(pi/cutoff)
          dfcutikdxi=temp1*drikdxi
          dfcutikdyi=temp1*drikdyi
          dfcutikdzi=temp1*drikdzi
          dfcutikdxj=0.0d0
          dfcutikdyj=0.0d0
          dfcutikdzj=0.0d0
          dfcutikdxk=-1.d0*dfcutikdxi
          dfcutikdyk=-1.d0*dfcutikdyi
          dfcutikdzk=-1.d0*dfcutikdzi

          rjk=(lstb(j,1)-lstb(k,1))**2 +&
              (lstb(j,2)-lstb(k,2))**2 +&
              (lstb(j,3)-lstb(k,3))**2
          rjk=dsqrt(rjk)

          if(rjk.le.0.5d0) then
            write(*,*)'Error rjk .le. 0.5',ntype
!            stop
          endif

          if(rjk.le.cutoff) then
          drjkdxj=(lstb(j,1)-lstb(k,1))/rjk
          drjkdyj=(lstb(j,2)-lstb(k,2))/rjk
          drjkdzj=(lstb(j,3)-lstb(k,3))/rjk
          drjkdxk=-1.d0*drjkdxj
          drjkdyk=-1.d0*drjkdyj
          drjkdzk=-1.d0*drjkdzj
          drjkdxi=0.0d0
          drjkdyi=0.0d0
          drjkdzi=0.0d0

          fcutjk=0.5d0*(dcos(pi*rjk/cutoff)+1.d0)
          temp1=0.5d0*(-dsin(pi*rjk/cutoff))*(pi/cutoff)
          dfcutjkdxj=temp1*drjkdxj
          dfcutjkdyj=temp1*drjkdyj
          dfcutjkdzj=temp1*drjkdzj
          dfcutjkdxk=-1.d0*dfcutjkdxj
          dfcutjkdyk=-1.d0*dfcutjkdyj
          dfcutjkdzk=-1.d0*dfcutjkdzj
          dfcutjkdxi=0.0d0
          dfcutjkdyi=0.0d0
          dfcutjkdzi=0.0d0

! costheta=(rjk**2 - rij**2 -rik**2)/(-2.d0*rij*rik)
            f=rjk**2 - rij**2 -rik**2
            g=-2.d0*rij*rik
            if(lforce) then
               dfdxi=-2.d0*rij*drijdxi - 2.d0*rik*drikdxi
               dfdyi=-2.d0*rij*drijdyi - 2.d0*rik*drikdyi
               dfdzi=-2.d0*rij*drijdzi - 2.d0*rik*drikdzi

               dfdxj=2.d0*rjk*drjkdxj - 2.d0*rij*drijdxj
               dfdyj=2.d0*rjk*drjkdyj - 2.d0*rij*drijdyj
               dfdzj=2.d0*rjk*drjkdzj - 2.d0*rij*drijdzj

               dfdxk=2.d0*rjk*drjkdxk - 2.d0*rik*drikdxk
               dfdyk=2.d0*rjk*drjkdyk - 2.d0*rik*drikdyk
               dfdzk=2.d0*rjk*drjkdzk - 2.d0*rik*drikdzk

               dgdxi=-2.d0*(drijdxi*rik + rij*drikdxi)
               dgdyi=-2.d0*(drijdyi*rik + rij*drikdyi)
               dgdzi=-2.d0*(drijdzi*rik + rij*drikdzi)

               dgdxj=-2.d0*drijdxj*rik
               dgdyj=-2.d0*drijdyj*rik
               dgdzj=-2.d0*drijdzj*rik

               dgdxk=-2.d0*rij*drikdxk
               dgdyk=-2.d0*rij*drikdyk
               dgdzk=-2.d0*rij*drikdzk
            endif

            costheta=f/g
            if(costheta.lt.-1.0d0)then
              costheta=costheta+0.000000000001d0
            endif
            if(costheta.lt.-1.0d0)then
              write(*,'(a,f20.16)')'ERROR costheta ',costheta
              stop
            endif
            if(costheta.gt.1.0d0)then
              costheta=costheta-0.000000000001d0
            endif
            if(costheta.gt.1.0d0)then
              write(*,'(a,f20.16)')'ERROR costheta ',costheta
              stop
            endif

            if(lforce) then
               dcosthetadxi=(dfdxi*g - f*dgdxi)/g**2
               dcosthetadyi=(dfdyi*g - f*dgdyi)/g**2
               dcosthetadzi=(dfdzi*g - f*dgdzi)/g**2
               dcosthetadxj=(dfdxj*g - f*dgdxj)/g**2
               dcosthetadyj=(dfdyj*g - f*dgdyj)/g**2
               dcosthetadzj=(dfdzj*g - f*dgdzj)/g**2
               dcosthetadxk=(dfdxk*g - f*dgdxk)/g**2
               dcosthetadyk=(dfdyk*g - f*dgdyk)/g**2
               dcosthetadzk=(dfdzk*g - f*dgdzk)/g**2
            endif

            theta=acos(costheta)
            theta=4.d0*theta

            if(lforce) then
              temp1=-4.d0/sqrt(1-costheta**2)
              dthetadxi=temp1*dcosthetadxi
              dthetadyi=temp1*dcosthetadyi
              dthetadzi=temp1*dcosthetadzi
              dthetadxj=temp1*dcosthetadxj
              dthetadyj=temp1*dcosthetadyj
              dthetadzj=temp1*dcosthetadzj
              dthetadxk=temp1*dcosthetadxk
              dthetadyk=temp1*dcosthetadyk
              dthetadzk=temp1*dcosthetadzk
            endif

            costheta=cos(theta)
            costheta=costheta+1.d0 ! avoid negative values
            if(lforce) then
               temp1=-1.d0*sin(theta)
               dcosthetadxi=temp1*dthetadxi
               dcosthetadyi=temp1*dthetadyi
               dcosthetadzi=temp1*dthetadzi
               dcosthetadxj=temp1*dthetadxj
               dcosthetadyj=temp1*dthetadyj
               dcosthetadzj=temp1*dthetadzj
               dcosthetadxk=temp1*dthetadxk
               dcosthetadyk=temp1*dthetadyk
               dcosthetadzk=temp1*dthetadzk
            endif

! calculation of exponential term
      expxyz=dexp(-alpha*(rij**2+rik**2+rjk**2))
      temp1=-alpha*2.0d0*expxyz
      dexpxyzdxi=(rij*drijdxi+rik*drikdxi+rjk*drjkdxi)*temp1
      dexpxyzdyi=(rij*drijdyi+rik*drikdyi+rjk*drjkdyi)*temp1
      dexpxyzdzi=(rij*drijdzi+rik*drikdzi+rjk*drjkdzi)*temp1
      dexpxyzdxj=(rij*drijdxj+rik*drikdxj+rjk*drjkdxj)*temp1
      dexpxyzdyj=(rij*drijdyj+rik*drikdyj+rjk*drjkdyj)*temp1
      dexpxyzdzj=(rij*drijdzj+rik*drikdzj+rjk*drjkdzj)*temp1
      dexpxyzdxk=(rij*drijdxk+rik*drikdxk+rjk*drjkdxk)*temp1
      dexpxyzdyk=(rij*drijdyk+rik*drikdyk+rjk*drjkdyk)*temp1
      dexpxyzdzk=(rij*drijdzk+rik*drikdzk+rjk*drjkdzk)*temp1

      xx(ii,1,i)=xx(ii,1,i)+costheta*expxyz*fcutij*fcutik*fcutjk

      if(lforce) then
! dxxii/dx_i
        temp1=(+dcosthetadxi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxi)
        temp2=(+dcosthetadxj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxj)
        temp3=(+dcosthetadxk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxk)
        dxdy(ii,i,i,1)=dxdy(ii,i,i,1)+temp1
        dxdy(ii,i,n,1)=dxdy(ii,i,n,1)+temp2
        dxdy(ii,i,m,1)=dxdy(ii,i,m,1)+temp3
        strs(1,1,ii,i)=strs(1,1,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,1,ii,i)=strs(2,1,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,1,ii,i)=strs(3,1,ii,i)+deltazj*temp2&
              +deltazk*temp3
! dxxii/dy_i
        temp1=(+dcosthetadyi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyi)
        temp2=(+dcosthetadyj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyj)
        temp3=(+dcosthetadyk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyk)
        dxdy(ii,i,i,2)=dxdy(ii,i,i,2)+temp1
        dxdy(ii,i,n,2)=dxdy(ii,i,n,2)+temp2
        dxdy(ii,i,m,2)=dxdy(ii,i,m,2)+temp3
        strs(1,2,ii,i)=strs(1,2,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,2,ii,i)=strs(2,2,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,2,ii,i)=strs(3,2,ii,i)+deltazj*temp2&
              +deltazk*temp3
!! dxxii/dz_i
        temp1=(+dcosthetadzi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzi)
        temp2=(+dcosthetadzj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzj)
        temp3=(+dcosthetadzk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzk)
        dxdy(ii,i,i,3)=dxdy(ii,i,i,3)+temp1
        dxdy(ii,i,n,3)=dxdy(ii,i,n,3)+temp2
        dxdy(ii,i,m,3)=dxdy(ii,i,m,3)+temp3
        strs(1,3,ii,i)=strs(1,3,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,3,ii,i)=strs(2,3,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,3,ii,i)=strs(3,3,ii,i)+deltazj*temp2&
              +deltazk*temp3
      endif

            endif ! rjk .le. cutoff
            endif ! rik .le. cutoff
          endif ! j .ne. k
        enddo ! k
        endif ! rij .le. cutoff
      enddo ! j
      enddo ! i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      elseif(ntype.eq.14) then ! three-body term with global cutoff
      do i=1,nat
        if(ltest) write(*,*)'calculating three body term for atom ',i

      do j=lsta(1,i),lsta(2,i) ! loop over all neighbors of this atom
        n=lstc(j)
        rij=lstb(j,4)

        if(rij.le.cutoff)then

        deltaxj=-1.d0*(xyz(i,1)-lstb(j,1))
        deltayj=-1.d0*(xyz(i,2)-lstb(j,2))
        deltazj=-1.d0*(xyz(i,3)-lstb(j,3))
        drijdxi=-deltaxj/rij
        drijdyi=-deltayj/rij
        drijdzi=-deltazj/rij
        drijdxj=-1.d0*drijdxi
        drijdyj=-1.d0*drijdyi
        drijdzj=-1.d0*drijdzi
        drijdxk=0.0d0
        drijdyk=0.0d0
        drijdzk=0.0d0

        fcutij=0.5d0*(dcos(pi*rij/cutoff)+1.d0)
        temp1=0.5d0*(-dsin(pi*rij/cutoff))*(pi/cutoff)
        dfcutijdxi=temp1*drijdxi
        dfcutijdyi=temp1*drijdyi
        dfcutijdzi=temp1*drijdzi
        dfcutijdxj=-1.d0*dfcutijdxi
        dfcutijdyj=-1.d0*dfcutijdyi
        dfcutijdzj=-1.d0*dfcutijdzi
        dfcutijdxk=0.0d0
        dfcutijdyk=0.0d0
        dfcutijdzk=0.0d0

        do k=lsta(1,i),lsta(2,i) ! loop over all neighbors of this atom
          m=lstc(k)

!          if(j.ne.k) then ! neighbors j and k must be different
          if(k.gt.j) then ! neighbors j and k must be different

          rik=lstb(k,4)

          if(rik.le.cutoff)then
          deltaxk=-1.d0*(xyz(i,1)-lstb(k,1))
          deltayk=-1.d0*(xyz(i,2)-lstb(k,2))
          deltazk=-1.d0*(xyz(i,3)-lstb(k,3))
          drikdxi=-deltaxk/rik
          drikdyi=-deltayk/rik
          drikdzi=-deltazk/rik
          drikdxk=-1.d0*drikdxi
          drikdyk=-1.d0*drikdyi
          drikdzk=-1.d0*drikdzi
          drikdxj=0.0d0
          drikdyj=0.0d0
          drikdzj=0.0d0

          fcutik=0.5d0*(dcos(pi*rik/cutoff)+1.d0)
          temp1=0.5d0*(-dsin(pi*rik/cutoff))*(pi/cutoff)
          dfcutikdxi=temp1*drikdxi
          dfcutikdyi=temp1*drikdyi
          dfcutikdzi=temp1*drikdzi
          dfcutikdxj=0.0d0
          dfcutikdyj=0.0d0
          dfcutikdzj=0.0d0
          dfcutikdxk=-1.d0*dfcutikdxi
          dfcutikdyk=-1.d0*dfcutikdyi
          dfcutikdzk=-1.d0*dfcutikdzi

          rjk=(lstb(j,1)-lstb(k,1))**2 +&
              (lstb(j,2)-lstb(k,2))**2 +&
              (lstb(j,3)-lstb(k,3))**2
          rjk=dsqrt(rjk)

          if(rjk.le.0.5d0) then
            write(*,*)'Error rjk .le. 0.5',ntype
!            stop
          endif

          if(rjk.le.cutoff) then
          drjkdxj=(lstb(j,1)-lstb(k,1))/rjk
          drjkdyj=(lstb(j,2)-lstb(k,2))/rjk
          drjkdzj=(lstb(j,3)-lstb(k,3))/rjk
          drjkdxk=-1.d0*drjkdxj
          drjkdyk=-1.d0*drjkdyj
          drjkdzk=-1.d0*drjkdzj
          drjkdxi=0.0d0
          drjkdyi=0.0d0
          drjkdzi=0.0d0

          fcutjk=0.5d0*(dcos(pi*rjk/cutoff)+1.d0)
          temp1=0.5d0*(-dsin(pi*rjk/cutoff))*(pi/cutoff)
          dfcutjkdxj=temp1*drjkdxj
          dfcutjkdyj=temp1*drjkdyj
          dfcutjkdzj=temp1*drjkdzj
          dfcutjkdxk=-1.d0*dfcutjkdxj
          dfcutjkdyk=-1.d0*dfcutjkdyj
          dfcutjkdzk=-1.d0*dfcutjkdzj
          dfcutjkdxi=0.0d0
          dfcutjkdyi=0.0d0
          dfcutjkdzi=0.0d0

! costheta=(rjk**2 - rij**2 -rik**2)/(-2.d0*rij*rik)
            f=rjk**2 - rij**2 -rik**2
            g=-2.d0*rij*rik
            if(lforce) then
               dfdxi=-2.d0*rij*drijdxi - 2.d0*rik*drikdxi
               dfdyi=-2.d0*rij*drijdyi - 2.d0*rik*drikdyi
               dfdzi=-2.d0*rij*drijdzi - 2.d0*rik*drikdzi

               dfdxj=2.d0*rjk*drjkdxj - 2.d0*rij*drijdxj
               dfdyj=2.d0*rjk*drjkdyj - 2.d0*rij*drijdyj
               dfdzj=2.d0*rjk*drjkdzj - 2.d0*rij*drijdzj

               dfdxk=2.d0*rjk*drjkdxk - 2.d0*rik*drikdxk
               dfdyk=2.d0*rjk*drjkdyk - 2.d0*rik*drikdyk
               dfdzk=2.d0*rjk*drjkdzk - 2.d0*rik*drikdzk

               dgdxi=-2.d0*(drijdxi*rik + rij*drikdxi)
               dgdyi=-2.d0*(drijdyi*rik + rij*drikdyi)
               dgdzi=-2.d0*(drijdzi*rik + rij*drikdzi)

               dgdxj=-2.d0*drijdxj*rik
               dgdyj=-2.d0*drijdyj*rik
               dgdzj=-2.d0*drijdzj*rik

               dgdxk=-2.d0*rij*drikdxk
               dgdyk=-2.d0*rij*drikdyk
               dgdzk=-2.d0*rij*drikdzk
            endif

            costheta=f/g
            if(costheta.lt.-1.0d0)then
              costheta=costheta+0.000000000001d0
            endif
            if(costheta.lt.-1.0d0)then
              write(*,'(a,f20.16)')'ERROR costheta ',costheta
              stop
            endif
            if(costheta.gt.1.0d0)then
              costheta=costheta-0.000000000001d0
            endif
            if(costheta.gt.1.0d0)then
              write(*,'(a,f20.16)')'ERROR costheta ',costheta
              stop
            endif

            if(lforce) then
               dcosthetadxi=(dfdxi*g - f*dgdxi)/g**2
               dcosthetadyi=(dfdyi*g - f*dgdyi)/g**2
               dcosthetadzi=(dfdzi*g - f*dgdzi)/g**2
               dcosthetadxj=(dfdxj*g - f*dgdxj)/g**2
               dcosthetadyj=(dfdyj*g - f*dgdyj)/g**2
               dcosthetadzj=(dfdzj*g - f*dgdzj)/g**2
               dcosthetadxk=(dfdxk*g - f*dgdxk)/g**2
               dcosthetadyk=(dfdyk*g - f*dgdyk)/g**2
               dcosthetadzk=(dfdzk*g - f*dgdzk)/g**2
            endif

            theta=acos(costheta)
            theta=8.d0*theta

            if(lforce) then
              temp1=-8.d0/sqrt(1-costheta**2)
              dthetadxi=temp1*dcosthetadxi
              dthetadyi=temp1*dcosthetadyi
              dthetadzi=temp1*dcosthetadzi
              dthetadxj=temp1*dcosthetadxj
              dthetadyj=temp1*dcosthetadyj
              dthetadzj=temp1*dcosthetadzj
              dthetadxk=temp1*dcosthetadxk
              dthetadyk=temp1*dcosthetadyk
              dthetadzk=temp1*dcosthetadzk
            endif

            costheta=cos(theta)
            costheta=costheta+1.d0 ! avoid negative values
            if(lforce) then
               temp1=-1.d0*sin(theta)
               dcosthetadxi=temp1*dthetadxi
               dcosthetadyi=temp1*dthetadyi
               dcosthetadzi=temp1*dthetadzi
               dcosthetadxj=temp1*dthetadxj
               dcosthetadyj=temp1*dthetadyj
               dcosthetadzj=temp1*dthetadzj
               dcosthetadxk=temp1*dthetadxk
               dcosthetadyk=temp1*dthetadyk
               dcosthetadzk=temp1*dthetadzk
            endif

! calculation of exponential term
      expxyz=dexp(-alpha*(rij**2+rik**2+rjk**2))
      temp1=-alpha*2.0d0*expxyz
      dexpxyzdxi=(rij*drijdxi+rik*drikdxi+rjk*drjkdxi)*temp1
      dexpxyzdyi=(rij*drijdyi+rik*drikdyi+rjk*drjkdyi)*temp1
      dexpxyzdzi=(rij*drijdzi+rik*drikdzi+rjk*drjkdzi)*temp1
      dexpxyzdxj=(rij*drijdxj+rik*drikdxj+rjk*drjkdxj)*temp1
      dexpxyzdyj=(rij*drijdyj+rik*drikdyj+rjk*drjkdyj)*temp1
      dexpxyzdzj=(rij*drijdzj+rik*drikdzj+rjk*drjkdzj)*temp1
      dexpxyzdxk=(rij*drijdxk+rik*drikdxk+rjk*drjkdxk)*temp1
      dexpxyzdyk=(rij*drijdyk+rik*drikdyk+rjk*drjkdyk)*temp1
      dexpxyzdzk=(rij*drijdzk+rik*drikdzk+rjk*drjkdzk)*temp1

      xx(ii,1,i)=xx(ii,1,i)+costheta*expxyz*fcutij*fcutik*fcutjk

      if(lforce) then
! dxxii/dx_i
        temp1=(+dcosthetadxi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxi)
        temp2=(+dcosthetadxj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxj)
        temp3=(+dcosthetadxk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdxk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdxk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdxk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdxk)
        dxdy(ii,i,i,1)=dxdy(ii,i,i,1)+temp1
        dxdy(ii,i,n,1)=dxdy(ii,i,n,1)+temp2
        dxdy(ii,i,m,1)=dxdy(ii,i,m,1)+temp3
        strs(1,1,ii,i)=strs(1,1,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,1,ii,i)=strs(2,1,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,1,ii,i)=strs(3,1,ii,i)+deltazj*temp2&
              +deltazk*temp3
! dxxii/dy_i
        temp1=(+dcosthetadyi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyi)
        temp2=(+dcosthetadyj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyj)
        temp3=(+dcosthetadyk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdyk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdyk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdyk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdyk)
        dxdy(ii,i,i,2)=dxdy(ii,i,i,2)+temp1
        dxdy(ii,i,n,2)=dxdy(ii,i,n,2)+temp2
        dxdy(ii,i,m,2)=dxdy(ii,i,m,2)+temp3
        strs(1,2,ii,i)=strs(1,2,ii,i)+deltaxj*temp2&
              +deltayk*temp3
        strs(3,2,ii,i)=strs(3,2,ii,i)+deltazj*temp2&
              +deltazk*temp3
!! dxxii/dz_i
        temp1=(+dcosthetadzi*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzi*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzi*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzi*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzi)
        temp2=(+dcosthetadzj*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzj*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzj*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzj*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzj)
        temp3=(+dcosthetadzk*expxyz*fcutij*fcutik*fcutjk&
              +costheta*dexpxyzdzk*fcutij*fcutik*fcutjk&
              +costheta*expxyz*dfcutijdzk*fcutik*fcutjk&
              +costheta*expxyz*fcutij*dfcutikdzk*fcutjk&
              +costheta*expxyz*fcutij*fcutik*dfcutjkdzk)
        dxdy(ii,i,i,3)=dxdy(ii,i,i,3)+temp1
        dxdy(ii,i,n,3)=dxdy(ii,i,n,3)+temp2
        dxdy(ii,i,m,3)=dxdy(ii,i,m,3)+temp3
        strs(1,3,ii,i)=strs(1,3,ii,i)+deltaxj*temp2&
              +deltaxk*temp3
        strs(2,3,ii,i)=strs(2,3,ii,i)+deltayj*temp2&
              +deltayk*temp3
        strs(3,3,ii,i)=strs(3,3,ii,i)+deltazj*temp2&
              +deltazk*temp3
      endif

            endif ! rjk .le. cutoff
            endif ! rik .le. cutoff
          endif ! j .ne. k
        enddo ! k
        endif ! rij .le. cutoff
      enddo ! j
      enddo ! i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      else
        write(*,*)'Unknown function type ',ntype,' in neural.in, ERROR'
        stop
      endif ! ntype

      enddo ! ii symmetry functions
      close(66)
!         do jj = 1,nn
!          do i = 1,nat
!             do j = 1,nat
!           if(j.ne.i)then
!      do ii = 1,nn
!             dxx(ii,i,j)=dxx(ii,i,j)+xx(ii,1,i)-xx(ii,1,j)
!             ddxxyy(ii,i,j,1)=dxdy(ii,i,i,1)-dxdy(ii,j,j,1)
!             ddxxyy(ii,i,j,2)=dxdy(ii,i,i,2)-dxdy(ii,j,j,2)
!             ddxxyy(ii,i,j,3)=dxdy(ii,i,i,3)-dxdy(ii,j,j,3)
!             fij = ddxxyy(ii,i,j,1)**2+ddxxyy(ii,i,j,2)**2+ddxxyy(ii,i,j,3)**2
!             fij = sqrt(fij)
!             write(*,'(4f20.15)')dxx(ii,i,j),xx(ii,1,i)-xx(ii,1,j),xx(ii,1,i),xx(ii,1,j)
!          enddo
!          endif
!             write(*,*)xx(ii,1,i)-xx(ii,1,j)
!          write(*,*)dxdy(ii,i,i,1),dxdy(ii,i,i,2),dxdy(ii,i,i,3)
!             enddo
!          enddo
!      enddo
      return
      end

!*********************************************************

      subroutine translate(nat,lattice,xyz,ldebug)
      implicit none

      integer i,j,k
      integer nat
      real*8 lattice(3,3)
      real*8 xyz(nat,3)
      real*8 axb(3),bxc(3),cxa(3)
      real*8 n1(3),n2(3),n3(3),n4(3),n5(3),n6(3)
      real*8 d1,d2,d3,d4,d5,d6
      real*8 distance1,distance2
      real*8 vorzeichen
      logical ldebug

      if(ldebug)then
        write(6,*)'lattice vectors (A)'
        do i=1,3
          write(6,'(3f14.6)')lattice(i,1),lattice(i,2),lattice(i,3)
        enddo
      endif

! calculation of the vector products
      axb(1)=lattice(1,2)*lattice(2,3)-lattice(1,3)*lattice(2,2)      
      axb(2)=lattice(1,3)*lattice(2,1)-lattice(1,1)*lattice(2,3)      
      axb(3)=lattice(1,1)*lattice(2,2)-lattice(1,2)*lattice(2,1)      

      bxc(1)=lattice(2,2)*lattice(3,3)-lattice(2,3)*lattice(3,2)      
      bxc(2)=lattice(2,3)*lattice(3,1)-lattice(2,1)*lattice(3,3)      
      bxc(3)=lattice(2,1)*lattice(3,2)-lattice(2,2)*lattice(3,1)      

      cxa(1)=lattice(3,2)*lattice(1,3)-lattice(3,3)*lattice(1,2)      
      cxa(2)=lattice(3,3)*lattice(1,1)-lattice(3,1)*lattice(1,3)      
      cxa(3)=lattice(3,1)*lattice(1,2)-lattice(3,2)*lattice(1,1)      
      
! plane (0,0,0),(a,0,0),(0,b,0)
! a x b through origin
      n1(1)=axb(1)/dsqrt(axb(1)**2+axb(2)**2+axb(3)**2)
      n1(2)=axb(2)/dsqrt(axb(1)**2+axb(2)**2+axb(3)**2)
      n1(3)=axb(3)/dsqrt(axb(1)**2+axb(2)**2+axb(3)**2)
      d1=0.d0 ! Plane through origin

! plane (0,0,c),(a,0,c),(0,b,c)
! parallel to a x b through c
      n2(1)=axb(1)/dsqrt(axb(1)**2+axb(2)**2+axb(3)**2)
      n2(2)=axb(2)/dsqrt(axb(1)**2+axb(2)**2+axb(3)**2)
      n2(3)=axb(3)/dsqrt(axb(1)**2+axb(2)**2+axb(3)**2)
      d2=-n2(1)*lattice(3,1)-n2(2)*lattice(3,2)-n2(3)*lattice(3,3)  

! plane (0,0,0),(0,b,0),(0,0,c)
! b x c through origin
      n3(1)=bxc(1)/dsqrt(bxc(1)**2+bxc(2)**2+bxc(3)**2)
      n3(2)=bxc(2)/dsqrt(bxc(1)**2+bxc(2)**2+bxc(3)**2)
      n3(3)=bxc(3)/dsqrt(bxc(1)**2+bxc(2)**2+bxc(3)**2)
      d3=0.d0

! plane (a,0,0),(a,b,0),(a,0,c)
! parallel to b x c trough a
      n4(1)=bxc(1)/dsqrt(bxc(1)**2+bxc(2)**2+bxc(3)**2)
      n4(2)=bxc(2)/dsqrt(bxc(1)**2+bxc(2)**2+bxc(3)**2)
      n4(3)=bxc(3)/dsqrt(bxc(1)**2+bxc(2)**2+bxc(3)**2)
      d4=-n4(1)*lattice(1,1)-n4(2)*lattice(1,2)-n4(3)*lattice(1,3)

! plane (0,0,0),(0,0,c),(a,0,0)
! c x a through origin
      n5(1)=cxa(1)/dsqrt(cxa(1)**2+cxa(2)**2+cxa(3)**2)
      n5(2)=cxa(2)/dsqrt(cxa(1)**2+cxa(2)**2+cxa(3)**2)
      n5(3)=cxa(3)/dsqrt(cxa(1)**2+cxa(2)**2+cxa(3)**2)
      d5=0.0d0

! plane (0,b,0),(0,b,c),(a,b,0)
! parallel to c x a through b
      n6(1)=cxa(1)/dsqrt(cxa(1)**2+cxa(2)**2+cxa(3)**2)
      n6(2)=cxa(2)/dsqrt(cxa(1)**2+cxa(2)**2+cxa(3)**2)
      n6(3)=cxa(3)/dsqrt(cxa(1)**2+cxa(2)**2+cxa(3)**2)
      d6=-n6(1)*lattice(2,1)-n6(2)*lattice(2,2)-n6(3)*lattice(2,3)

      do i=1,nat
 10     continue
        if(ldebug) then
          write(6,'(a,3f14.6)')'atom',xyz(i,1),xyz(i,2),xyz(i,3)
        endif

! calculate distance to plane 1          
        distance1=n1(1)*xyz(i,1)+n1(2)*xyz(i,2)+n1(3)*xyz(i,3)+d1
!        write(6,'(a,f14.6)')'distance to plane 1 ',distance1

! calculate distance to plane 2          
        distance2=n2(1)*xyz(i,1)+n2(2)*xyz(i,2)+n2(3)*xyz(i,3)+d2
!        write(6,'(a,f14.6)')'distance to plane 2 ',distance2
  
        vorzeichen=n1(1)*lattice(3,1)+n1(2)*lattice(3,2)&
                  +n1(3)*lattice(3,3)
        vorzeichen=vorzeichen/abs(vorzeichen) 
!        write(6,*)'Vorzeichen',vorzeichen 
        if((distance1.lt.-0.00001d0).and.(distance2.lt.-0.00001d0))then
          xyz(i,1)=xyz(i,1)+lattice(3,1)*vorzeichen
          xyz(i,2)=xyz(i,2)+lattice(3,2)*vorzeichen
          xyz(i,3)=xyz(i,3)+lattice(3,3)*vorzeichen
!          write(6,*)'Point is below plane 1'
!           write(6,*)'Atom is shifted into cell'
          goto 10
      elseif((distance1.gt.0.00001d0).and.(distance2.gt.0.00001d0))then
          xyz(i,1)=xyz(i,1)-lattice(3,1)*vorzeichen
          xyz(i,2)=xyz(i,2)-lattice(3,2)*vorzeichen
          xyz(i,3)=xyz(i,3)-lattice(3,3)*vorzeichen
!          write(6,*)'Point is above plane 2'
!           write(6,*)'Atom is shifted into cell'
          goto 10
        else
!          write(6,*)'Point is between planes 1 and 2'
        endif

! calculate distance to plane 3          
        distance1=n3(1)*xyz(i,1)+n3(2)*xyz(i,2)+n3(3)*xyz(i,3)+d3
!        write(6,'(a,f14.6)')'distance to plane 3 ',distance1

! calculate distance to plane 4          
        distance2=n4(1)*xyz(i,1)+n4(2)*xyz(i,2)+n4(3)*xyz(i,3)+d4
!        write(6,'(a,f14.6)')'distance to plane 4 ',distance2

        vorzeichen=n1(1)*lattice(3,1)+n1(2)*lattice(3,2)&
                  +n1(3)*lattice(3,3)
        vorzeichen=vorzeichen/abs(vorzeichen) 
!        write(6,*)'Vorzeichen',vorzeichen 
        if((distance1.lt.-0.00001d0).and.(distance2.lt.-0.00001d0))then
          xyz(i,1)=xyz(i,1)+lattice(1,1)*vorzeichen
          xyz(i,2)=xyz(i,2)+lattice(1,2)*vorzeichen
          xyz(i,3)=xyz(i,3)+lattice(1,3)*vorzeichen
!          write(6,*)'Point is below plane 3'
!           write(6,*)'Atom is shifted into cell'
          goto 10
      elseif((distance1.gt.0.00001d0).and.(distance2.gt.0.00001d0))then
          xyz(i,1)=xyz(i,1)-lattice(1,1)*vorzeichen
          xyz(i,2)=xyz(i,2)-lattice(1,2)*vorzeichen
          xyz(i,3)=xyz(i,3)-lattice(1,3)*vorzeichen
!          write(6,*)'Point is above plane 4'
!           write(6,*)'Atom is shifted into cell'
          goto 10
        else
!          write(6,*)'Point is between planes 3 and 4'
        endif

! calculate distance to plane 5          
        distance1=n5(1)*xyz(i,1)+n5(2)*xyz(i,2)+n5(3)*xyz(i,3)+d5
!        write(6,'(a,f14.6)')'distance to plane 5 ',distance1

! calculate distance to plane 6          
        distance2=n6(1)*xyz(i,1)+n6(2)*xyz(i,2)+n6(3)*xyz(i,3)+d6
!        write(6,'(a,f14.6)')'distance to plane 6 ',distance2

        vorzeichen=n1(1)*lattice(3,1)+n1(2)*lattice(3,2)&
                  +n1(3)*lattice(3,3)
        vorzeichen=vorzeichen/abs(vorzeichen)
!        write(6,*)'Vorzeichen',vorzeichen 
        if((distance1.lt.-0.00001d0).and.(distance2.lt.-0.00001d0))then
          xyz(i,1)=xyz(i,1)+lattice(2,1)*vorzeichen
          xyz(i,2)=xyz(i,2)+lattice(2,2)*vorzeichen
          xyz(i,3)=xyz(i,3)+lattice(2,3)*vorzeichen
!          write(6,*)'Point is below plane 5'
!           write(6,*)'Atom is shifted into cell'
          goto 10
      elseif((distance1.gt.0.00001d0).and.(distance2.gt.0.00001d0))then
          xyz(i,1)=xyz(i,1)-lattice(2,1)*vorzeichen
          xyz(i,2)=xyz(i,2)-lattice(2,2)*vorzeichen
          xyz(i,3)=xyz(i,3)-lattice(2,3)*vorzeichen
!          write(6,*)'Point is above plane 6'
!           write(6,*)'Atom is shifted into cell'
          goto 10
        else
!          write(6,*)'Point is between planes 5 and 6'
        endif




      enddo



      return
      end

!*********************************************************

      subroutine neighbor(nat,listdim,lsta,lstb,lstc,&
                 cutoff,lattice,xyz,lperiodic,ldebug,e,volume)
      implicit none
      integer nat,listdim,i,j,k
      integer indc, iat,jat
      integer lsta(2,nat)
      integer lstc(listdim)
      real*8 lstb(listdim,4)
      integer na,nb,nc
      integer n1,n2,n3
      real*8 xtemp,ytemp,ztemp
      real*8 alat,blat,clat
      real*8 xyz(nat,3)
      real*8 lattice(3,3)
      real*8 cutoff,rr
      real*8 xrel1,xrel2,xrel3,rr1
      real*8 dummy,e,volume,a1,a2,b1,b2,x0,n
      real*8 axb(3),axc(3),bxc(3) ! vector products
      real*8 absaxb,absaxc,absbxc
      real*8 proja,projb,projc
      logical lperiodic,ldebug

      if(lperiodic) then
! determine how often we have to multiply the cell in which direction
      na = 0
      nb = 0
      nc = 0
! calculate the norm vectors of the 3 planes
      axb(1)=lattice(1,2)*lattice(2,3)-lattice(1,3)*lattice(2,2)
      axb(2)=lattice(1,3)*lattice(2,1)-lattice(1,1)*lattice(2,3)
      axb(3)=lattice(1,1)*lattice(2,2)-lattice(1,2)*lattice(2,1)
      absaxb=sqrt(axb(1)**2+axb(2)**2+axb(3)**2)
      axb(1)=axb(1)/absaxb
      axb(2)=axb(2)/absaxb
      axb(3)=axb(3)/absaxb
      axc(1)=lattice(1,2)*lattice(3,3)-lattice(1,3)*lattice(3,2)
      axc(2)=lattice(1,3)*lattice(3,1)-lattice(1,1)*lattice(3,3)
      axc(3)=lattice(1,1)*lattice(3,2)-lattice(1,2)*lattice(3,1)
      absaxc=sqrt(axc(1)**2+axc(2)**2+axc(3)**2)
      axc(1)=axc(1)/absaxc
      axc(2)=axc(2)/absaxc
      axc(3)=axc(3)/absaxc
      bxc(1)=lattice(2,2)*lattice(3,3)-lattice(2,3)*lattice(3,2)
      bxc(2)=lattice(2,3)*lattice(3,1)-lattice(2,1)*lattice(3,3)
      bxc(3)=lattice(2,1)*lattice(3,2)-lattice(2,2)*lattice(3,1)
      absbxc=sqrt(bxc(1)**2+bxc(2)**2+bxc(3)**2)
      bxc(1)=bxc(1)/absbxc
      bxc(2)=bxc(2)/absbxc
      bxc(3)=bxc(3)/absbxc
! calculate the projections
      proja=lattice(1,1)*bxc(1)+lattice(1,2)*bxc(2)+lattice(1,3)*bxc(3)
      projb=lattice(2,1)*axc(1)+lattice(2,2)*axc(2)+lattice(2,3)*axc(3)
      projc=lattice(3,1)*axb(1)+lattice(3,2)*axb(2)+lattice(3,3)*axb(3)
      proja=abs(proja)
      projb=abs(projb)
      projc=abs(projc)
! determine how often we have to multiply the cell in which direction
      na = 0
      nb = 0
      nc = 0
      do while(dble(na)*proja.le.cutoff)
        na=na+1
      enddo
      do while(dble(nb)*projb.le.cutoff)
        nb=nb+1
      enddo
      do while(dble(nc)*projc.le.cutoff)
        nc=nc+1
      enddo
!      write(*,*)'proja,projb,projc',proja,projb,projc
!      write(*,*)'na,nb,nc',na,nb,nc
      endif

      volume=volume*(0.529177**3)
        a1=18144.5+8177.14*volume
        a2=651
        b1=4.42063+0.0934796*volume
        b2=1.41915-0.018669*volume
        x0=0.682428-0.029776*volume
        n=3.73178-0.0150275*volume
e=0
! determine lstb
      if(lperiodic)then
      indc=0
      do i=1,nat
        lsta(1,i)=indc+1
        do j=1,nat
          do n1=-na,na
            do n2=-nb,nb
              do n3=-nc,nc
! avoid interaction of atom with itself
         if((n1.eq.0).and.(n2.eq.0).and.(n3.eq.0).and.(i.eq.j)) then
         else
           xtemp=xyz(j,1)&
      +dble(n1)*lattice(1,1)+dble(n2)*lattice(2,1)+dble(n3)*lattice(3,1)
           ytemp=xyz(j,2)&
      +dble(n1)*lattice(1,2)+dble(n2)*lattice(2,2)+dble(n3)*lattice(3,2)
           ztemp=xyz(j,3)&
      +dble(n1)*lattice(1,3)+dble(n2)*lattice(2,3)+dble(n3)*lattice(3,3)
!           write(*,'(i4,3f14.6)')j,xtemp,ytemp,ztemp
           xrel1=xtemp-xyz(i,1)
           xrel2=ytemp-xyz(i,2)
           xrel3=ztemp-xyz(i,3)
           rr=xrel1**2+xrel2**2+xrel3**2
           if(rr.le.cutoff**2) then
             indc=indc+1
             if(indc.gt.listdim)then
                write(*,*)'Redimension lstb'
                stop
             endif
             lstb(indc,1)=xtemp
             lstb(indc,2)=ytemp
             lstb(indc,3)=ztemp
             lstb(indc,4)=sqrt(rr)
             rr1=sqrt(rr)
             e=e+a1*exp(-b1*rr1)/rr1-a2*cos(2*b2*(rr1-x0))/rr1**n
!             write(*,'(3i4,4f14.6)')indc,i,j,xtemp,ytemp,ztemp,sqrt(rr)
             lstc(indc)=j ! identification of atom
           endif
         endif
              enddo ! n3
            enddo ! n2
          enddo ! n1
        enddo ! j
        lsta(2,i)=indc
!        write(*,*)'atom',i,'from',lsta(1,i),'to',lsta(2,i)
      enddo ! i
      else ! .not. lperiodic
      indc=0
      do i=1,nat
        lsta(1,i)=indc+1
        do j=1,nat
         if(i.eq.j) then
         else
           xtemp=xyz(j,1)
           ytemp=xyz(j,2)
           ztemp=xyz(j,3)
           xrel1=xtemp-xyz(i,1)
           xrel2=ytemp-xyz(i,2)
           xrel3=ztemp-xyz(i,3)
           rr=xrel1**2+xrel2**2+xrel3**2
           if(rr.le.cutoff**2) then
             indc=indc+1
             if(indc.gt.listdim)then
                write(*,*)'Redimension lstb'
                stop
             endif
             lstb(indc,1)=xtemp
             lstb(indc,2)=ytemp
             lstb(indc,3)=ztemp
             lstb(indc,4)=sqrt(rr)
             lstc(indc)=j ! identification of atom
           endif
         endif
        enddo ! j
        lsta(2,i)=indc
!        write(*,*)'atom',i,'from',lsta(1,i),'to',lsta(2,i)
      enddo ! i
      endif
      return
      end

!*********************************************************
      subroutine neighborold(nat,lsta,lstb,cutoff,lattice,xxyyzz,ldebug)
      implicit none
      integer nat,i
      integer indc, iat,jat
      integer lsta(2,nat)
      integer lstb(24096)
      real*8 xxyyzz(27*nat,3)
      real*8 lattice(3,3)
      real*8 cutoff,rr
      real*8 xrel1,xrel2,xrel3
      real*8 dummy
      logical ldebug


! check if cell is larger than cutoff
      do i=1,3
        dummy=lattice(i,1)**2+lattice(i,2)**2+lattice(i,3)**2
        if(dummy.le.cutoff**2) then
          write(6,*)'Warning: cutoff is larger than lattice vector',i
          stop
        endif
      enddo

      indc=0
      do iat=1,nat
!        write(6,*)'starting atom ',iat
! starting position
        lsta(1,iat)=indc+1
        do jat=1,27*nat
          if(jat.ne.iat) then
            xrel1=xxyyzz(jat,1)-xxyyzz(iat,1)
            xrel2=xxyyzz(jat,2)-xxyyzz(iat,2)
            xrel3=xxyyzz(jat,3)-xxyyzz(iat,3)
            rr=xrel1**2+xrel2**2+xrel3**2
            if(rr.le.cutoff**2) then
              indc=indc+1
!              write(6,*)'close atom ',jat,rr
! check if estimated dimension is ok
              if(indc.gt.24096)then
                write(6,*)'Redimension lstb'
                stop
              endif
! nearest neighbor found
              lstb(indc)=jat
            endif
          endif
        enddo ! jat
! ending position
        lsta(2,iat)=indc
!        write(6,*)'lsta(2,iat)',lsta(2,iat)
      enddo ! iat

!      if(ldebug) then
      if(.false.) then
        write(6,*)'Neighbor list'
        do iat=1,nat
          write(6,*)'Atom ',iat,lsta(1,iat),lsta(2,iat)
          do jat=lsta(1,iat),lsta(2,iat)
            write(6,*)lstb(jat)
          enddo
        enddo
      endif

      return
      end

!********************************************************

      subroutine multiply(nat,lattice,xyz,xxyyzz,ldebug)
      implicit none
      integer i,j,k
      integer nat
      integer icount
      real*8 lattice(3,3)
      real*8 xyz(nat,3)
      real*8 xxyyzz(27*nat,3)
      logical ldebug

! copy original atoms
      do i=1,nat
        xxyyzz(i,1)=xyz(i,1)
        xxyyzz(i,2)=xyz(i,2)
        xxyyzz(i,3)=xyz(i,3)
      enddo

! a direction
      icount=0
      do i=nat+1,2*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)+lattice(1,1)
        xxyyzz(i,2)=xyz(icount,2)+lattice(1,2)
        xxyyzz(i,3)=xyz(icount,3)+lattice(1,3)
      enddo
! -a direction
      icount=0
      do i=2*nat+1,3*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)-lattice(1,1)
        xxyyzz(i,2)=xyz(icount,2)-lattice(1,2)
        xxyyzz(i,3)=xyz(icount,3)-lattice(1,3)
      enddo
! b direction
      icount=0
      do i=3*nat+1,4*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)+lattice(2,1)
        xxyyzz(i,2)=xyz(icount,2)+lattice(2,2)
        xxyyzz(i,3)=xyz(icount,3)+lattice(2,3)
      enddo
! -b direction
      icount=0
      do i=4*nat+1,5*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)-lattice(2,1)
        xxyyzz(i,2)=xyz(icount,2)-lattice(2,2)
        xxyyzz(i,3)=xyz(icount,3)-lattice(2,3)
      enddo
! c direction
      icount=0
      do i=5*nat+1,6*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)+lattice(3,1)
        xxyyzz(i,2)=xyz(icount,2)+lattice(3,2)
        xxyyzz(i,3)=xyz(icount,3)+lattice(3,3)
      enddo
! -c direction
      icount=0
      do i=6*nat+1,7*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)-lattice(3,1)
        xxyyzz(i,2)=xyz(icount,2)-lattice(3,2)
        xxyyzz(i,3)=xyz(icount,3)-lattice(3,3)
      enddo
! +a+b
      icount=0
      do i=7*nat+1,8*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)+lattice(1,1)+lattice(2,1)
        xxyyzz(i,2)=xyz(icount,2)+lattice(1,2)+lattice(2,2)
        xxyyzz(i,3)=xyz(icount,3)+lattice(1,3)+lattice(2,3)
      enddo
! +a-b
      icount=0
      do i=8*nat+1,9*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)+lattice(1,1)-lattice(2,1)
        xxyyzz(i,2)=xyz(icount,2)+lattice(1,2)-lattice(2,2)
        xxyyzz(i,3)=xyz(icount,3)+lattice(1,3)-lattice(2,3)
      enddo
! -a+b
      icount=0
      do i=9*nat+1,10*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)-lattice(1,1)+lattice(2,1)
        xxyyzz(i,2)=xyz(icount,2)-lattice(1,2)+lattice(2,2)
        xxyyzz(i,3)=xyz(icount,3)-lattice(1,3)+lattice(2,3)
      enddo
! -a-b
      icount=0
      do i=10*nat+1,11*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)-lattice(1,1)-lattice(2,1)
        xxyyzz(i,2)=xyz(icount,2)-lattice(1,2)-lattice(2,2)
        xxyyzz(i,3)=xyz(icount,3)-lattice(1,3)-lattice(2,3)
      enddo
! +a+c
      icount=0
      do i=11*nat+1,12*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)+lattice(1,1)+lattice(3,1)
        xxyyzz(i,2)=xyz(icount,2)+lattice(1,2)+lattice(3,2)
        xxyyzz(i,3)=xyz(icount,3)+lattice(1,3)+lattice(3,3)
      enddo
! +a-c
      icount=0
      do i=12*nat+1,13*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)+lattice(1,1)-lattice(3,1)
        xxyyzz(i,2)=xyz(icount,2)+lattice(1,2)-lattice(3,2)
        xxyyzz(i,3)=xyz(icount,3)+lattice(1,3)-lattice(3,3)
      enddo
! -a+c
      icount=0
      do i=13*nat+1,14*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)-lattice(1,1)+lattice(3,1)
        xxyyzz(i,2)=xyz(icount,2)-lattice(1,2)+lattice(3,2)
        xxyyzz(i,3)=xyz(icount,3)-lattice(1,3)+lattice(3,3)
      enddo
! -a-c
      icount=0
      do i=14*nat+1,15*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)-lattice(1,1)-lattice(3,1)
        xxyyzz(i,2)=xyz(icount,2)-lattice(1,2)-lattice(3,2)
        xxyyzz(i,3)=xyz(icount,3)-lattice(1,3)-lattice(3,3)
      enddo
! +b+c
      icount=0
      do i=15*nat+1,16*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)+lattice(2,1)+lattice(3,1)
        xxyyzz(i,2)=xyz(icount,2)+lattice(2,2)+lattice(3,2)
        xxyyzz(i,3)=xyz(icount,3)+lattice(2,3)+lattice(3,3)
      enddo
! +b-c
      icount=0
      do i=16*nat+1,17*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)+lattice(2,1)-lattice(3,1)
        xxyyzz(i,2)=xyz(icount,2)+lattice(2,2)-lattice(3,2)
        xxyyzz(i,3)=xyz(icount,3)+lattice(2,3)-lattice(3,3)
      enddo
! -b+c
      icount=0
      do i=17*nat+1,18*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)-lattice(2,1)+lattice(3,1)
        xxyyzz(i,2)=xyz(icount,2)-lattice(2,2)+lattice(3,2)
        xxyyzz(i,3)=xyz(icount,3)-lattice(2,3)+lattice(3,3)
      enddo
! -b-c
      icount=0
      do i=18*nat+1,19*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)-lattice(2,1)-lattice(3,1)
        xxyyzz(i,2)=xyz(icount,2)-lattice(2,2)-lattice(3,2)
        xxyyzz(i,3)=xyz(icount,3)-lattice(2,3)-lattice(3,3)
      enddo
! +a+b+c
      icount=0
      do i=19*nat+1,20*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)+lattice(1,1)+lattice(2,1)+lattice(3,1)
        xxyyzz(i,2)=xyz(icount,2)+lattice(1,2)+lattice(2,2)+lattice(3,2)
        xxyyzz(i,3)=xyz(icount,3)+lattice(1,3)+lattice(2,3)+lattice(3,3)
      enddo
! +a+b-c
      icount=0
      do i=20*nat+1,21*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)+lattice(1,1)+lattice(2,1)-lattice(3,1)
        xxyyzz(i,2)=xyz(icount,2)+lattice(1,2)+lattice(2,2)-lattice(3,2)
        xxyyzz(i,3)=xyz(icount,3)+lattice(1,3)+lattice(2,3)-lattice(3,3)
      enddo
! +a-b+c
      icount=0
      do i=21*nat+1,22*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)+lattice(1,1)-lattice(2,1)+lattice(3,1)
        xxyyzz(i,2)=xyz(icount,2)+lattice(1,2)-lattice(2,2)+lattice(3,2)
        xxyyzz(i,3)=xyz(icount,3)+lattice(1,3)-lattice(2,3)+lattice(3,3)
      enddo
! -a+b+c
      icount=0
      do i=22*nat+1,23*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)-lattice(1,1)+lattice(2,1)+lattice(3,1)
        xxyyzz(i,2)=xyz(icount,2)-lattice(1,2)+lattice(2,2)+lattice(3,2)
        xxyyzz(i,3)=xyz(icount,3)-lattice(1,3)+lattice(2,3)+lattice(3,3)
      enddo
! +a-b-c
      icount=0
      do i=23*nat+1,24*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)+lattice(1,1)-lattice(2,1)-lattice(3,1)
        xxyyzz(i,2)=xyz(icount,2)+lattice(1,2)-lattice(2,2)-lattice(3,2)
        xxyyzz(i,3)=xyz(icount,3)+lattice(1,3)-lattice(2,3)-lattice(3,3)
      enddo
! -a+b-c
      icount=0
      do i=24*nat+1,25*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)-lattice(1,1)+lattice(2,1)-lattice(3,1)
        xxyyzz(i,2)=xyz(icount,2)-lattice(1,2)+lattice(2,2)-lattice(3,2)
        xxyyzz(i,3)=xyz(icount,3)-lattice(1,3)+lattice(2,3)-lattice(3,3)
      enddo
! -a-b+c
      icount=0
      do i=25*nat+1,26*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)-lattice(1,1)-lattice(2,1)+lattice(3,1)
        xxyyzz(i,2)=xyz(icount,2)-lattice(1,2)-lattice(2,2)+lattice(3,2)
        xxyyzz(i,3)=xyz(icount,3)-lattice(1,3)-lattice(2,3)+lattice(3,3)
      enddo
! -a-b-c
      icount=0
      do i=26*nat+1,27*nat
        icount=icount+1
        xxyyzz(i,1)=xyz(icount,1)-lattice(1,1)-lattice(2,1)-lattice(3,1)
        xxyyzz(i,2)=xyz(icount,2)-lattice(1,2)-lattice(2,2)-lattice(3,2)
        xxyyzz(i,3)=xyz(icount,3)-lattice(1,3)-lattice(2,3)-lattice(3,3)
      enddo
!      if(ldebug) then
!        write(6,*)'Multiplied coordinates:'
!        do i=1,27*nat
!          write(6,'(i4,3f14.8)')i,xxyyzz(i,1),xxyyzz(i,2),xxyyzz(i,3)
!        enddo
!      endif

      return
      end

!******************************************************************



