program dumper

   use ieee_arithmetic
   implicit none

   character(len=999) :: fname
   integer, external :: iargc
   character(len=999) :: buf
   character(len=16) dataset, dataver
   character(len=64) datamod
   character(len=132) comment
   character(len=8) axtz,pmap,datum
   character(len=12) daten
   character(len=4) utmhem
   character(len=8) clab
   character(len=8), allocatable, dimension(:) :: clabs

   character(len=16), parameter :: prgname = 'calmet_dump'

   integer:: idh0, isec0, idh1, isec1
   integer :: novar
   character(len=8), dimension(100) :: ovar


   integer :: idbg, jdbg,kdbg
   integer :: nclab, ncheader
   logical :: lnoisy

   integer :: ncom
   integer ::ibyr,ibmo,ibdy,ibhr,ibsec,ieyr,iemo,iedy,iehr,iesec,irlg,irtype,nx,ny,nz,&
           nssta,nusta,npsta,nowsta,nlu,iwat1,iwat2,iutmzn
   real :: dgrid, xorigr,yorigr,feast,fnorth,rnlat0,relon0,xlat1,xlat2
   logical :: lcalgrd
   real, allocatable, dimension(:) :: zfacem,xssta,yssta,xusta,yusta,xpsta,ypsta
   real, allocatable, dimension(:,:) :: z0,elev,xlai,ustar,zi,el,wstar,rmm,tempk,rho,qsw
   integer, allocatable, dimension(:,:) :: ilandu,nears,ipgt,irh,ipcode
   real, allocatable, dimension(:,:,:) :: uu,vv,ww,tt

   integer,dimension(99) :: idum
   integer :: n,i,k,it,ios

   ! get file name

   lnoisy=.true.
   idbg=1
   jdbg=1
   kdbg=1
   idh0= 190000000
   idh1= 209936624
   isec0=0
   isec1=3600
   novar = -99
   !n = iargc()
   n = command_argument_count()
   if (n == 1) then
      !call getarg(1, fname)
      call get_command_argument(1, fname)
   elseif (n == 3) then
      !call getarg(1, fname)
      call get_command_argument(1, fname)
      read(fname, *, iostat=ios) idbg
      if (ios /= 0) then
         print*,'arg1 parse error: ', trim(fname)
         call usage
      endif
      !call getarg(2, fname)
      call get_command_argument(2, fname)
      read(fname, *, iostat=ios) jdbg
      if (ios /= 0) then
         print*,'arg2 parse error: ', trim(fname)
         call usage
      endif
      !call getarg(3, fname)
      call get_command_argument(3, fname)
   else
      call usage
   endif
   if (trim(fname) == '-') then
      ! read input from stdin
      call userin
   endif
   if (lnoisy) print*,trim(fname)
   open(11,file=fname,form='unformatted', status='old')


   ! scan entire file first, read for clabs
   call scn

   !print*,nclab,ncheader
   !print*,clabs(1:20)

   allocate(z0(nx,ny),ilandu(nx,ny),elev(nx,ny),xlai(nx,ny),nears(nx,ny))
   allocate(uu(nx,ny,nz),vv(nx,ny,nz),ww(nx,ny,nz),tt(nx,ny,nz))
   allocate(ipgt(nx,ny),ustar(nx,ny),zi(nx,ny),el(nx,ny),wstar(nx,ny),rmm(nx,ny), &
      tempk(nx,ny),rho(nx,ny),qsw(nx,ny),irh(nx,ny),ipcode(nx,ny))
   rewind(11)

   call dmp

contains 

   subroutine scn
      integer :: ic, ios, i
      ! header, metadata part
      read(11)dataset, dataver, datamod
      read(11)ncom
      if (lnoisy) print*,ncom
      if (ncom <= 10) then
         do i=1,ncom
            read(11) comment
            if (lnoisy) print*,comment
         enddo
      else
         do i=1,ncom
            read(11) comment
            if (lnoisy) then
               if (i <= 6 .or. i >=(ncom-2)) then
                  print*,comment
               elseif (i==7) then
                       print*
                       print*,'   ...   '
                       print*
               endif
            endif
         enddo
      endif

      ! header, scalar part
      read(11)ibyr,ibmo,ibdy,ibhr,ibsec,ieyr,iemo,iedy,iehr,iesec,axtz,&
              irlg,irtype,nx,ny,nz,dgrid,xorigr,yorigr,&
              nssta,nusta,npsta,nowsta,nlu,iwat1,iwat2,lcalgrd,pmap,datum,daten,&
              feast,fnorth,utmhem,iutmzn,rnlat0,relon0,xlat1,xlat2

      if (lnoisy) then
         print*
         print*,'i,jdbg  :', idbg,jdbg
         print*,'begin   :', ibyr,ibmo,ibdy,ibhr,ibsec
         print*,'end     :', ieyr,iemo,iedy,iehr,iesec
         print*,'irlg    :', irlg
         print*,'irtype  :', irlg
         print*,'nx/ny/nz:', nx,ny,nz
         print*,'xor/yor/dx:', xorigr,yorigr,dgrid
         print*
      endif

      ! i dont know why i cannot preduct # of records...
      !allocate(clabs(nx*ny*irlg*4+100))
      !allocate(clabs(99999))
      allocate(clabs(999999))
      !print*,nx*ny*irlg*4+100

      nclab = 0
      do while(.true.)
         read(11,iostat=ios)clab,idum(1:5)
         if (ios/=0) exit
         !print*,clab,idum(1:5)

         nclab = nclab + 1
         clabs(nclab) = clab
      enddo
      !print*,nclab

      do ic=1,nclab
         if (clabs(ic)(1:5)=='U-LEV') then
            ncheader = ic - 1
            exit
         endif
      enddo

   end subroutine

   subroutine dmp
      integer :: ic,ityp, ilab
      character(len=80) :: fmtf, fmti
      read(11)dataset, dataver, datamod
      read(11)ncom
      do i=1,ncom
         read(11) comment
      enddo

      ! header, scalar part
      read(11)ibyr,ibmo,ibdy,ibhr,ibsec,ieyr,iemo,iedy,iehr,iesec,axtz,&
              irlg,irtype,nx,ny,nz,dgrid,xorigr,yorigr,&
              nssta,nusta,npsta,nowsta,nlu,iwat1,iwat2,lcalgrd,pmap,datum,daten,&
              feast,fnorth,utmhem,iutmzn,rnlat0,relon0,xlat1,xlat2

      fmtf = "(a8,',',2(i9,',',i4,','),3(i4,','),e14.7)"
      fmti = "(a8,',',2(i9,',',i4,','),3(i4,','),i14)"
      ! data
      !print*,ncheader
      do ic=1,ncheader
         clab=clabs(ic)
         !print*,ic,clab
         if (clab == 'ZFACE') then
            read(11) clab, idum(1:4), zfacem
            if(lnoisy) print*, clab,idum(1:4),zfacem
         elseif (clab == 'XSSTA') then
            read(11)clab,idum(1:4),xssta
            if(lnoisy) print*, clab,idum(1:4)
         elseif (clab == 'YSSTA') then
            read(11)clab,idum(1:4),yssta
            if(lnoisy) print*, clab,idum(1:4)
         elseif (clab == 'XUSTA') then
            read(11)clab,idum(1:4),xusta
            if(lnoisy) print*, clab,idum(1:4)
         elseif (clab == 'YUSTA') then
            read(11)clab,idum(1:4),yusta
            if(lnoisy) print*, clab,idum(1:4)
         elseif (clab == 'XPSTA') then
            read(11)clab,idum(1:4),xpsta
            if(lnoisy) print*, clab,idum(1:4)
         elseif (clab == 'YPSTA') then
            read(11)clab,idum(1:4),ypsta
            if(lnoisy) print*, clab,idum(1:4)
         elseif (clab == 'Z0') then
            read(11)clab,idum(1:4),z0
            if(lnoisy) print*,clab, idum(1:4),z0(idbg,jdbg)
            if (novar == 0 .or. any(clab == ovar)) then
               write(*,fmtf) clab, 0, 0, 0, 0,idbg,jdbg,0, z0(idbg,jdbg)
            endif
         elseif (clab == 'ILANDU') then
            read(11)clab,idum(1:4),ilandu
            if(lnoisy) print*,clab, idum(1:4),ilandu(idbg,jdbg)
            if (novar == 0 .or. any(clab == ovar)) then
               write(*,fmti) clab, 0, 0, 0, 0,idbg,jdbg,0, ilandu(idbg,jdbg)
            endif
         elseif (clab == 'ELEV') then
            read(11)clab,idum(1:4),elev
            if(lnoisy) print*,clab, idum(1:4),elev(idbg,jdbg)
            if (novar == 0 .or. any(clab == ovar)) then
               write(*,fmtf) clab, 0, 0, 0, 0,idbg,jdbg,0, elev(idbg,jdbg)
            endif
         elseif (clab == 'XLAI') then
            read(11)clab,idum(1:4),xlai
            if(lnoisy) print*,clab, idum(1:4),xlai(idbg,jdbg)
            if (novar == 0 .or. any(clab == ovar)) then
               write(*,fmtf) clab, 0, 0, 0, 0,idbg,jdbg,0, xlai(idbg,jdbg)
            endif
         elseif (clab == 'NEARS') then
            read(11)clab,idum(1:4),nears
            if(lnoisy) print*,clab, idum(1:4),nears(idbg,jdbg)
            if (novar == 0 .or. any(clab == ovar)) then
               write(*,fmti) clab, 0, 0, 0, 0,idbg,jdbg,0, nears(idbg,jdbg)
            endif
         else
            print*, 'unknown clab (header):', clab
            stop
         endif



      enddo
      do ic=ncheader+1, nclab
         clab=clabs(ic)
         !print*,ic,clab
         if (clab(1:5) == 'U-LEV') then
            read(clab(6:8),*)k
            read(11)clab,idum(1:4),uu(:,:,k)
            if(lnoisy) print*, clab,idum(1:4),uu(idbg,jdbg,k)
            if (novar == 0 .or. any(clab == ovar)) then
               write(clab(6:8), '("   ")')
               if (k==kdbg .and. &
                  (idum(1) > idh0 .or. (idum(1) == idh0 .and. idum(2) >= isec0)) .and. &
                  (idum(3) < idh1 .or. (idum(3) == idh1 .and. idum(4) <= isec1)) ) then
                  write(*,fmtf) clab, idum(1:4), idbg,jdbg,k,uu(idbg,jdbg,k)
               endif
            endif
         elseif (clab(1:5) == 'V-LEV') then
            read(clab(6:8),*)k
            read(11)clab,idum(1:4),vv(:,:,k)
            if(lnoisy) print*, clab,idum(1:4),vv(idbg,jdbg,k)
            if (novar == 0 .or. any(clab == ovar)) then
               write(clab(6:8), '("   ")')
               if (k==kdbg .and. &
                  (idum(1) > idh0 .or. (idum(1) == idh0 .and. idum(2) >= isec0)) .and. &
                  (idum(3) < idh1 .or. (idum(3) == idh1 .and. idum(4) <= isec1)) ) then
                  write(*,fmtf) clab, idum(1:4), idbg,jdbg,k,vv(idbg,jdbg,k)
               endif
            endif
         elseif (clab(1:5) == 'WFACE') then
            read(clab(6:8),*)k
            read(11)clab,idum(1:4),ww(:,:,k)
            if(lnoisy) print*, clab,idum(1:4),ww(idbg,jdbg,k)
            if (novar == 0 .or. any(clab == ovar)) then
               write(clab(6:8), '("   ")')
               if (k==kdbg .and. &
                  (idum(1) > idh0 .or. (idum(1) == idh0 .and. idum(2) >= isec0)) .and. &
                  (idum(3) < idh1 .or. (idum(3) == idh1 .and. idum(4) <= isec1)) ) then
                  write(*,fmtf) clab, idum(1:4), idbg,jdbg,k,ww(idbg,jdbg,k)
               endif
            endif
         elseif (clab(1:5) == 'T-LEV') then
            read(clab(6:8),*)k
            read(11)clab,idum(1:4),tt(:,:,k)
            if(lnoisy) print*, clab,idum(1:4),tt(idbg,jdbg,k)
            if (novar == 0 .or. any(clab == ovar)) then
               write(clab(6:8), '("   ")')
               if (k==kdbg .and. &
                  (idum(1) > idh0 .or. (idum(1) == idh0 .and. idum(2) >= isec0)) .and. &
                  (idum(3) < idh1 .or. (idum(3) == idh1 .and. idum(4) <= isec1)) ) then
                  write(*,fmtf) clab, idum(1:4), idbg,jdbg,k,tt(idbg,jdbg,k)
               endif
            endif
         elseif (clab(1:5) == 'IPGT') then
            read(11)clab,idum(1:4),ipgt
            k = 1
            if(lnoisy) print*, clab,idum(1:4),ipgt(idbg,jdbg)
            if (novar == 0 .or. any(clab == ovar)) then
               if (k==kdbg .and. &
                  (idum(1) > idh0 .or. (idum(1) == idh0 .and. idum(2) >= isec0)) .and. &
                  (idum(3) < idh1 .or. (idum(3) == idh1 .and. idum(4) <= isec1)) ) then
                  write(*,fmti) clab, idum(1:4), idbg,jdbg,k,ipgt(idbg,jdbg)
               endif
            endif
         elseif (clab(1:5) == 'USTAR') then
            read(11)clab,idum(1:4),ustar
            k = 1
            if(lnoisy) print*, clab,idum(1:4),ustar(idbg,jdbg)
            if (novar == 0 .or. any(clab == ovar)) then
               if (k==kdbg .and. &
                  (idum(1) > idh0 .or. (idum(1) == idh0 .and. idum(2) >= isec0)) .and. &
                  (idum(3) < idh1 .or. (idum(3) == idh1 .and. idum(4) <= isec1)) ) then
                  write(*,fmtf) clab, idum(1:4), idbg,jdbg,k,ustar(idbg,jdbg)
               endif
            endif
         elseif (clab(1:5) == 'ZI') then
            read(11)clab,idum(1:4),zi
            k = 1
            if(lnoisy) print*, clab,idum(1:4),zi(idbg,jdbg)
            if (novar == 0 .or. any(clab == ovar)) then
               if (k==kdbg .and. &
                  (idum(1) > idh0 .or. (idum(1) == idh0 .and. idum(2) >= isec0)) .and. &
                  (idum(3) < idh1 .or. (idum(3) == idh1 .and. idum(4) <= isec1)) ) then
                  write(*,fmtf) clab, idum(1:4), idbg,jdbg,k,zi(idbg,jdbg)
               endif
            endif
         elseif (clab(1:5) == 'EL') then
            read(11)clab,idum(1:4),el
            k = 1
            if(lnoisy) print*, clab,idum(1:4),el(idbg,jdbg)
            if (novar == 0 .or. any(clab == ovar)) then
               if (k==kdbg .and. &
                  (idum(1) > idh0 .or. (idum(1) == idh0 .and. idum(2) >= isec0)) .and. &
                  (idum(3) < idh1 .or. (idum(3) == idh1 .and. idum(4) <= isec1)) ) then
                  write(*,fmtf) clab, idum(1:4), idbg,jdbg,k,el(idbg,jdbg)
               endif
            endif
         elseif (clab(1:5) == 'WSTAR') then
            read(11)clab,idum(1:4),wstar
            k = 1
            if(lnoisy) print*, clab,idum(1:4),wstar(idbg,jdbg)
            if (novar == 0 .or. any(clab == ovar)) then
               if (k==kdbg .and. &
                  (idum(1) > idh0 .or. (idum(1) == idh0 .and. idum(2) >= isec0)) .and. &
                  (idum(3) < idh1 .or. (idum(3) == idh1 .and. idum(4) <= isec1)) ) then
                  write(*,fmtf) clab, idum(1:4), idbg,jdbg,k,wstar(idbg,jdbg)
               endif
            endif
         elseif (clab(1:5) == 'RMM') then
            read(11)clab,idum(1:4),rmm
            k = 1
            if(lnoisy) print*, clab,idum(1:4),rmm(idbg,jdbg)
            if (novar == 0 .or. any(clab == ovar)) then
               if (k==kdbg .and. &
                  (idum(1) > idh0 .or. (idum(1) == idh0 .and. idum(2) >= isec0)) .and. &
                  (idum(3) < idh1 .or. (idum(3) == idh1 .and. idum(4) <= isec1)) ) then
                  write(*,fmtf) clab, idum(1:4), idbg,jdbg,k,rmm(idbg,jdbg)
               endif
            endif
         elseif (clab(1:5) == 'TEMPK') then
            read(11)clab,idum(1:4),tempk
            k = 1
            if(lnoisy) print*, clab,idum(1:4),tempk(idbg,jdbg)
            if (novar == 0 .or. any(clab == ovar)) then
               if (k==kdbg .and. &
                  (idum(1) > idh0 .or. (idum(1) == idh0 .and. idum(2) >= isec0)) .and. &
                  (idum(3) < idh1 .or. (idum(3) == idh1 .and. idum(4) <= isec1)) ) then
                  write(*,fmtf) clab, idum(1:4), idbg,jdbg,k,tempk(idbg,jdbg)
               endif
            endif
         elseif (clab(1:5) == 'RHO') then
            read(11)clab,idum(1:4),rho
            k = 1
            if(lnoisy) print*, clab,idum(1:4),rho(idbg,jdbg)
            if (novar == 0 .or. any(clab == ovar)) then
               if (k==kdbg .and. &
                  (idum(1) > idh0 .or. (idum(1) == idh0 .and. idum(2) >= isec0)) .and. &
                  (idum(3) < idh1 .or. (idum(3) == idh1 .and. idum(4) <= isec1)) ) then
                  write(*,fmtf) clab, idum(1:4), idbg,jdbg,k,rho(idbg,jdbg)
               endif
            endif
         elseif (clab(1:5) == 'QSW') then
            read(11)clab,idum(1:4),qsw
            k = 1
            if(lnoisy) print*, clab,idum(1:4),qsw(idbg,jdbg)
            if (novar == 0 .or. any(clab == ovar)) then
               if (k==kdbg .and. &
                  (idum(1) > idh0 .or. (idum(1) == idh0 .and. idum(2) >= isec0)) .and. &
                  (idum(3) < idh1 .or. (idum(3) == idh1 .and. idum(4) <= isec1)) ) then
                  write(*,fmtf) clab, idum(1:4), idbg,jdbg,k,qsw(idbg,jdbg)
               endif
            endif
         elseif (clab(1:5) == 'IRH') then
            read(11)clab,idum(1:4),irh
            k = 1
            if(lnoisy) print*, clab,idum(1:4),irh(idbg,jdbg)
            if (novar == 0 .or. any(clab == ovar)) then
               if (k==kdbg .and. &
                  (idum(1) > idh0 .or. (idum(1) == idh0 .and. idum(2) >= isec0)) .and. &
                  (idum(3) < idh1 .or. (idum(3) == idh1 .and. idum(4) <= isec1)) ) then
                  write(*,fmti) clab, idum(1:4), idbg,jdbg,k,irh(idbg,jdbg)
               endif
            endif
         elseif (clab(1:6) == 'IPCODE') then
            read(11)clab,idum(1:4),ipcode
            k = 1
            if(lnoisy) print*, clab,idum(1:4),ipcode(idbg,jdbg)
            if (novar == 0 .or. any(clab == ovar)) then
               if (k==kdbg .and. &
                  (idum(1) > idh0 .or. (idum(1) == idh0 .and. idum(2) >= isec0)) .and. &
                  (idum(3) < idh1 .or. (idum(3) == idh1 .and. idum(4) <= isec1)) ) then
                  write(*,fmti) clab, idum(1:4), idbg,jdbg,k,ipcode(idbg,jdbg)
               endif
            endif
         else
            print*, 'unknown clab (data):', ic, clab
            stop
         endif

         
      enddo

   end subroutine

   subroutine userin
      character(len=999) s
      lnoisy = .false.
      write(0,'(a)', advance='no')'fname: '
      read(*, '(a)') fname
      write(0, *) trim(fname)

      write(0,'(a)', advance='no')'i,j: '
      s = ''
      read(*, '(a)', iostat=ios) s
      if (ios /= 0 .or. s == '') then
         write(0, *) 'not specified, using default'
      else
         read(s, *, iostat=ios) idbg, jdbg
         if (ios /= 0 ) then
            write(0,*) 'cannot parse ijdx: ', trim(s)
            stop
         endif
      endif
      write(0, *) idbg, jdbg

      write(0,'(a)', advance='no') 'btime (yyyyjjjhh ssss): '
      s = ''
      read(*, '(a)', iostat=ios) s
      if (ios /= 0 .or. s == '') then
         write(0, *) 'not specified, using default'
      else
         read(s, *, iostat=ios) idh0, isec0
         if (ios /= 0 ) then
            write(0,*) 'cannot parse btime: ', trim(s)
            stop
         endif
      endif
      write(0,*) idh0, isec0

      write(0,'(a)', advance='no') 'etime (yyyyjjjhh ssss): '
      s = ''
      read(*, '(a)', iostat=ios) s
      if (ios /= 0 .or. s == '') then
         write(0, *) 'not specified, using default'
      else
         read(s, *, iostat=ios) idh1,isec1
         if (ios /= 0 ) then
            write(0,*) 'cannot parse etime: ', trim(s)
            stop
         endif
      endif
      write(0,*) idh1, isec1

      novar = 0
      write(0,'(a)', advance='no') 'variables (blank to all): '
      do while(.true.)
        read(*, '(a)', iostat=ios) s
        if (ios/=0) then
           s= ''
        endif
        if (s=='') then
           if (novar==0) then
              write(0, *) 'not specified, output all var'
           else
              write(0, *) 'finished, n=', novar
           endif
           exit
        endif
        novar=novar+1
        ovar(novar) = trim(s)
        write(0, '(a)', advance='no') ' more variable (blank to end): '
      enddo

      ! filename
      ! i,j cell to dump
      ! time to dump
      ! variable to dump

   end subroutine

   subroutine usage()
      print '(99a)', 'usage:  ', trim(prgname) , ' [ idbg jdbg ] _m3d_file_'
      print '(99a)', '        ', trim(prgname) , ' -'
      print*
      print*, 'first method dumps values of all variables to screan, more of debugging output'
      print*, 'second method prompt further input from STDIN, and create csv format output'
      stop
   end subroutine

end program

! vim: et sw=3
