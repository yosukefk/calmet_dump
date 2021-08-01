program dumper

   use ieee_arithmetic
   implicit none

   character(len=999) :: fname, oname
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
   logical :: lncf


   integer :: ncom
   integer ::ibyr,ibmo,ibdy,ibhr,ibsec,ieyr,iemo,iedy,iehr,iesec,irlg,irtype,nx,ny,nz,&
           nssta,nusta,npsta,nowsta,nlu,iwat1,iwat2,iutmzn, iwfcod
   real :: dgrid, xorigr,yorigr,feast,fnorth,rnlat0,relon0,xlat1,xlat2
   logical :: lcalgrd
   real, allocatable, dimension(:) :: zfacem,xssta,yssta,xusta,yusta,xpsta,ypsta
   real, allocatable, dimension(:,:) :: z0,elev,xlai,ustar,zi,el,wstar,rmm,tempk,rho,qsw
   integer, allocatable, dimension(:,:) :: ilandu,nears,ipgt,irh,ipcode
   real, allocatable, dimension(:,:,:) :: uu,vv,ww,tt

   integer,dimension(99) :: idum
   integer :: n,i,k,it,ios,ipos

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
   lncf = .false.
   !n = iargc()
   n = command_argument_count()
   if (n == 1) then
      !call getarg(1, fname)
      call get_command_argument(1, fname)

   elseif (n <= 3) then
      !call getarg(1, fname)
      call get_command_argument(1, fname)

      print*,'"' // fname // '"'
      print*,'"' // fname(1:1) // '"'
      print*,fname(1:1)=='-'

      ! case 1, -ncf fname oname
      if (fname(1:1) == '-') then
         if (fname(1:2)=='-n' .or. fname(1:2)=='-N' .or. fname(1:3)=='--n' .or. fname(1:3) =='--N') then
            lncf = .true.
         else
            print*,'unknown option: ', trim(fname)
            call usage
         endif

         call get_command_argument(2, fname)
         if (n > 2) then
            call get_command_argument(3, oname)
         else
            ipos=len_trim(fname)-3
            print*,'ipos',ipos
            print*,fname(ipos:ipos)
            if (fname(ipos:ipos)=='.') then
               oname = fname(1:ipos) // 'nc'
            else
               oname = fname(1:len_trim(fname)) // '.nc'
            endif
         endif

      else

         ! case 2 idx jdx
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
      endif
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

   allocate(zfacem(nz+1))
   allocate(z0(nx,ny),ilandu(nx,ny),elev(nx,ny),xlai(nx,ny),nears(nx,ny))
   allocate(uu(nx,ny,nz),vv(nx,ny,nz),ww(nx,ny,nz+1),tt(nx,ny,nz))
   allocate(ipgt(nx,ny),ustar(nx,ny),zi(nx,ny),el(nx,ny),wstar(nx,ny),rmm(nx,ny), &
      tempk(nx,ny),rho(nx,ny),qsw(nx,ny),irh(nx,ny),ipcode(nx,ny))
   rewind(11)

   if (lncf) then
      print*,'here'
      call ncf(oname)
   else
      call dmp
   endif

contains 

   subroutine ncf(oname)
      ! https://www.unidata.ucar.edu/software/netcdf/examples/programs/
      use netcdf
      implicit none
      character(len=999) :: oname
      integer :: ncid
      character(len=*), parameter :: lvl_name = 'Level'
      character(len=*), parameter :: zface_name = 'Zface'
      character(len=*), parameter :: yyy_name = 'Y'
      character(len=*), parameter :: xxx_name = 'X'
      character(len=*), parameter :: lon_name = 'longitude'
      character(len=*), parameter :: lat_name = 'latitude'
      character(len=*), parameter :: tim_name = 'Time'
      character(len=*), parameter :: tnchar_name = 'DateStrLen'

      integer :: lvl_dimid, yyy_dimid, xxx_dimid, tim_dimid, zface_dimid, tnchar_dimid, lon_varid, lat_varid

      integer :: zface_varid, xxx_varid, yyy_varid, tim_varid

      character(len=*), parameter :: units = 'units'
      character(len=*), parameter :: xxx_units = 'm'
      character(len=*), parameter :: yyy_units = 'm'
      character(len=*), parameter :: lon_units = 'degree_east'
      character(len=*), parameter :: lat_units = 'degree_north'

      integer :: z0_varid, lu_varid, elev_varid, lai_varid

      character(len=*), parameter :: uu_name = 'U'
      character(len=*), parameter :: vv_name = 'V'
      character(len=*), parameter :: ww_name = 'W'
      character(len=*), parameter :: tt_name = 'T'
      integer :: uu_varid, vv_varid, ww_varid, tt_varid
      integer :: ipgt_varid, ust_varid, zi_varid, el_varid, wst_varid
      integer :: rmm_varid, tempk_varid, rho_varid, qsw_varid, irh_varid, ipcode_varid

      character(len=*), parameter :: description = 'description'

      real, dimension(:), allocatable :: xxx
      real, dimension(:), allocatable :: yyy
      real, dimension(:,:), allocatable :: lon
      real, dimension(:,:), allocatable :: lat

      integer :: ic, ii, jj
      character(len=8) :: clab0
      character(len=19) :: tstamp, lst_tstamp
      integer :: tdx, lst_tdx
      integer :: status
      real :: rdum, tmsone
      integer :: idum1
      character(len=12) :: caction
      real(kind=8), dimension(9) :: vecti, vecto
      character(len=8) :: pmapo
      character(len=4) :: utmhem2
      integer :: iutmzn2


      print*,oname

      read(11)dataset, dataver, datamod
      read(11)ncom
      do i=1,ncom
         read(11) comment
      enddo

      ! header, scalar part
      read(11)ibyr,ibmo,ibdy,ibhr,ibsec,ieyr,iemo,iedy,iehr,iesec,axtz,&
              irlg,irtype,nx,ny,nz,dgrid,xorigr,yorigr,&
              iwfcod,nssta,nusta,npsta,nowsta,nlu,iwat1,iwat2,lcalgrd,pmap,datum,daten,&
              feast,fnorth,utmhem,iutmzn,rnlat0,relon0,xlat1,xlat2

      pmapo = 'LL'

      tmsone =  1.
      print*,axtz
      print*, xorigr, yorigr
      print*,pmap, datum, rnlat0, relon0
      print*,pmap, iutmzn, tmsone, xlat1, xlat2, rnlat0, relon0, feast, fnorth
      print*
      print*, vecti
      print*
      print*, vecto
      call globe1(pmap, iutmzn, tmsone,xlat1, xlat2, rnlat0, relon0,  &
                  feast, fnorth, &
                  pmapo, idum1, rdum, rdum, rdum, rdum, rdum, &
                  rdum, rdum, &
                  caction, vecti, vecto)

      ! create file
      call check( nf90_create(oname, nf90_clobber, ncid) )

      ! define dims
      call check( nf90_def_dim( ncid, tim_name, nf90_unlimited, tim_dimid) )
      call check( nf90_def_dim( ncid, tnchar_name, 19, tnchar_dimid) )
      call check( nf90_def_dim( ncid, lvl_name, nz, lvl_dimid) )
      call check( nf90_def_dim( ncid, yyy_name, ny, yyy_dimid) )
      call check( nf90_def_dim( ncid, xxx_name, nx, xxx_dimid) )
      call check( nf90_def_dim( ncid, zface_name, nz+1, zface_dimid) )

      ! define coord vars
      call check( nf90_def_var( ncid, lon_name, nf90_real, (/xxx_dimid, yyy_dimid/), lon_varid) )
      call check( nf90_def_var( ncid, lat_name, nf90_real, (/xxx_dimid, yyy_dimid/), lat_varid) )
      call check( nf90_def_var( ncid, xxx_name, nf90_real, xxx_dimid, xxx_varid))
      call check( nf90_def_var( ncid, yyy_name, nf90_real, yyy_dimid, yyy_varid))
      call check( nf90_def_var( ncid, tim_name, nf90_char, (/tnchar_dimid, tim_dimid/),  tim_varid))

      call check( nf90_put_att( ncid, xxx_varid, units, xxx_units) )
      call check( nf90_put_att( ncid, yyy_varid, units, yyy_units) )
      call check( nf90_put_att( ncid, lon_varid, units, lon_units) )
      call check( nf90_put_att( ncid, lat_varid, units, lat_units) )

      allocate(xxx(nx), yyy(ny), lon(nx, ny), lat(nx,ny))
      do ii = 1, nx
         xxx(ii) = xorigr + (ii-1) * dgrid + .5 * dgrid
      enddo
      do jj = 1, ny
         yyy(jj) = yorigr + (jj-1) * dgrid + .5 * dgrid
      enddo

      do ii=1,nx
         do jj=1,ny
            call globe(6, caction, datum, vecti, datum, vecto, &
               xxx(ii)*.001, yyy(jj)*.001, lon(ii,jj), lat(ii,jj), iutmzn2, utmhem2)
         enddo
      enddo


      do ic=1,ncheader
         clab=clabs(ic)
         print*,clab
         if (clab=='ZFACE') then
            call check( nf90_def_var( ncid, trim(clab), nf90_real, zface_dimid, zface_varid))
            call check( nf90_put_att( ncid, zface_varid, units, 'm') )
            call check( nf90_put_att( ncid, zface_varid, description, 'heights of cell faces') )
         elseif (clab=='Z0') then
            call check( nf90_def_var( ncid, trim((clab)), nf90_real, (/xxx_dimid, yyy_dimid/), z0_varid))
            call check( nf90_put_att( ncid, z0_varid, units, 'm') )
            call check( nf90_put_att( ncid, z0_varid, description, 'surface roughness lengths') )
         elseif (clab=='ILANDU') then
            call check( nf90_def_var( ncid, trim((clab)), nf90_int, (/xxx_dimid, yyy_dimid/), lu_varid))
            call check( nf90_put_att( ncid, lu_varid, description, 'land use category') )
         elseif (clab=='ELEV') then
            call check( nf90_def_var( ncid, trim((clab)), nf90_real, (/xxx_dimid, yyy_dimid/), elev_varid))
            call check( nf90_put_att( ncid, elev_varid, units, 'm') )
            call check( nf90_put_att( ncid, elev_varid, description, 'terrain elevations') )
         elseif (clab=='XLAI') then
            call check( nf90_def_var( ncid, trim((clab)), nf90_real, (/xxx_dimid, yyy_dimid/), lai_varid))
            call check( nf90_put_att( ncid, lai_varid, units, 'm/m') )
            call check( nf90_put_att( ncid, lai_varid, description, 'leaf area index') )
         endif

      enddo
      clab0=''
      uu_varid = 0
      vv_varid = 0
      ww_varid = 0
      tt_varid = 0
      ipgt_varid = 0
      ust_varid = 0
      zi_varid = 0
      el_varid = 0
      wst_varid = 0
      rmm_varid = 0
      tempk_varid = 0
      rho_varid = 0
      qsw_varid = 0
      irh_varid = 0
      ipcode_varid = 0
      do ic=ncheader+1, nclab
         clab=clabs(ic)
         print*,clab
         if (ic == ncheader+1) then
            clab0 = clab
         else
            if (clab==clab0) then
               exit
            endif
         endif



         if (clab(1:5)=='U-LEV') then
            !call check(  nf90_def_var( ncid, uu_name, nf90_real, (/xxx_dimid, yyy_dimid, lvl_dimid, tim_dimid/), uu_varid))
            status =  nf90_def_var( ncid, uu_name, nf90_real, (/xxx_dimid, yyy_dimid, lvl_dimid, tim_dimid/), uu_varid)
            !status = nf90_def_var( ncid, uu_name, nf90_real, (/tim_dimid, lvl_dimid, yyy_dimid, xxx_dimid/), uu_varid)
            call check( nf90_put_att( ncid, uu_varid, units, 'm/s') )
            call check( nf90_put_att( ncid, uu_varid, description, 'U-component of the winds') )
            
         elseif (clab(1:5)=='V-LEV') then
            status =  nf90_def_var( ncid, vv_name, nf90_real, (/xxx_dimid, yyy_dimid, lvl_dimid, tim_dimid/), vv_varid)
            call check( nf90_put_att( ncid, vv_varid, units, 'm/s') )
            call check( nf90_put_att( ncid, vv_varid, description, 'V-component of the winds') )
            
         elseif (clab(1:5)=='WFACE') then
            status =  nf90_def_var( ncid, ww_name, nf90_real, (/xxx_dimid, yyy_dimid, zface_dimid, tim_dimid/), ww_varid)
            call check( nf90_put_att( ncid, ww_varid, units, 'm/s') )
            call check( nf90_put_att( ncid, ww_varid, description, 'W-component of the winds') )
         elseif (clab(1:5)=='T-LEV') then
            status =  nf90_def_var( ncid, tt_name, nf90_real, (/xxx_dimid, yyy_dimid, lvl_dimid, tim_dimid/), tt_varid)
            call check( nf90_put_att( ncid, tt_varid, units, 'K') )
            call check( nf90_put_att( ncid, tt_varid, description, 'air temperature') )
         elseif (clab=='IPGT') then
            call check( nf90_def_var( ncid, trim((clab)), nf90_int, (/xxx_dimid, yyy_dimid, tim_dimid/), ipgt_varid) )
            call check( nf90_put_att( ncid, ipgt_varid, description, 'PGT stability class') )
         elseif (clab=='USTAR') then
            call check( nf90_def_var( ncid, trim(clab), nf90_real, (/xxx_dimid, yyy_dimid, tim_dimid/), ust_varid) )
            call check( nf90_put_att( ncid, ust_varid, units, 'm/s') )
            call check( nf90_put_att( ncid, ust_varid, description, 'surface friction velocity') )
         elseif (clab=='ZI') then
            call check( nf90_def_var( ncid, trim(clab), nf90_real, (/xxx_dimid, yyy_dimid, tim_dimid/), zi_varid) )
            call check( nf90_put_att( ncid, zi_varid, units, 'm') )
            call check( nf90_put_att( ncid, zi_varid, description, 'mixing height') )
         elseif (clab=='EL') then
            call check( nf90_def_var( ncid, trim(clab), nf90_real, (/xxx_dimid, yyy_dimid, tim_dimid/), el_varid) )
            call check( nf90_put_att( ncid, el_varid, units, 'm') )
            call check( nf90_put_att( ncid, el_varid, description, 'Monin-Obukhov length') )
         elseif (clab=='WSTAR') then
            call check( nf90_def_var( ncid, trim(clab), nf90_real, (/xxx_dimid, yyy_dimid, tim_dimid/), wst_varid) )
            call check( nf90_put_att( ncid, wst_varid, units, 'm/s') )
            call check( nf90_put_att( ncid, wst_varid, description, 'convective velocity scale') )
         elseif (clab=='RMM') then
            call check( nf90_def_var( ncid, trim(clab), nf90_real, (/xxx_dimid, yyy_dimid, tim_dimid/), rmm_varid) )
            call check( nf90_put_att( ncid, rmm_varid, units, 'mm/hr') )
            call check( nf90_put_att( ncid, rmm_varid, description, 'precipitation rate') )
         elseif (clab=='TEMPK') then
            call check( nf90_def_var( ncid, trim(clab), nf90_real, (/xxx_dimid, yyy_dimid, tim_dimid/), tempk_varid) )
            call check( nf90_put_att( ncid, tempk_varid, units, 'K') )
            call check( nf90_put_att( ncid, tempk_varid, description, 'near-surface temperature') )
         elseif (clab=='RHO') then
            call check( nf90_def_var( ncid, trim(clab), nf90_real, (/xxx_dimid, yyy_dimid, tim_dimid/), rho_varid) )
            call check( nf90_put_att( ncid, rho_varid, units, 'kg/m3') )
            call check( nf90_put_att( ncid, rho_varid, description, 'near-surface air density') )
         elseif (clab=='QSW') then
            call check( nf90_def_var( ncid, trim(clab), nf90_real, (/xxx_dimid, yyy_dimid, tim_dimid/), qsw_varid) )
            call check( nf90_put_att( ncid, qsw_varid, units, 'W/m2') )
            call check( nf90_put_att( ncid, qsw_varid, description, 'short-wave solar radiation') )
         elseif (clab=='IRH') then
            call check( nf90_def_var( ncid, trim(clab), nf90_int, (/xxx_dimid, yyy_dimid, tim_dimid/), irh_varid) )
            call check( nf90_put_att( ncid, irh_varid, units, '%') )
            call check( nf90_put_att( ncid, irh_varid, description, 'near-surface relative humidity') )
         elseif (clab=='IPCODE') then
            call check( nf90_def_var( ncid, trim(clab), nf90_int, (/xxx_dimid, yyy_dimid, tim_dimid/), ipcode_varid) )
            call check( nf90_put_att( ncid, ipcode_varid, description, 'precipitation type code') )
         endif
      enddo

      call check( nf90_enddef(ncid) )

      call check( nf90_put_var( ncid, xxx_varid, xxx ) )
      call check( nf90_put_var( ncid, yyy_varid, yyy ) )
      call check( nf90_put_var( ncid, lon_varid, lon ) )
      call check( nf90_put_var( ncid, lat_varid, lat ) )

      do ic=1,ncheader
         clab=clabs(ic)
         if (clab=='ZFACE') then
            read(11) clab, idum(1:4), zfacem
            call check( nf90_put_var( ncid, zface_varid, zfacem ) )
         elseif (clab=='Z0') then
            read(11)clab,idum(1:4),z0
            call check( nf90_put_var( ncid, z0_varid, z0 ) )
         elseif (clab=='ILANDU') then
            read(11)clab,idum(1:4),ilandu
            call check( nf90_put_var( ncid, lu_varid, ilandu ) )
         elseif (clab=='ELEV') then
            read(11)clab,idum(1:4),elev
            call check( nf90_put_var( ncid, elev_varid, elev ) )
         elseif (clab=='XLAI') then
            read(11)clab,idum(1:4),xlai
            call check( nf90_put_var( ncid, lai_varid, xlai ) )
         endif
      enddo

      lst_tstamp = ''
      lst_tdx = 0
      do ic=ncheader+1, nclab
         clab=clabs(ic)
         print*,clab
         if (clab(1:5) == 'U-LEV') then
            read(clab(6:8),*)k
            read(11)clab,idum(1:4),uu(:,:,k)
            !print*,clab, idum(1:2)
            call cmtime_to_tstamp(idum(1), idum(2), tstamp)

            if (tstamp /= lst_tstamp) then
               !newtime
               print*,'newtime'
               tdx = lst_tdx + 1

               if (lst_tstamp /= '') then
                  print*, 'tim', tstamp
                  call check (nf90_put_var( ncid, tim_varid, tstamp, start=(/1, lst_tdx/)))
                  ! write values
                  if (uu_varid /= 0) then
                     call check (nf90_put_var( ncid, uu_varid, uu, start=(/1,1,1, lst_tdx/)) )
                  endif
                  if (vv_varid /= 0) then
                     call check (nf90_put_var( ncid, vv_varid, vv, start=(/1,1,1, lst_tdx/)) )
                  endif
                  if (ww_varid /= 0) then
                     call check (nf90_put_var( ncid, ww_varid, ww, start=(/1,1,1, lst_tdx/)) )
                  endif
                  if (tt_varid /= 0) then
                     call check (nf90_put_var( ncid, tt_varid, tt, start=(/1,1,1, lst_tdx/)) )
                  endif
                  if (ipgt_varid /= 0) then
                     call check (nf90_put_var( ncid, ipgt_varid, ipgt, start=(/1,1, lst_tdx/)) )
                  endif
                  if (ust_varid /= 0) then
                     call check (nf90_put_var( ncid, ust_varid, ustar, start=(/1,1, lst_tdx/)) )
                  endif
                  if (zi_varid /= 0) then
                     call check (nf90_put_var( ncid, zi_varid, zi, start=(/1,1, lst_tdx/)) )
                  endif
                  if (el_varid /= 0) then
                     call check (nf90_put_var( ncid, el_varid, el, start=(/1,1, lst_tdx/)) )
                  endif
                  if (wst_varid /= 0) then
                     call check (nf90_put_var( ncid, wst_varid, wstar, start=(/1,1, lst_tdx/)) )
                  endif
                  if (rmm_varid /= 0) then
                     call check (nf90_put_var( ncid, rmm_varid, rmm, start=(/1,1, lst_tdx/)) )
                  endif
                  if (tempk_varid /= 0) then
                     call check (nf90_put_var( ncid, tempk_varid, tempk, start=(/1,1, lst_tdx/)) )
                  endif
                  if (rho_varid /= 0) then
                     call check (nf90_put_var( ncid, rho_varid, rho, start=(/1,1, lst_tdx/)) )
                  endif
                  if (qsw_varid /= 0) then
                     call check (nf90_put_var( ncid, qsw_varid, qsw, start=(/1,1, lst_tdx/)) )
                  endif
                  if (irh_varid /= 0) then
                     call check (nf90_put_var( ncid, irh_varid, irh, start=(/1,1, lst_tdx/)) )
                  endif
                  if (ipcode_varid /= 0) then
                     call check (nf90_put_var( ncid, ipcode_varid, ipcode, start=(/1,1, lst_tdx/)) )
                  endif
               endif
               lst_tstamp = tstamp
               lst_tdx = tdx
            endif

         elseif (clab(1:5) == 'V-LEV') then
            read(clab(6:8),*)k
            read(11)clab,idum(1:4),vv(:,:,k)
            !print*,clab, idum(1:2)
         elseif (clab(1:5) == 'WFACE') then
            read(clab(6:8),*)k
            read(11)clab,idum(1:4),ww(:,:,k)
            !print*,clab, idum(1:2)
         elseif (clab(1:5) == 'T-LEV') then
            read(clab(6:8),*)k
            read(11)clab,idum(1:4),tt(:,:,k)
         elseif (clab(1:5) == 'IPGT') then
            read(11)clab,idum(1:4),ipgt
         elseif (clab(1:5) == 'USTAR') then
            read(11)clab,idum(1:4),ustar
         elseif (clab(1:5) == 'ZI') then
            read(11)clab,idum(1:4),zi
         elseif (clab(1:5) == 'EL') then
            read(11)clab,idum(1:4),el
         elseif (clab(1:5) == 'WSTAR') then
            read(11)clab,idum(1:4),wstar
         elseif (clab(1:5) == 'RMM') then
            read(11)clab,idum(1:4),rmm
         elseif (clab(1:5) == 'TEMPK') then
            read(11)clab,idum(1:4),tempk
         elseif (clab(1:5) == 'RHO') then
            read(11)clab,idum(1:4),rho
         elseif (clab(1:5) == 'QSW') then
            read(11)clab,idum(1:4),qsw
         elseif (clab(1:5) == 'IRH') then
            read(11)clab,idum(1:4),irh
         elseif (clab(1:6) == 'IPCODE') then
            read(11)clab,idum(1:4),ipcode
         else
            print*,'sssss'
            stop
         endif
      enddo

      ! dump the last tstep
      if (lst_tstamp /= '') then
         print*, 'tim', tstamp
         call check (nf90_put_var( ncid, tim_varid, tstamp, start=(/1, lst_tdx/)))
         ! write values
         if (uu_varid /= 0) then
            call check (nf90_put_var( ncid, uu_varid, uu, start=(/1,1,1, lst_tdx/)) )
         endif
         if (vv_varid /= 0) then
            call check (nf90_put_var( ncid, vv_varid, vv, start=(/1,1,1, lst_tdx/)) )
         endif
         if (ww_varid /= 0) then
            call check (nf90_put_var( ncid, ww_varid, ww, start=(/1,1,1, lst_tdx/)) )
         endif
         if (tt_varid /= 0) then
            call check (nf90_put_var( ncid, tt_varid, tt, start=(/1,1,1, lst_tdx/)) )
         endif
         if (ipgt_varid /= 0) then
            call check (nf90_put_var( ncid, ipgt_varid, ipgt, start=(/1,1, lst_tdx/)) )
         endif
         if (ust_varid /= 0) then
            call check (nf90_put_var( ncid, ust_varid, ustar, start=(/1,1, lst_tdx/)) )
         endif
         if (zi_varid /= 0) then
            call check (nf90_put_var( ncid, zi_varid, zi, start=(/1,1, lst_tdx/)) )
         endif
         if (el_varid /= 0) then
            call check (nf90_put_var( ncid, el_varid, el, start=(/1,1, lst_tdx/)) )
         endif
         if (wst_varid /= 0) then
            call check (nf90_put_var( ncid, wst_varid, wstar, start=(/1,1, lst_tdx/)) )
         endif
         if (rmm_varid /= 0) then
            call check (nf90_put_var( ncid, rmm_varid, rmm, start=(/1,1, lst_tdx/)) )
         endif
         if (tempk_varid /= 0) then
            call check (nf90_put_var( ncid, tempk_varid, tempk, start=(/1,1, lst_tdx/)) )
         endif
         if (rho_varid /= 0) then
            call check (nf90_put_var( ncid, rho_varid, rho, start=(/1,1, lst_tdx/)) )
         endif
         if (qsw_varid /= 0) then
            call check (nf90_put_var( ncid, qsw_varid, qsw, start=(/1,1, lst_tdx/)) )
         endif
         if (irh_varid /= 0) then
            call check (nf90_put_var( ncid, irh_varid, irh, start=(/1,1, lst_tdx/)) )
         endif
         if (ipcode_varid /= 0) then
            call check (nf90_put_var( ncid, ipcode_varid, ipcode, start=(/1,1, lst_tdx/)) )
         endif
      endif



      
      call check( nf90_close( ncid))
   end subroutine

   subroutine cmtime_to_tstamp(idh, isec, tstamp)
      implicit none
      integer, intent(in) :: idh, isec
      character(len=19), intent(out) :: tstamp
      integer, dimension(13), parameter :: jdays_leap = (/ &
         0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 &
         /)
      integer, dimension(13), parameter :: jdays_regular = (/ &
         0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 &
         /)
      integer :: iyr, ijul, ihr,  imn, isc, imo, idy

      iyr = (idh / 100000)
      ijul = ( idh - iyr * 100000 ) / 100
      ihr = mod(idh , 100)
      imn = isec / 60
      isc = mod(isec , 60)

      if ( mod(iyr , 4) == 9) then
         imo=findloc(ijul < jdays_leap, .true., dim=1) - 1
         idy = ijul - jdays_leap(imo)
      else
         imo=findloc(ijul < jdays_regular, .true., dim=1) - 1
         idy = ijul - jdays_regular(imo)
      endif


      write(tstamp, "(i4, '-', i2.2, '-', i2.2,'_', i2.2, ':', i2.2, ':', i2.2 )") iyr, imo, idy, ihr, imn, isc

   end subroutine

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
              iwfcod,nssta,nusta,npsta,nowsta,nlu,iwat1,iwat2,lcalgrd,pmap,datum,daten,&
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
      character(len=80) :: fmtf, fmti, fmta
      read(11)dataset, dataver, datamod
      read(11)ncom
      do i=1,ncom
         read(11) comment
      enddo

      ! header, scalar part
      read(11)ibyr,ibmo,ibdy,ibhr,ibsec,ieyr,iemo,iedy,iehr,iesec,axtz,&
              irlg,irtype,nx,ny,nz,dgrid,xorigr,yorigr,&
              iwfcod,nssta,nusta,npsta,nowsta,nlu,iwat1,iwat2,lcalgrd,pmap,datum,daten,&
              feast,fnorth,utmhem,iutmzn,rnlat0,relon0,xlat1,xlat2

      fmtf = "(a8,',',2(i9,',',i4,','),3(i4,','),e14.7)"
      fmti = "(a8,',',2(i9,',',i4,','),3(i4,','),i14)"
      fmta = "(a8,',',2(a9,',',a4,','),3(a4,','),a14)"
      ! data
      !print*,ncheader
      if (.not.lnoisy) write(*,fmta) 'variable', 'bdh', 'bsec', 'edh', 'esec','i','j','k', 'value'
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
                  (idum(1) < idh1 .or. (idum(1) == idh1 .and. idum(2) <= isec1)) ) then
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
                  (idum(1) < idh1 .or. (idum(1) == idh1 .and. idum(2) <= isec1)) ) then
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
                  (idum(1) < idh1 .or. (idum(1) == idh1 .and. idum(2) <= isec1)) ) then
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
                  (idum(1) < idh1 .or. (idum(1) == idh1 .and. idum(2) <= isec1)) ) then
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
                  (idum(1) < idh1 .or. (idum(1) == idh1 .and. idum(2) <= isec1)) ) then
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
                  (idum(1) < idh1 .or. (idum(1) == idh1 .and. idum(2) <= isec1)) ) then
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
                  (idum(1) < idh1 .or. (idum(1) == idh1 .and. idum(2) <= isec1)) ) then
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
                  (idum(1) < idh1 .or. (idum(1) == idh1 .and. idum(2) <= isec1)) ) then
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
                  (idum(1) < idh1 .or. (idum(1) == idh1 .and. idum(2) <= isec1)) ) then
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
                  (idum(1) < idh1 .or. (idum(1) == idh1 .and. idum(2) <= isec1)) ) then
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
                  (idum(1) < idh1 .or. (idum(1) == idh1 .and. idum(2) <= isec1)) ) then
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
                  (idum(1) < idh1 .or. (idum(1) == idh1 .and. idum(2) <= isec1)) ) then
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
                  (idum(1) < idh1 .or. (idum(1) == idh1 .and. idum(2) <= isec1)) ) then
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
                  (idum(1) < idh1 .or. (idum(1) == idh1 .and. idum(2) <= isec1)) ) then
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
                  (idum(1) < idh1 .or. (idum(1) == idh1 .and. idum(2) <= isec1)) ) then
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
         ! update default for endtime
         idh1 = idh0
         isec1 = isec0
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
   subroutine check(status)
      use netcdf
      IMPLICIT NONE
      !---------------------------------------------------------------------
      integer(4), intent ( in) :: status
      !
      !- End of header -----------------------------------------------------
      !---------------------------------------------------------------------
      ! Subroutine Body
      !---------------------------------------------------------------------

      if(status /= nf90_noerr) then 
         print *, trim(nf90_strerror(status)) 
         stop "Stopped" 
      end if 
   end subroutine check


   Pure Function to_upper (str) Result (string) 
      !   ==============================
      !   Changes a string to upper case
      !   ============================== 
      Implicit None 
      Character(*), Intent(In) :: str 
      Character(LEN(str))      :: string
     
      Integer :: ic, i 
      
      Character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
      Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz' 
      
      !   Capitalize each letter if it is lowecase 
      string = str 
      do i = 1, LEN_TRIM(str) 
        ic = INDEX(low, str(i:i)) 
        if (ic > 0) string(i:i) = cap(ic:ic) 
     end do 
  End Function to_upper
   Pure Function to_lower (str) Result (string) 
      !   ==============================
      !   Changes a string to upper case
      !   ============================== 
      Implicit None 
      Character(*), Intent(In) :: str 
      Character(LEN(str))      :: string
     
      Integer :: ic, i 
      
      Character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
      Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz' 
      
      !   Capitalize each letter if it is lowecase 
      string = str 
      do i = 1, LEN_TRIM(str) 
        ic = INDEX(cap, str(i:i)) 
        if (ic > 0) string(i:i) = low(ic:ic) 
     end do 
  End Function to_lower

end program

! vim: et sw=3
