      character*500 InFNam, OutFnam
      real*4        mjd, pa, avg, sigma, AsceMin, AsceMax,
     +              DescMin, DescMax
      real*8        dMJD, nullval, sum0, sum1, sumsq0, sumsq1
      integer*4     nargs, iargc, nOut, FileID, status,
     +              readwrite, blocksize, nrows, ncols, nPA, nAsce,
     +              nMJD, nIncd, hdutype, k, felem, nelems, stat(4)
      integer       n, n0, n1
      byte          incd, asce
      logical*4     dbg
      logical       anynull
c      
      data nOut/0/, dbg/.false./, n0,n1/2*0/,
     +     sum0,sum1,sumsq0,sumsq1/4*0.0d0/,
     +     AsceMin,DescMin/2*9.9e9/, AsceMax,DescMax/2*-9.9e9/
c     
c-----------------------------------------------------------------------
c
      nargs = iargc()
      if (nargs .lt. 2) then
        print *,'fbt vsn 1.1  B71002'
        print *,'usage: fbt infile outfile'
        print *
        print *,
     +    'where: infile  is an unWISE NEO3 FITS binary table file'
        print *,'       outfile is the ASCII-text file for PSFavr8'
        stop
      end if
      call getarg(1,InFNam)
      if (Access(InFNam(1:lnblnk(InFNam)),' ') .ne. 0) then
        print *,'File not found: ',InFNam(1:lnblnk(InFNam))
        call exit(64)
      end if
      call getarg(2,OutFnam)
      dbg = nargs .gt. 2
c-----------------------------------------------------------------------
      status = 0
      call ftgiou(FileID,status)         ! get unit number
      readwrite = 0                      ! open FITS file
      call ftopen(FileID,InFNam,readwrite,blocksize,status)
      if (status .ne. 0) then
        write(6,'(a)') 'ERROR: Could not open '//trim(InFNam)
        call exit(64)
      endif
      call ftmahd(FileID,2,hdutype,status) ! move to header #2
      if (dbg) print *,'hdutype for nhdu-2, status:     ',
     +                  hdutype, status
      status = 0
      call ftgnrw(FileID,nrows,status)     ! get #rows
      if (status .ne. 0) go to 3001
      if (dbg) print *,'No. of rows in table and status:', nrows, status
      call ftgncl(FileID,ncols,status)     ! get #cols
      if (status .ne. 0) go to 3002
      if (dbg) print *,'No. of cols in table and status:', ncols, status
      call ftgcno(FileID, .false., 'pa      ', nPA,  status)    ! row# for PA
      if (status .ne. 0) go to 3003
      if (dbg) print *,'Col. no and status for pa:      ', nPA, status
      call ftgcno(FileID, .false., 'ascending', nAsce, status)  ! row# for ascending
      if (dbg) print *,'Col. no and status for ascending:     ', nAsce, status
      if (status .ne. 0) go to 3004
      call ftgcno(FileID, .false., 'mjd     ', nMJD, status)  ! row# for mjd
      if (dbg) print *,'Col. no and status for mjd:     ', nMJD, status
      if (status .ne. 0) go to 3005
      call ftgcno(FileID, .false., 'included', nIncd, status) ! row# for included
      if (status .ne. 0) go to 3006
      if (dbg) print *,'Col. no and status for included:', nIncd, status
c      
      open (12, file = OutFnam)
      felem  = 1
      nelems = 1
c      
      do 100 n = 1, nrows
c      
        k = 1
        call ftgcve(FileID,nPA,n,felem,nelems,0.0d0,pa,anynull,status)
        if (status .ne. 0) go to 3000
        stat(1) = status
        k = 2
        call ftgcvb(FileID,nAsce,n,felem,nelems,0.0d0,asce,anynull,status)
        if (status .ne. 0) go to 3000
        stat(2) = status
        k = 3
        call ftgcvd(FileID,nMJD,n,felem,nelems,0.0d0,dMJD,anynull,status)
        if (status .ne. 0) go to 3000
        mjd = dMJD
        stat(3) = status
        k = 4
        call ftgcvb(FileID,nIncd,n,felem,nelems,0,incd,anynull,status)
        if (status .ne. 0) go to 3000
        stat(4) = status
c        
        if (dbg .and. (n .eq. 1)) then
          print *,'n,felem,pa,status:        ', n, felem, pa,   stat(1)
          print *,'n,felem,ascending,status: ', n, felem, asce, stat(2)
          print *,'n,felem,mjd,status:       ', n, felem, dMJD, stat(3)
          print *,'n,felem,incd,status:      ', n, felem, incd, stat(4)
        end if
c
        if (incd .eq. 0) go to 100       
        write (12, '(I1,2F13.2)') asce, pa, mjd
        nOut = nOut + 1
        if (asce .eq. 0) then
          n0 = n0 + 1
          sum0   = sum0   + pa
          sumsq0 = sumsq0 + pa**2
          if (pa .lt. DescMin) DescMin = pa
          if (pa .gt. DescMax) DescMax = pa
        else
          n1 = n1 + 1
          sum1   = sum1   + pa
          sumsq1 = sumsq1 + pa**2       
          if (pa .lt. AsceMin) AsceMin = pa
          if (pa .gt. AsceMax) AsceMax = pa
        end if
100   continue      
c     
      print *,'No. of rows read:     ', nrows
      print *,'No. of lines written: ', nOut
      call ftclos(FileID, status)
      call ftfiou(FileID, status)
c
      print *,'No. ascending:',n1,'; no. descending:',n0
      if (n1 .gt. 0) then
        avg   = sum1/dfloat(n1)
        sigma = dsqrt(dabs(sumsq1/dfloat(n1) - dble(avg)**2))
        print *,'Mean and sigma of PA ascending: ', avg, sigma
        print *,'Min and Max ascending PA:       ', AsceMin, AsceMax
        if (sigma .gt. 5.0) then
          print *,'WARNING: sigma is excessively high!'        
          print *,'         some epochs may have both scan directions'
        end if
      end if
      if (n0 .gt. 0) then
        avg   = sum0/dfloat(n0)
        sigma = dsqrt(dabs(sumsq0/dfloat(n0) - dble(avg)**2))
        print *,'Mean and sigma of PA descending:', avg, sigma
        print *,'Min and Max descending PA:      ', DescMin, DescMax
        if (sigma .gt. 5.0) then
          print *,'WARNING: sigma is excessively high!'        
          print *,'         some epochs may have both scan directions'
        end if
      end if
c
      stop
c
3000  print *,'ERROR reading table row no.', n,
     +        '; parameter', k, '; status =',status
      call exit(64)
3001  print *,'ERROR: unable to get @rows; status =',status
      call exit(64)
3002  print *,'ERROR: unable to get #cols; status =',status
      call exit(64)
3003  print *,'ERROR: pa not found; status =',status
      call exit(64)
3004  print *,'ERROR: ascending not found; status =',status
      call exit(64)
3005  print *,'ERROR: mjd not found; status =',status
      call exit(64)
3006  print *,'ERROR: included not found; status =',status
      call exit(64)
c      
      end
c      
c=======================================================================
c
      subroutine Cel2Ec(RA, Dec, Long, Lat)
c
      real*4 RA, Dec, Long, Lat, SOb, Cob, X, Y, Z, d2r,
     +       cRA, cDec, sRA, sDec, X2, Y2
c
c   Obliquity(2015) in J2000: 23.43734105     
c
      data d2r/1.745329252e-2/, cOb, sOb/0.9174956, 0.39777459/
c
c-----------------------------------------------------------------------
c
      cRA   = cos(d2r*RA)
      cDec  = cos(d2r*Dec)
      sRA   = sin(d2r*RA)
      sDec  = sin(d2r*Dec)
c
      X =  sDec
      Y = -cDec*sRA
      Z =  cDec*cRA
c
      X2 =  X*Cob + Y*Sob
      Y2 = -X*Sob + Y*Cob
c     Z2 =  Z
c
      Lat  = asin(X2)/d2r
      Long = atan2(-Y2,Z)/d2r
      if (Long .lt. 0.0) Long = Long + 360.0
c
      return
      end
