C******************************************************************************
C!    ssplfit.f
C******************************************************************************

      subroutine ssplfit(coefc, coefd, nobs, m, origlon, origlat, 
     $     	predlon, predlat, npred, vpred)


C
C     Purpose:
C     ========
C     
C        fit the values at {(predlon(i), predlat(i))}, 1 <= i <= npred by 
C        using the information from the subroutine sspline.
C
C
C     Arguments:
C     ==========
C     
C        coefc: input, double, dimension nobs by 1, the estimated 
C               coefficient c.
C        coefd: input, double, scalar, the estimated coefficient d.
C        origlon: input, double, dimension nobs by 1, the original longitudes.
C        origlat: input, double, dimension nobs by 1, the original latitudes.
C        nobs: input, integer, the dimension of vector coefc.
C	 m: input, integer, the order of smoothing used in fitting the sspline.
C        predlon: input, double, dimension npred by 1, the longitudes on which
C             the prediction are to be made.
C        predlat: input, double, dimension npred by 1, the latitudes on which
C             the prediction are to be made.
C        npred: input, integer, the dimension of vector lon and lat.
C
C        vpred: output, double, dimension npred by 1, the predicted values.
C
C
C     ***Note***: all the longitudes and latitudes are in the radiant measure.
C

      parameter (pi=3.1415926d0)

C
C     declaration of the arguments
C
      integer nobs, npred, m
      double precision coefd, coefc(nobs), origlon(nobs), origlat(nobs)
      double precision predlon(npred), predlat(npred), vpred(npred)

C
C     declaration of the local variables
C
      double precision zfunc, zval, qval, tmpr

      double precision A, C, W, q1, q2, q3, q4
      double precision q5, q6, q7, q8, q9, q10

      integer factorial, fac1, fac2


      if (m .lt. 1 .or. m .gt. 10) then
         stop 'the order of smoothing should be between 1 and 10!'
      end if
      
C
C     calculate the two factorials relatated to the order m
C
      fac1 = factorial(m)
      fac2 = factorial(m + 1)

C
C     calculating the predictted values
C

C
C     case 1: m = 1
C
      if (m .eq. 1) then
         do 20, i = 1, npred
            vpred(i) = coefd
            do 10, j = 1, nobs
               zval = zfunc(predlon(i), predlat(i), origlon(j), 
     $              origlat(j))

               if (dabs(zval - 1.0d0) .le. 1e-7) then
                  qval = 2.0d0 / DBLE(m + 2)
               else
                  W = (1.0d0 - zval) / 2.0d0
                  C = 2.0d0 * dsqrt(W)
                  A = dlog(1.0d0 + 1.0d0 / dsqrt(W))
                  qval = q1(A, C, W)
               end if

               tmpr = (qval/fac1 - 1.0d0/fac2) / (2.0d0 * pi)
               vpred(i) = vpred(i) + coefc(j) * tmpr
 10         continue
 20      continue
C
C     case 2: m = 2
C
      else if (m .eq. 2) then
         do 40, i = 1, npred
            vpred(i) = coefd
            do 30, j = 1, nobs
               zval = zfunc(predlon(i), predlat(i), origlon(j), 
     $              origlat(j))

               if (dabs(zval - 1.0d0) .le. 1e-7) then
                  qval = 2.0d0 / DBLE(m + 2)
               else
                  W = (1.0d0 - zval) / 2.0d0
                  C = 2.0d0 * dsqrt(W)
                  A = dlog(1.0d0 + 1.0d0 / dsqrt(W))
                  qval = q2(A, C, W)
               end if

               tmpr = (qval/fac1 - 1.0d0/fac2) / (2.0d0 * pi)
               vpred(i) = vpred(i) + coefc(j) * tmpr
 30         continue
 40      continue
C
C     case 3: m = 3
C
      else if (m .eq. 3) then
         do 60, i = 1, npred
            vpred(i) = coefd
            do 50, j = 1, nobs
               zval = zfunc(predlon(i), predlat(i), origlon(j), 
     $              origlat(j))

               if (dabs(zval - 1.0d0) .le. 1e-7) then
                  qval = 2.0d0 / DBLE(m + 2)
               else
                  W = (1.0d0 - zval) / 2.0d0
                  C = 2.0d0 * dsqrt(W)
                  A = dlog(1.0d0 + 1.0d0 / dsqrt(W))
                  qval = q3(A, C, W)
               end if

               tmpr = (qval/fac1 - 1.0d0/fac2) / (2.0d0 * pi)
               vpred(i) = vpred(i) + coefc(j) * tmpr
 50         continue
 60      continue
C
C     case 4: m = 4
C
      else if (m .eq. 4) then
         do 80, i = 1, npred
            vpred(i) = coefd
            do 70, j = 1, nobs
               zval = zfunc(predlon(i), predlat(i), origlon(j), 
     $              origlat(j))

               if (dabs(zval - 1.0d0) .le. 1e-7) then
                  qval = 2.0d0 / DBLE(m + 2)
               else
                  W = (1.0d0 - zval) / 2.0d0
                  C = 2.0d0 * dsqrt(W)
                  A = dlog(1.0d0 + 1.0d0 / dsqrt(W))
                  qval = q4(A, C, W)
               end if

               tmpr = (qval/fac1 - 1.0d0/fac2) / (2.0d0 * pi)
               vpred(i) = vpred(i) + coefc(j) * tmpr
 70         continue
 80      continue
C
C     case 5: m = 5
C
      else if (m .eq. 5) then
         do 100, i = 1, npred
            vpred(i) = coefd
            do 90, j = 1, nobs
               zval = zfunc(predlon(i), predlat(i), origlon(j), 
     $              origlat(j))

               if (dabs(zval - 1.0d0) .le. 1e-7) then
                  qval = 2.0d0 / DBLE(m + 2)
               else
                  W = (1.0d0 - zval) / 2.0d0
                  C = 2.0d0 * dsqrt(W)
                  A = dlog(1.0d0 + 1.0d0 / dsqrt(W))
                  qval = q5(A, C, W)
               end if

               tmpr = (qval/fac1 - 1.0d0/fac2) / (2.0d0 * pi)
               vpred(i) = vpred(i) + coefc(j) * tmpr
 90         continue
 100     continue
C
C     case 6: m = 6
C
      else if (m .eq. 6) then
         do 120, i = 1, npred
            vpred(i) = coefd
            do 110, j = 1, nobs
               zval = zfunc(predlon(i), predlat(i), origlon(j), 
     $              origlat(j))

               if (dabs(zval - 1.0d0) .le. 1e-7) then
                  qval = 2.0d0 / DBLE(m + 2)
               else
                  W = (1.0d0 - zval) / 2.0d0
                  C = 2.0d0 * dsqrt(W)
                  A = dlog(1.0d0 + 1.0d0 / dsqrt(W))
                  qval = q6(A, C, W)
               end if

               tmpr = (qval/fac1 - 1.0d0/fac2) / (2.0d0 * pi)
               vpred(i) = vpred(i) + coefc(j) * tmpr
 110        continue
 120     continue
C
C     case 7: m = 7
C
      else if (m .eq. 7) then
         do 140, i = 1, npred
            vpred(i) = coefd
            do 130, j = 1, nobs
               zval = zfunc(predlon(i), predlat(i), origlon(j), 
     $              origlat(j))

               if (dabs(zval - 1.0d0) .le. 1e-7) then
                  qval = 2.0d0 / DBLE(m + 2)
               else
                  W = (1.0d0 - zval) / 2.0d0
                  C = 2.0d0 * dsqrt(W)
                  A = dlog(1.0d0 + 1.0d0 / dsqrt(W))
                  qval = q7(A, C, W)
               end if

               tmpr = (qval/fac1 - 1.0d0/fac2) / (2.0d0 * pi)
               vpred(i) = vpred(i) + coefc(j) * tmpr
 130        continue
 140     continue
C
C     case 8: m = 8
C
      else if (m .eq. 8) then
         do 160, i = 1, npred
            vpred(i) = coefd
            do 150, j = 1, nobs
               zval = zfunc(predlon(i), predlat(i), origlon(j), 
     $              origlat(j))

               if (dabs(zval - 1.0d0) .le. 1e-7) then
                  qval = 2.0d0 / DBLE(m + 2)
               else
                  W = (1.0d0 - zval) / 2.0d0
                  C = 2.0d0 * dsqrt(W)
                  A = dlog(1.0d0 + 1.0d0 / dsqrt(W))
                  qval = q8(A, C, W)
               end if

               tmpr = (qval/fac1 - 1.0d0/fac2) / (2.0d0 * pi)
               vpred(i) = vpred(i) + coefc(j) * tmpr
 150        continue
 160     continue
C
C     case 9: m = 9
C
      else if (m .eq. 9) then
         do 180, i = 1, npred
            vpred(i) = coefd
            do 170, j = 1, nobs
               zval = zfunc(predlon(i), predlat(i), origlon(j), 
     $           origlat(j))

               if (dabs(zval - 1.0d0) .le. 1e-7) then
                  qval = 2.0d0 / DBLE(m + 2)
               else
                  W = (1.0d0 - zval) / 2.0d0
                  C = 2.0d0 * dsqrt(W)
                  A = dlog(1.0d0 + 1.0d0 / dsqrt(W))
                  qval = q9(A, C, W)
               end if

               tmpr = (qval/fac1 - 1.0d0/fac2) / (2.0d0 * pi)
               vpred(i) = vpred(i) + coefc(j) * tmpr
 170        continue
 180     continue
C
C     case 10: m = 10
C
      else if (m .eq. 10) then
         do 200, i = 1, npred
            vpred(i) = coefd
            do 190, j = 1, nobs
               zval = zfunc(predlon(i), predlat(i), origlon(j), 
     $              origlat(j))

               if (dabs(zval - 1.0d0) .le. 1e-7) then
                  qval = 2.0d0 / DBLE(m + 2)
               else
                  W = (1.0d0 - zval) / 2.0d0
                  C = 2.0d0 * dsqrt(W)
                  A = dlog(1.0d0 + 1.0d0 / dsqrt(W))
                  qval = q10(A, C, W)
               end if

               tmpr = (qval/fac1 - 1.0d0/fac2) / (2.0d0 * pi)
               vpred(i) = vpred(i) + coefc(j) * tmpr
 190        continue
 200     continue
      end if

      return
      end









