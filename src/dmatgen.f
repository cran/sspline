C*************************************************************************
C!    dmatgen.f
C*************************************************************************

      double precision function zfunc(lon1, lat1, lon2, lat2) 

C
C     Purpose:
C     ========
C
C         define function to calculate the cosine value of two points
C         P1(lon1, lat1), P2(lon2, lat2) on the sphere.
C
C
C     Arguments:
C     ==========
C
C         lon1: double, longitude of the first point
C         lat1: double, latitude of the first point
C         lon2: double, longitude of the second point
C         lat2: double, latitude of the second point
C
C
C     Value:
C     ======
C
C        zfunc: double, the cosine value.
C
C
C     ***Note***: the longitude and latitude are in the radiant measure.
C

C
C     declaration of the arguments
C
      double precision lon1, lat1, lon2, lat2
      
C
C     declaration of the local variables
C
      double precision cosval(4), sinval(4)

C
C     calculate the temporary variables
C
      cosval(1) = dcos(lon1)
      cosval(2) = dcos(lat1)
      cosval(3) = dcos(lon2)
      cosval(4) = dcos(lat2)

      sinval(1) = dsin(lon1)
      sinval(2) = dsin(lat1)
      sinval(3) = dsin(lon2)
      sinval(4) = dsin(lat2)

C
C     calcualte the return value of the z function:
C     
C        zfunc = cos(lat1)*cos(lon1)*cos(lat2)*cos(lon2) +
C                cos(lat1)*sin(lon1)*cos(lat2)*sin(lon2) +
C                sin(lat1)*sin(lat2)
C
C        lon ~ longitude, lat ~ latitude 
C
      zfunc = cosval(2)*cosval(1)*cosval(4)*cosval(3) +
     $        cosval(2)*sinval(1)*cosval(4)*sinval(3) +
     $        sinval(2)*sinval(4)


      end

C###########################################################################

      integer function factorial(n)

C
C     Purpose:
C     ========
C
C        caculate the factorial n!
C
C     Arguments:
C        n: input, integer
C     
C        factorial: ouput, integer
C

      integer n

      factorial = 0

      if (n .lt. 0) then
         print *, 'Error occurs when calculating factorial!'
      else
         factorial = 1
         if (n .gt. 2) then
            do 10, i = 2, n
               factorial = i * factorial
 10         continue
         end if
      end if

      end
      
C########################################################################

      subroutine dmatgen(lon, lat, n, m, R)

C
C     Purpose:
C     ========
C
C        define the subroutine to generate the matrix R, where 
C        R(i, j) = R(P_i, P_j), and P_i, P_j are two points on 
C        the sphere. 
C        
C
C     Arguments:
C     ==========
C
C          n: input, integer, dimension of vector lon, lat
C          m: input, integer, order of smoothing
C
C        lon: input, double array, dimension n by 1, longitute
C        lat: input, double array, dimension n by 1, latitude
C
C        R: output, double array, returned values in R(1:n, 1:n)
C
C
C     ***Note***: the longitutes and latitudes are in the radiant measure.
C

      parameter (pi=3.1415926d0)

C
C     declaration of the arguments
C
      integer m, n
      double precision lon(n), lat(n), R(n, n)
      double precision zfunc, zval, qval

      double precision A, C, W, q1, q2, q3, q4
      double precision q5, q6, q7, q8, q9, q10

      integer factorial, fac1, fac2

      fac1 = factorial(m)
      fac2 = factorial(m + 1)

C
C     calculate the matrix R
C

C
C     case 1: m = 1
C
      if (m .eq. 1) then
         do 20, i = 1, n
            do 10, j = 1, i
               zval = zfunc(lon(i), lat(i), lon(j), lat(j))

               if (dabs(zval - 1.0d0) .le. 1e-7) then
                  qval = 2.0d0 / DBLE(m + 2)
               else
                  W = (1.0d0 - zval) / 2.0d0
                  C = 2.0d0 * dsqrt(W)
                  A = dlog(1.0d0 + 1.0d0 / dsqrt(W))
                  qval = q1(A, C, W)
               end if

               R(i, j) = (qval/fac1 - 1.0d0/fac2) / (2.0d0 * pi)
               R(j, i) = R(i, j)
 10         continue
 20      continue
C
C     case 2: m = 2
C
      else if (m .eq. 2) then
         do 40, i = 1, n
            do 30, j = 1, i
               zval = zfunc(lon(i), lat(i), lon(j), lat(j))

               if (dabs(zval - 1.0d0) .le. 1e-7) then
                  qval = 2.0d0 / DBLE(m + 2)
               else
                  W = (1.0d0 - zval) / 2.0d0
                  C = 2.0d0 * dsqrt(W)
                  A = dlog(1.0d0 + 1.0d0 / dsqrt(W))
                  qval = q2(A, C, W)
               end if

               R(i, j) = (qval/fac1 - 1.0d0/fac2) / (2.0d0 * pi)
               R(j, i) = R(i, j)
 30         continue
 40      continue
C
C     case 3: m = 3
C
      else if (m .eq. 3) then
         do 60, i = 1, n
            do 50, j = 1, i
               zval = zfunc(lon(i), lat(i), lon(j), lat(j))

               if (dabs(zval - 1.0d0) .le. 1e-7) then
                  qval = 2.0d0 / DBLE(m + 2)
               else
                  W = (1.0d0 - zval) / 2.0d0
                  C = 2.0d0 * dsqrt(W)
                  A = dlog(1.0d0 + 1.0d0 / dsqrt(W))
                  qval = q3(A, C, W)
               end if

               R(i, j) = (qval/fac1 - 1.0d0/fac2) / (2.0d0 * pi)
               R(j, i) = R(i, j)
 50         continue
 60      continue
C
C     case 4: m = 4
C
      else if (m .eq. 4) then
         do 80, i = 1, n
            do 70, j = 1, i
               zval = zfunc(lon(i), lat(i), lon(j), lat(j))

               if (dabs(zval - 1.0d0) .le. 1e-7) then
                  qval = 2.0d0 / DBLE(m + 2)
               else
                  W = (1.0d0 - zval) / 2.0d0
                  C = 2.0d0 * dsqrt(W)
                  A = dlog(1.0d0 + 1.0d0 / dsqrt(W))
                  qval = q4(A, C, W)
               end if

               R(i, j) = (qval/fac1 - 1.0d0/fac2) / (2.0d0 * pi)
               R(j, i) = R(i, j)
 70         continue
 80      continue
C
C     case 5: m = 5
C
      else if (m .eq. 5) then
         do 100, i = 1, n
            do 90, j = 1, i
               zval = zfunc(lon(i), lat(i), lon(j), lat(j))

               if (dabs(zval - 1.0d0) .le. 1e-7) then
                  qval = 2.0d0 / DBLE(m + 2)
               else
                  W = (1.0d0 - zval) / 2.0d0
                  C = 2.0d0 * dsqrt(W)
                  A = dlog(1.0d0 + 1.0d0 / dsqrt(W))
                  qval = q5(A, C, W)
               end if

               R(i, j) = (qval/fac1 - 1.0d0/fac2) / (2.0d0 * pi)
               R(j, i) = R(i, j)
 90         continue
 100     continue
C
C     case 6: m = 6
C
      else if (m .eq. 6) then
         do 120, i = 1, n
            do 110, j = 1, i
               zval = zfunc(lon(i), lat(i), lon(j), lat(j))

               if (dabs(zval - 1.0d0) .le. 1e-7) then
                  qval = 2.0d0 / DBLE(m + 2)
               else
                  W = (1.0d0 - zval) / 2.0d0
                  C = 2.0d0 * dsqrt(W)
                  A = dlog(1.0d0 + 1.0d0 / dsqrt(W))
                  qval = q6(A, C, W)
               end if

               R(i, j) = (qval/fac1 - 1.0d0/fac2) / (2.0d0 * pi)
               R(j, i) = R(i, j)
 110        continue
 120     continue
C
C     case 7: m = 7
C
      else if (m .eq. 7) then
         do 140, i = 1, n
            do 130, j = 1, i
               zval = zfunc(lon(i), lat(i), lon(j), lat(j))

               if (dabs(zval - 1.0d0) .le. 1e-7) then
                  qval = 2.0d0 / DBLE(m + 2)
               else
                  W = (1.0d0 - zval) / 2.0d0
                  C = 2.0d0 * dsqrt(W)
                  A = dlog(1.0d0 + 1.0d0 / dsqrt(W))
                  qval = q7(A, C, W)
               end if

               R(i, j) = (qval/fac1 - 1.0d0/fac2) / (2.0d0 * pi)
               R(j, i) = R(i, j)
 130        continue
 140     continue
C
C     case 8: m = 8
C
      else if (m .eq. 8) then
         do 160, i = 1, n
            do 150, j = 1, i
               zval = zfunc(lon(i), lat(i), lon(j), lat(j))

               if (dabs(zval - 1.0d0) .le. 1e-7) then
                  qval = 2.0d0 / DBLE(m + 2)
               else
                  W = (1.0d0 - zval) / 2.0d0
                  C = 2.0d0 * dsqrt(W)
                  A = dlog(1.0d0 + 1.0d0 / dsqrt(W))
                  qval = q8(A, C, W)
               end if

               R(i, j) = (qval/fac1 - 1.0d0/fac2) / (2.0d0 * pi)
               R(j, i) = R(i, j)
 150        continue
 160     continue
C
C     case 9: m = 9
C
      else if (m .eq. 9) then
         do 180, i = 1, n
            do 170, j = 1, i
               zval = zfunc(lon(i), lat(i), lon(j), lat(j))

               if (dabs(zval - 1.0d0) .le. 1e-7) then
                  qval = 2.0d0 / DBLE(m + 2)
               else
                  W = (1.0d0 - zval) / 2.0d0
                  C = 2.0d0 * dsqrt(W)
                  A = dlog(1.0d0 + 1.0d0 / dsqrt(W))
                  qval = q9(A, C, W)
               end if

               R(i, j) = (qval/fac1 - 1.0d0/fac2) / (2.0d0 * pi)
               R(j, i) = R(i, j)
 170        continue
 180     continue
C
C     case 10: m = 10
C
      else if (m .eq. 10) then
         do 200, i = 1, n
            do 190, j = 1, i
               zval = zfunc(lon(i), lat(i), lon(j), lat(j))

               if (dabs(zval - 1.0d0) .le. 1e-7) then
                  qval = 2.0d0 / DBLE(m + 2)
               else
                  W = (1.0d0 - zval) / 2.0d0
                  C = 2.0d0 * dsqrt(W)
                  A = dlog(1.0d0 + 1.0d0 / dsqrt(W))
                  qval = q10(A, C, W)
               end if

               R(i, j) = (qval/fac1 - 1.0d0/fac2) / (2.0d0 * pi)
               R(j, i) = R(i, j)
 190        continue
 200     continue
      end if

      return
      end

