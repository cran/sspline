C***************************************************************************
C!   sspline.f
C***************************************************************************

      subroutine sspline(lon, lat, y, n, m, smth, lambda, gcv, varhat,
     $     coefc, coefd, yhat, Q, ldq, S, tau, work, lwork, d, e)


C
C     Purpose:
C     ========
C
C        The Fortran interface for solving spherical spline problem (can be
C        called by other Fortran programs/subroutines or R .Fortran 
C        interface.)
C
C     Arguments:
C     ==========
C
C        lon: input, double array, dimension n by 1, the longitudes.
C        lat: input, double array, dimension n by 1, the latitudes.
C          y: input, double array, dimension n by 1, the observations.
C          n: input, integer, the dimension of lon, lat and y.
C          m: input, integer, the order of smoothing on the spherical
C             function.
C       smth: input, integer, the way of choosing the smoothing parameter,
C                0: using the gcv return as the smoothing parameter
C                1: using the smoothing parameter specified in lambda.
C             more: to be implemented...
C     lambda: input/output, double, used to store the optimal lambda
C             returned by gcv or the prespecified lambda (if smth = 1).
C
C     varhat: output, double, the estimated variance.
C        gcv: output, double, the optimal gcv value.
C      coefc: output, double array, dimension n by 1, the coefficients 
C             in the estimated function.
C      coefd: output, double, the coefficients in the estimated
C             function.
C       yhat: ouput/workspace, double array, the predicted values at
C             the observed points.
C          Q: workspace, double array, dimension ldq by n, used for
C             storing with the kernel matrix Q.
C        ldq: input, integer, leading dimension of the working arrry Q.
C          S: workspace, double array, dimension n by 1.
C        tau: workspace, double array, dimension n by 1, used in the matrix
C             decomposition.
C      lwork: input, integer, the leading dimension of working array work.
C       work: workspace, double array, dimension lwork by 1, used in the
C             matrix decomposition.
C       d, e: workspace, double array, dimension n by 1, used in the matrix
C             decomposition.
C
C     ***Note***: all the longitudes and latitudes are in the radiant measure.
C


C
C     declaration of the arguments
C
      integer n, m, smth, ldq, lwork
      double precision lon(n), lat(n), y(n)
      double precision yhat(n), S(n), Q(ldq, n)
      double precision tau(n), work(lwork)
      double precision d(n), e(n), coefc(n)
      double precision coefd, lambda, varhat, gcv

C
C     declaration of the local variables
C
      double precision  low, up, prec, trace


C
C     generate the matrix Q from the longitudes and the latitudes.
C
      if (m .lt. 1 .or. m .gt. 10) then
         stop 'the order of smoothing should be between 1 and 10!'
      end if

      call dmatgen(lon, lat, n, m, Q)

C
C     save the y vector into yhat for later use
C
      call dcopy(n, y, 1, yhat, 1)

C
C     set up the problem: find the orthogonal matrix, eigenvalues and etc. 
C
      call dgcvst(Q, ldq, n, y, S, d, e, tau, work, lwork)

C
C********************* beginning of the gcv searching **********************
C
      if (smth .eq. 0) then
C     
C     set the range of searching for lambda
C
         prec = 1.0d0
         do 10 while (1.0d0 + prec > 1.0d0)
            prec = prec/2.0d0
 10      continue
      
         low = 2.0d0 * prec
         up  = 1.0d0 / low 
C
C     search for the lambda
C
         call gcvsch(y(2), d, n, low, up, lambda, gcv, varhat)

      else

         if (lambda < 0.0d0) then
            stop 'The smoothing parameter lambda must be positive!'
         end if

         call dgcvev(y(2), d, n, dlog(lambda), gcv)

C
C     calculate the estimated variance:
C
C        varhat = gcv * lambda * tr((T + n*lambda*I)^{-1})
C
         trace = 0.0d0

         do 20, i = 1, n-1
            trace = trace + lambda / (d(i) + dble(n)*lambda)
 20      continue

         varhat = gcv * trace
         
      end if
C
C************************ end of the gcv searching ************************
C

C
C     solve the equation using the lambda found in the above step or using
C     the prespecified lambda.
C
      call deqsv(Q, ldq, n, y, d, lambda, coefc, coefd, S, 
     $           work, lwork, e)

C
C     calculate the estimated y value: yhat = y - n*lambda*coefc
C
      call daxpy(n, -dble(n)*lambda, coefc, 1, yhat, 1)


      return
      end















