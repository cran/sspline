C**************************************************************************
C!    dgcvev.f
C**************************************************************************


      subroutine dgcvev(z, eigen, n, loglam, gcv)


C
C     Purpose:
C     ========
C     
C       calculate the gcv value at lambda = exp(loglam), the formula is
C
C                       1/n * z'*diag(eigen + n*lambda, n)^(-2)*z
C          V(lambda) = ----------------------------------------------
C                       [1/n * tr(diag(eigen + n*lambda, n)^(-1))]^2
C
C     Arguments:
C     ==========
C     
C       z: input, double array, dimension (n - 1) by 1
C       eigen: input, double array, dimension (n - 1) by 1
C       n: input, integer, dimension of z (and eigen) + 1
C       loglam: input, double, logarithm of lambda (base e)
C
C       gcv: output, double, the gcv value at exp(loglam)
C


C
C     declaration of the arguments
C
      integer  n
      double precision  z(n-1), eigen(n-1)
      double precision  loglam, gcv

C     
C     declaration of the local variables
C
      double precision  lambda, numer, denom


      lambda = dexp(loglam)

      numer = 0.0d0
      do 10, i = 1, n-1
         numer = numer + (z(i)/(eigen(i) + n*lambda))**2
 10   continue

      denom = 0.0d0
      do 20, i = 1, n-1
         denom = denom + 1.0d0/(eigen(i) + n*lambda)
 20   continue
      denom = denom**2

      gcv = dble(n) * numer/denom


      return
      end







