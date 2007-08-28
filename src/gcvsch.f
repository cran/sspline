C***************************************************************************
C!   gcvsch.f
C***************************************************************************


      subroutine gcvsch(x, eigen, n, low, up, lambda, gcv, varhat)


C
C     Purpose:
C     ========
C     
C       seaching for the lambda minimizing the gcv score V(lambda) using
C       the gold section method.
C
C                        1/n * x'*diag(eigen + n*lambda, n)^(-2)*x
C          V(lambda) = ----------------------------------------------
C                       [1/n * tr(diag(eigen + n*lambda, n)^(-1))]^2
C
C     Arguments:
C     ==========
C     
C       x: input, double array, length n - 1, on entry, x = U'*F2'*y,
C          where y is the observation vector, U is the orthogonal
C          matrix of the eigenvalue-eignevector decompostion of 
C          F2'*Q*F2 (nonegative definite).
C       eigen: input, double array, length n - 1, on entry, the 
C              eigenvalues of the matrix F2'*Q*F2.
C       n: input, integer, (dimension of vector x and eigen) + 1
C
C       low: input, double, the lower bound for searching gcv lambda
C       up:  input, double, the upper bound for searching gcv lambda
C       lam: output, double, the optimal gcv lambda
C       gcv: output, double, the gcv score at the optimal lambda
C       varhat: output, double, the estimated variance
C

C
C     declaration of the arguments
C
      integer n
      double precision low, up, lambda, gcv 
      double precision x(n-1), eigen(n-1)
      double precision varhat

C
C     declaration of the local variables
C
      double precision  olow, oup, mlow, mup
      double precision  lgcv, rgcv, trace


C
C     convert the range of searching to log scale
C
      olow = dlog(low)
      oup  = dlog(up)
      mlow = oup  - (dsqrt(5.0d0) - 1.0d0)/2.0d0 * (oup - olow)
      mup  = olow + (dsqrt(5.0d0) - 1.0d0)/2.0d0 * (oup - olow)

C
C     search for the lambda, first calculate the gcv scores at the
C     two inner end points.
C
      call dgcvev(x, eigen, n, mlow, lgcv)
      call dgcvev(x, eigen, n, mup,  rgcv)

C--------------------(main loop starts here) -------------------------     

      do 10 while (dabs(mup - mlow) .ge. 1.0d-7)
         if (lgcv .lt. rgcv) then
            oup  = mup
            mup  = mlow
            rgcv = lgcv
            mlow = oup  - (dsqrt(5.0d0) - 1.0d0)/2.0d0 * (oup - olow)
            call dgcvev(x, eigen, n, mlow, lgcv)
         else
            olow = mlow
            mlow = mup
            lgcv = rgcv
            mup  = olow + (dsqrt(5.0d0) - 1.0d0)/2.0d0 * (oup - olow)
            call dgcvev(x, eigen, n, mup, rgcv)
         endif
 10   continue

C-----------------------(main loop ends)-----------------------------

C
C     calculate the optimal lambda and its gcv score.
C
      lambda = dexp((mup + mlow)/2.0d0)
      call dgcvev(x, eigen, n, (mup + mlow)/2.0d0, gcv)
      
C
C     calculate the estimated variance:
C
C        varhat = gcv * lambda * tr((T + n*lambda*I)^{-1})
C
      trace = 0.0d0

      do 20, i = 1, n-1
         trace = trace + lambda / (eigen(i) + dble(n)*lambda)
 20   continue

      varhat = gcv * trace


      return
      end

      


