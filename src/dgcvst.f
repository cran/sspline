C*************************************************************************
C!   dgcvst.f								     
C*************************************************************************


      subroutine dgcvst(Q, ldq, n, y, S, d, e, tau, work, lwork)


C
C     Purpose:
C     ========
C
C       compute the gcv function for the following problem
C
C           (Q + n*lambda*I)*c + S*d = y
C            S' * c = 0
C
C           where S = (1, 1, ..., 1)'
C
C       the gcv function is
C     
C                        n * (z'*(T + n*lambda*I)^(-2)*z)
C           V(lambda) = ----------------------------------
C                        [trace((T + n*lambda*I)^(-1))]^2
C
C
C     Arguments:
C     ==========
C
C        Q: input/output, double array, dimension ldq by n
C           on entry: the kernel matrix at points t1,t2,...,tn
C           on exit: Q(1:n-1, 1) = F2'*Q*F1, Q(1:n-1, 2:n) = U,
C           where F1, F2, U are the orthonormal matrices for the 
C           Q-R decomposition of all 1 vector S, and the orthogonal 
C           matrix for reducing the symmetric matrix F2'*Q*F2 to 
C           diagonal matrix respectively.
C        ldq: input, integer, leading dimension of Q
C        n: input, integer, number of observations.
C
C        y: input/output, double array, dimension n by 1, 
C           on entry: the observation vector
C           on exit: y(1) = F1'*y; y(2:n) = (U' * F2' * y), where
C           U and F2 are the orthogonal and orthonormal matrices
C           specified above, S = (F1 F2)*(R 0)', F2'*Q*F2 = U*D*U'.
C
C
C        d: output/workspace, double array, dimension (n - 1) by 1
C           on exit: eigenvalues of the symmetric matrix F2'*Q*F2 
C        e: output/workspace, double array, dimension (n - 2) by 1
C           on exit: off diagonal elements of the tridiagonal matrix
C           corresponding to the symmetric matrix (F2' * Q * F2)
C        tau: output/workspace, double array, dimension (n - 2) by 1
C             on exit: coefficients of elementary reflectors in the
C             reduction of (F2' * Q * F2) to tridiagonal matrix.
C
C
C        S: workspace, double array, dimension n by 1, on exit, S(1)
C           is stored with the R in the Q-R decomposition of S.
C        work: workspace, double array, dimension lwork by 1
C        lwork: input, integer, the dimension of work
C


C
C     declaration of the arguments
C
      integer  n, ldq, lwork
      double precision  d(n-1), e(n-2), S(n), y(n), tau(*) 
      double precision  work(lwork), Q(ldq, n)

C
C     declaration of the local variables
C
      integer  info, jpvt(1)

      jpvt(1) = 0

C
C     construct an all 1 vector S (length n)
C
      do 10, i = 1, n
         S(i) = 1
 10   continue

C
C     do the Q-R decomposition on the all 1 vector S (n by 1)
C
      call dgeqp3(n, 1, S, n, jpvt, tau, work, lwork, info)


C
C     get the matrix-vector product ((F1 F2)' * y) in y, where
C     (F1 F2) is the orthogonal matrix in the Q-R decomposition
C     of the n by 1 matrix S
C
      call dormqr('L', 'T', n, 1, 1, S, n, tau, y, n,
     $            work, lwork, info)


C
C     get the matrix-matrix product (F1 F2)' * Q * (F1 F2) in Q
C
      call dormqr('L', 'T', n, n, 1, S, n, tau, Q, ldq,
     $            work, lwork, info)

      call dormqr('R', 'N', n, n, 1, S, n, tau, Q, ldq,
     $            work, lwork, info)


C
C     move the elments of Q one row up to get (F2' * Q * F1) and
C     (F2' * Q * F2) in the first (n - 1) rows of Q
C
      do 30, i = 1, n - 1
         do 20, j = 1, n
            Q(i, j) = Q(i+1, j)
 20      continue
 30   continue


C      
C     reduce the symmetric matrix (F2' * Q * F2) to tridiagonol form
C
      call dsytrd('U', n-1, Q(1, 2), ldq, d, e, tau, work, lwork, info)


C
C     generate the orthogonal matrix in the above reduction for 
C     the next step
C
      call dorgtr('U', n-1, Q(1, 2), ldq, tau, work, lwork, info)

C
C     get the eigenvalues and eigenvectors of the original symmetric
C     matrix (d is overwritten with the eigenvalues, Q is overwritten
C     with the orthogonal matrix.)
C
      call dsteqr('V', n-1, d, e, Q(1, 2), ldq, work, info)


C
C     if the eigenvalue is too small (less than 1e-10), truncate it to 0.
C
      do 40, i = 1, n-1
         if (d(i) .le. 1.0d-10) then
            d(i) = 0.0d0
         endif
 40   continue

      
C
C     get the matrix-vector product (U' * F2' * y) in y, where U is the
C     orthogonal matrix in the above eigenvalue-eigenvector decomposition
C
      call dgemv('T', n-1, n-1, 1.0d0, Q(1, 2), ldq, y(2), 1,
     $          0.0d0, S(2), 1)

      call dcopy(n-1, S(2), 1, y(2), 1)      


      return
      end







