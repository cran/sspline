C***************************************************************************
C!   deqsv.f
C***************************************************************************

      subroutine deqsv(Q, ldq, n, x, eigen, lambda, coefc, coefd, 
     $                 S, work1, lwork1, work2)

C
C     Purpose:
C     ========
C
C         Solve the equation based on the information returned by dgcvst
C         and calculate the estimated y value.
C
C     Arguments:
C     ==========
C
C       Q: input, double array, dimension ldq by n, on entry
C          Q(1:n-1, 1) is stored with F2'*Q*F1,
C          Q(1:n-1, 2:n) is stored with the orthogonal matrix of
C          the eigenvalue-eigenvector decomposition of F2'*Q*F2
C       ldq: input, integer, leading dimension of Q
C       n: input, integer, number of columns of Q
C       x: input, double array, dimension n by 1, on entry,
C          x(1) is stored with F1'*y, x(2:n) is stored with U'*F2'*y
C       eigen: input, double array, dimension (n-1) by 1, on entry,
C              the eigenvalues of F2'*Q*F2.
C       lambda: input, double, the optimal lambda chosen by GCV
C
C
C       coefc: output, double array, dimension n by 1, on exit, the
C              coefficient c
C       coefd: output, double, returned with coefficent d 
C
C
C       S: workspace, double array, dimension n by 1
C       work1: workspace, double array, dimension lwork by 1
C       lwork1: input, integer, the leading dimension of work
C       work2: workspace, double array, dimension n by 1
C

C
C     declaration of the arguments
C     
      integer n, ldq, lwork1
      double precision lambda, x(n), eigen(n-1), Q(ldq, n) 
      double precision coefd, coefc(n), S(n), work1(lwork1)
      double precision work2(n)

C
C     declaration of the local variables
C
      integer  info, jpvt(1)
      double precision  tau(1)


C
C     get the matrix-vector product (T + n*lambda*I)^{-1} * x in coefc, 
C     where T is the (n-1) dimensional diagonal matrix with the diagonal 
C     elements eigen, the product coefc has the dimension (n-1) by 1
C
      do 10, i = 1, n-1
         coefc(i) = x(i+1) / (eigen(i) + dble(n)*lambda)
 10   continue

C
C     get the matrix-vector product (U * coefc) (n-1 by 1) in work2,
C     note that U is stored in Q from the second column.
C
      call dgemv('N', n-1, n-1, 1.0d0, Q(1, 2), ldq, coefc, 1, 0.0d0, 
     $           work2, 1)

C
C     copy the product from work2 to coefc for later use
C
      call dcopy(n-1, work2, 1, coefc(2), 1)
      coefc(1) = 0.0d0

C
C     construct an all 1 vector S
C
      do 20, i = 1, n
         S(i) = 1.0d0
 20   continue

C
C     do the Q-R decomposition on the all 1 vector S (n by 1)
C
      jpvt(1) = 0

      call dgeqp3(n, 1, S, n, jpvt, tau, work1, lwork1, info)

C
C     get the product F2 * (U * (T + n*lambda*I)^{-1} * x) in coefc,
C     note the fact that (F1 F2) * (0.0 coefc(2:n))' = F2 * coefc(2:n)'
C
      call dormqr('L', 'N', n, 1, 1, S, n, tau, coefc, n,
     $            work1, lwork1, info)


C
C     get the inner product of F2'*Q*F1 and U*(T + n*lambda*I)^{-1}*x,  
C     and calculate coefd according to the following formula, note that
C     work2(1:n-1) is still stored with (U*(T + n*lambda*I)^{-1}*x).
C
C       coefd = R1^{-1}*(F1'*y - (F1'*Q*F2)*(U*(T + n*lambda*I)^{-1}*x))
C

      coefd = (x(1) - ddot(n-1, Q(1, 1), 1, work2, 1)) / S(1)


      return
      end











