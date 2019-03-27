function [A, RW, IER] = vtechLSHEP(M, N, X, F, RLSHEP)
% This subroutine computes a set of parameters defining a function that 
% interpolates N data values F(i) at scattered nodes X(i) in M dimensions. The 
% interpolant may be evaluated at an arbitrary point by the function LSHEPVAL.
% The interpolation scheme is a modified linear Shepard method, and 
% can use either the Shepard weighted least squares fit or robust 
% M-estimation for each local approximation.
%
% Input parameters:  
%   M is the dimension of the data.
%   N is the number of nodes.
%   X(M, N) contains the coordinates of the nodes.
%   F(N) contains the function values of the nodes.
% Output parameters:
%   A(M, N) contains the coefficients of the local fits.
%   RW(N) is an array containing the radius of influence for each point.
%   IER 
%   = 0, if no errors were encountered.
%   = 1, if N is too small relative to M.
%   = 2, if any least squares problem is rank deficient.
%   = 3, if the IRLS subroutine returns an error.
% Optional input:
%   RLSHEP specifies by its presence that robust M-estimation is to be used to 
%          compute the local approximation.

double precision;

IER1 = 0;
DIST(1:N-1) = 0.0;
IDIST(1:N-1) = 0.0;
TEMPVAR = min(N-1, floor((3*M+1)/2));
TEMPVAR = min(TEMPVAR, M);
ALPHA(1:TEMPVAR) = 0.0;
ALPHA_IN(1:TEMPVAR) = 0.0;
A(1:M, 1:N) = 0.0;
RW(1:N) = 0.0;

IER = 0;
% Check the validity of M and N
if  ( ( M < 1) || (N <= M+1) ) 
    IER = 1;
    return
end
% Set the number of nodes used in local linear fit
NP = min((N-1), floor((3 * M + 1)/2));
% Set the RCOND parameter used in DGELSS
RCOND = NP * eps(double(1.0));
%Calculate RW and A
DIAM = 0.0;
for K = 1:N
    while (true)
    J=0;    
    for I = 1:N
        if (I ~= K)
            J = J + 1;
            DIST(J) = dot(X(:, I) - X(:, K), X(:,I) - X(:,K));
            IDIST(J) = I;
        end
    end
    DIAM = max(DIAM, max(DIST(1:N-1)));
    [DIST, IDIST] = lshepsort(DIST, IDIST, NP, N-1);
    DIST(1:NP) = sqrt(DIST(1:NP));
    BETA_IN(1:NP) = F(IDIST(1:NP)) - F(K);
    for I = 1:NP
        ALPHA_IN(I, 1:M) = X(1:M, IDIST(I)) - X(1:M, K);
    end
    RW(K) = DIST(NP);
    if (nargin == 5)
        OMEGA(1:NP) = 1.0;
    else
        RP = 1.1 * RW(K);
        OMEGA(1:NP) = (RP-DIST(1:NP)) ./ (RP * DIST(1:NP));
    end
    BETA(1:NP) = OMEGA(1:NP) .* BETA_IN(1:NP);
    for I = 1: NP
        ALPHA(I, 1:M) = OMEGA(I) * ALPHA_IN(I, 1:M);
    end
    [ALPHAOUT, BETAOUT, S, RANK] = DGELSS(ALPHA(1:NP, 1:M), ...
        BETA(1:NP), RCOND);
    ALPHA = ALPHAOUT;
    BETA = BETAOUT';
    A(1:M, K) = BETA(1:M);
    if (RANK < M) 
        IER = 2;
    end
    if (nargin == 5)
        RES(1:NP) = (ALPHA_IN(1:NP, 1:M) * BETA(1:M)') - BETA_IN(1:NP)'; 
        MAR = MAD(RES(1:NP), NP);                                     
        if (MAR > eps) 
            FUNC = 1;
            ITER = 5;
            [MAR, OMEGA, BETA, ALPHA, RES, IER1] = IRLS(M, N, ... 
                MAR, ITER, FUNC, BETA_IN, ...
                ALPHA_IN, RES, NP);
            if (IER1 ~= 0)
                IER = 3;
                break;
            else
                A(1:M, K) = BETA(1:M);
            end
            MAR = MAD(RES(1:NP), NP);
            if (MAR > eps)
                FUNC = 2;
                ITER = 5;
                [MAR, OMEGA, BETA, ALPHA, RES, IER1] = IRLS(M, N, ... 
                    MAR,ITER, FUNC, BETA_IN, ...
                    ALPHA_IN, RES, NP);
                    if (IER1 ~= 0)
                        IER = 3;
                        break;
                    else
                        A(1:M, K) = BETA(1:M);
                    end
            end
         end
          for I = 1:NP
            if (OMEGA(I) < 0.8)
                  if (I == 1)
                      RW(K) = DIST(I) / 2.0;
                   else
                      RW(K) = ( DIST(I-1) + DIST(I) ) / 2.0;
                   end
                 break;
            end
          end
    end
            break;
    end
end
    RW(1:N) = min( (sqrt(DIAM) / 2.0), RW(1:N) );
return
end


function [dist, idist] = lshepsort(dist, idist, num, length)
%The subroutine SORT sorts the real array DIST of length LENGTH in
%ascending order for the smallest NUM elements.  IDIST stores the original
%label for each element in array DIST

double precision;
for i = 2:length
    for j = 1:min(i-1, num)
        if (dist(i) < dist(j))
            temp = dist(i);
            itemp = idist(i);
            dist(j+1:i) = dist(j:i-1);
            idist(j+1:i) = idist(j:i-1);
            dist(j) = temp;
            idist(j) = itemp;
        end
    end
end
end

function [S, OMEGA, BETA, ALPHA, RES, IER] = IRLS( M, N,  S, ITER, ...
    FUNC, BETA_IN, ALPHA_IN, RES, NP)
% IRLS returns the coefficients of the linear least squares fit at the current 
% node using robust M-estimation. For Tukey's bisquare influence function, 
% there is an additional convergence test.
%
% Input parameters:
%  M    is the dimension of the problem.
%  N    is the number of points.
%  S    is the scale constant.
%  ITER is the number of iterations to be used in the algorithm.
%  FUNC identifies which IRLS influence function to use here, 1 = Huber,
%       2 = Tukey's bisquare.
%  RES is the residuals from the linear fit.
%  BETA_IN holds distances.
%  ALPHA, ALPHA_IN are used for DGELSS calls.
% Output parameters:
%  OMEGA stores the final weights.
%  BETA stores the coefficients of the linear fit using robust M-estimation.
%  NP is the size of OMEGA and BETA.
%  IER 
%   = 0, if no errors were encountered.
%   = 1, if any least squares problem is rank deficient.
%   = 2, if the bisquare IRLS did not converge.
%   = 3, if the input parameters are invalid.
double precision;

TUNING = [1.0 3.0];
if ((FUNC ~= 1) && (FUNC ~=2))
    IER = 3;
    return;
end
S = TUNING(FUNC) * S;
%  Set the RCOND parameter used in DGELSS
RCOND = NP * eps(double(1.0));
if ( FUNC ==  2 )
   ERR_OLD = sum(1.0 - (max(0.0, 1.0 - (RES(1:NP) ./  S).^2 )).^3 );
end
for I = 1:ITER
% Update the weight from the previous linear fit residual.
   switch( FUNC )
% Huber function
       case(1) 
           for J = 1:NP
               if ( abs(RES(J)) < eps(double(1.0)) ) % if ( RES(J) == 0 )
                   OMEGA(J) = 1.0; 
               else
                   OMEGA(J) = sqrt(min( 1.0, S / abs(RES(J)) ));
               end
           end
% Bisquare function
       case(2)
           OMEGA(1:NP) = max( 0.0, 1.0 - (RES(1:NP) / S).^2);
   end
% Compute the linear fit with a new weight.
   BETA(1:NP) = OMEGA(1:NP) .* BETA_IN(1:NP);
   for J = 1:NP
      ALPHA(J, 1:M) = OMEGA(J) .* ALPHA_IN(J, 1:M);
   end
   [ALPHA(1:NP, 1:M), BETAOUT, JPVT, RANK] = ...
       DGELSS(ALPHA(1:NP, 1:M), BETA(1:NP), RCOND);
   BETA(1:M) = BETAOUT(1:M);
   if ( RANK < M )
      IER = 1;
      return;
   end
   RES(1:NP) = ALPHA_IN(1:NP, 1:M) * BETA(1:M)' - BETA_IN(1:NP)';
end
if ( FUNC ==  2 )
   ERR_NEW = sum( 1.0 - (max( 0.0, 1.0 - (RES(1:NP) ./ S).^2 )).^3 );
   if ( ERR_NEW > ERR_OLD )
      IER = 2;
      return;
   end
end
IER = 0;
return;
end

function VALUE = MAD( RES, NP)
% The function MAD computes the median absolute deviation (MAD) scaled estimate 
% for robust M-estimation from the residual array RES.
double precision;
  
FACTOR = 1.0 / 0.6745;
NMID = floor(NP / 2);
if( mod(NP, 2) == 0 )
   MED = (QUANTILE( RES(:), NMID ) + QUANTILE( RES(:), NMID + 1 )) / 2.0;
   VALUE = FACTOR * (QUANTILE( abs( RES(:) - MED ), NMID ) + ... 
        QUANTILE( abs( RES(:) - MED ), NMID + 1 )) / 2.0;
else
   MED = QUANTILE( RES(:), NMID + 1 );
   VALUE = FACTOR * QUANTILE( abs( RES(:) - MED ), NMID + 1 );
end
return
end


function VALUE = QUANTILE(A, K)
%QUANTILE returns the K-th smallest element in array A.
double precision;

Z = sort(A);
VALUE = Z(K);
return;

end


%function VALUE = QUANTILE( A, K )
%! Recursive function QUANTILE returns the K-th smallest element in array A.
%double precision;
%
%AK = A(K);
%J = size(find(A < AK), 2);
%if ( J >= K )
%   B = A < AK;
%   VALUE = QUANTILE( A(B), K );
%else
%   J = size(find( A > AK ), 2) + K - size(A, 2);
%   if ( J > 0 )
%       B = A > AK;
%       VALUE = QUANTILE( A(B), J );
%   else
%      VALUE = AK;
%   end
%end
%return
%end

function[AOUT, BOUT, SOUT, RANKOUT] = DGELSS(A, ...
    B, RCOND)

%Interim notes:
%  Work, Lwork have been eliminated from DGELSS calls.
% Purpose
%  =======
%
%  DGELSS computes the minimum norm solution to a real linear least
%  squares problem:
%
%  Minimize 2-norm(| b - A*x |).
%
%  using the singular value decomposition (SVD) of A. A is an M-by-N
%  matrix which may be rank-deficient.
%
%  Several right hand side vectors b and solution vectors x can be
%  handled in a single call; they are stored as the columns of the
%  M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix
%  X.
%
%  The effective rank of A is determined by treating as zero those
%  singular values which are less than RCOND times the largest singular
%  value.
%
%  Arguments
%  =========
%
%  M       (input) INTEGER
%          The number of rows of the matrix A. M >= 0.
%
%  N       (input) INTEGER
%          The number of columns of the matrix A. N >= 0.
%
%  NRHS    (input) INTEGER
%          The number of right hand sides, i.e., the number of columns
%          of the matrices B and X. NRHS >= 0.
%
%  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
%          On entry, the M-by-N matrix A.
%          On exit, the first min(m,n) rows of A are overwritten with
%          its right singular vectors, stored rowwise.
%
%  LDA     (input) INTEGER
%          The leading dimension of the array A.  LDA >= max(1,M).
%
%  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
%          On entry, the M-by-NRHS right hand side matrix B.
%          On exit, B is overwritten by the N-by-NRHS solution
%          matrix X.  If m >= n and RANK = n, the residual
%          sum-of-squares for the solution in the i-th column is given
%          by the sum of squares of elements n+1:m in that column.
%
%  LDB     (input) INTEGER
%          The leading dimension of the array B. LDB >= max(1,max(M,N)).
%
%  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
%          The singular values of A in decreasing order.
%          The condition number of A in the 2-norm = S(1)/S(min(m,n)).
%
%  RCOND   (input) DOUBLE PRECISION
%          RCOND is used to determine the effective rank of A.
%          Singular values S(i) <= RCOND*S(1) are treated as zero.
%          If RCOND < 0, machine precision is used instead.
%
%  RANK    (output) INTEGER
%          The effective rank of A, i.e., the number of singular values
%          which are greater than RCOND*S(1).
%
%  INFO    (output) INTEGER
%          = 0:  successful exit
%          < 0:  if INFO = -i, the i-th argument had an illegal value.
%          > 0:  the algorithm for computing the SVD failed to converge;
%                if INFO = i, i off-diagonal elements of an intermediate
%                bidiagonal form did not converge to zero.
%
%  =====================================================================

% This version of DGELSS uses the MATLAB SVD function to do a singular value
% decomposition of A.  Testing shows that this function finds the same
% values as the Fortran DGELSS.  It then uses this decomposition to compute
% the Moore-Penrose Pseudo-Inverse of A, pinv(A).  x = pinv(A)*b
% minimizes the 2-norm of |b - Ax|.  For the tolerance of the Moore-Penrose
% Pseuodoinverse, the same tolerance is used as in the original, RCOND
% times the largest singular value.

double precision;

temp = size(A);
M = temp(1);
N = temp(2);
B = B'; %  B is a row vector coming in, needs to be a column vector.
temp = size(B);
NRHS = temp(2);
[U, S, V] = svd(A);    
if (min(size(S)) > 1)
    SOUT = diag(S);
else
    SOUT = S;
end
%SOUT holds the singular values of A in descending order
S1 = max(max(SOUT));
if (RCOND <= 0)
    tol = max(size(A)) * eps(max(S));
else
    tol = RCOND * S1;
end
RANK = sum(SOUT > tol);
%Having taken the SVD of A, we can use that to find the pseudoinverse
% of A and use that to calculate the X that minimizes the l2-norm 
% (|b-Ax|), the least squares solution.

%Algorithm for pinv() function within Matlab
E = S';
S2 = E;
for i = 1:min(size(E))
    if E(i,i) > tol
        S2(i, i) = 1/E(i, i);
    else
        S2(i, i) = 0;
    end
end
pinvA = V * S2 * U';
X = pinvA * B;
%Residual sum of squares
 if (M >= N)
     if (RANK == N)
         %for i = 1:NRHS
         %    temp = B(:,i) - A*X(:,i);
         %    sumsq = 0;
         %    for j = 1:N
         %    sumsq = temp(j)^2 + sumsq;
         %    end
         %    sumsq = sqrt(sumsq);
         %    B(N+1, i) = sumsq;
         %    if (N+1 < M)
         %        B(N+2:M, i) = 0;
         %    end
         %end
     %else B(N+1:M, :) = 0;
     B(N+1) = norm(A(1:M,1:N)*X(1:N)-B(1:M));          % Fixes the code above that is commented out.
     end
 end
%end residual sum of squares
%B(N+1) = norm(B(1:M)-A(1:M,1:N)*X(1:N));  
B(1:N, 1:NRHS) = X(1:N, 1:NRHS);
BOUT = B;

RANKOUT = RANK;
%The first min(m,n) rows of A are overwritten with the right singular
%vectors, stored rowwise
V = V';
MINMN = min(M,N);
temp = size(V);
temp2 = temp(2);
A(1:MINMN, 1:temp2) = V(1:MINMN, 1:temp2);
AOUT = A;
end
