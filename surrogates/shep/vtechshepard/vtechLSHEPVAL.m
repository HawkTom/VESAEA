function [ LSHEPVAL, IER ] =  vtechLSHEPVAL( XP, M, N, X, F, A, RW )
% vtechLSHEPVAL returns the linear Shepard approximation at the point XP, using the 
% local linear approximations computed by LSHEP.
%
% Input parameters:
%  XP is the point at which the linear Shepard interpolant function 
%     approximation function is to be evaluated.
%  M is the dimension of the data.
%  N is the number of interpolation points.
%  X contains the interpolation nodes, by column.
%  F contains the function values of the interpolation nodes.
%  A is an M by N matrix containing the coefficients for linear nodal functions
%    returned by LSHEP.
%  RW contains the radius of influence about each interpolation node returned
%    by LSHEP.
% Output parameter
%  LSHEPVAL stores the linear Shepard approximation at the point XP.
%  IER 
%   = 0, normal returns.
%   = 1, if the point XP is outside the radius of influence RW(i) for all nodes,
%        in which case LSHEPVAL is computed using the original Shepard 
%        algorithm with the M+1 closest points.
%   = 2, if the hyperplane normals of the local approximations with positive
%        weights are significantly different. For a nonlinear underlying 
%        function f(x), e.g., quadratic f(x), very different normals are 
%        typical. For a piecewise linear underlying function f(x), IER = 2 
%        signals a potentially large error in LSHEPVAL, since local 
%        approximations from different facets of f(x) have been used.

double precision;

IER = 0;
I = 0;
J = 1;
D(1:M+1) = 0.0;
TEMP = 0.0;
for K = 1: N
   DIST = sqrt( dot( XP(:) - X(:, K), XP(:) - X(:, K) ) );
   if ( DIST < RW(K) )
      if ( RW(K) - DIST == RW(K) ) 
         LSHEPVAL = F(K);
         return
      end
      I = I + 1;
      W = ((RW(K) - DIST) / (RW(K) * DIST))^2;
      SW(I) = W;
      SWC(I) = dot( A(:, K), (XP(:) - X(:, K)) ) + F(K);
      SSLOPE(1:M, I) = A(1:M, K);
      SSLOPE(M+1, I) = -1.0;
   else
      W = 1.0 / DIST;
      if ( W > TEMP )
         D(J) = W;
         ID(J) = K;
% Determines location of minimum value for D(1:M+1).
         TEMP_MIN = D(1);
         for TEMP_COUNT = 1:M+1 
            if(D(TEMP_COUNT) < TEMP_MIN) 
               TEMP_MIN = D(TEMP_COUNT);
               J = TEMP_COUNT;
            end
         end
         TEMP = D(J);
      end
   end
end
% I = 0 iff the point XP is not within the radius RW(K) of X(:,K) for all K; 
% IER = 1.
if ( I == 0 )
   D = D.^2;								
   LSHEPVAL = dot( D / sum( D ), F(ID) );
   IER = 1; 
   return
else
   SW(1:I) = SW(1:I) / sum( SW(1:I) );
   LSHEPVAL = dot( SW(1:I), SWC(1:I) );
% Return IER = 2 iff the angle between two local approximation hyperplanes is 
% too large.
   SLOPE(1:M+1) = SSLOPE(1:M+1, 1:I) * SW(1:I)';
   SLOPE(:) = SLOPE(:) / sqrt( dot( SLOPE(:), SLOPE(:) ) );
   for K = 1: I
      SSLOPE(:, K) = SSLOPE(:, K) / sqrt( dot( SSLOPE(:, K), ...
           SSLOPE(:, K) ) );
      if ( dot( SLOPE(1:M+1), SSLOPE(1:M+1, K) ) < 0.9 )
         IER = 2;
         break
      end
   end
end
return
end
