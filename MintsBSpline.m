function [B] = MintsBSpline(n,tk,t)
%
% constructs the basis splines for a given knot vector, the basis
% associated with the i'th knot has support [tk_i,b], where b = n+1;
% total support of the BSpline set is [tk(1),tk(end)]
%
% INPUT:
% n = spline order
% tk = knot locations
% t = vector of sample points
%
% OUTPUT:
% B = matrix of splines B(:,k) = basis associated with k'th knot
%
% reference:
% Hetland, E. A., P. Muse, M. Simons, Y. N. Lin, P. S. Agram, and
% C. J. DiCaprio (2012), Multiscale InSAR Time Series (MInTS) analysis
% of surface deformation. J. Geophys. Res., 117, B02404,
% doi:10.1029/2011JB008731
%
% MInTS is free for use or modification, but please retain the
% above reference in any modified versions of this code. Please
% report any bugs or contribute modifications of the code to:
% Eric A. Hetland, Univ. Michigan, ehetland@alum.mit.edu
% Feb 2012
%

fcname = 'MintsBSpline';

if size(t,1)==1
  t = t';
end

B = zeros(length(t),length(tk)-1);
if n+1>=length(tk)
  fprintf(1, '%s: not enough knots for an order %d spline',fcane,n));
  return
end

i=0;
for j=0:length(tk)-2
  B(find((tk(j+1)<=t)&(t<tk(j+2))),j+1) = 1;
end

while i<n
  i=i+1;
  for j=0:length(tk)-2-i
    B(:,j+1) = (t - tk(j+1))./(tk(j+i+1)-tk(j+1)).*B(:,j+1) + ...
	(tk(j+i+2) - t)./(tk(j+i+2)-tk(j+2)).*B(:,j+2);
  end
end
B = B(:,1:length(tk)-1-n);


return

%
% update history:
%
% EA Hetland - Caltech 2007
% ehetland@alum.mit.edu
%
