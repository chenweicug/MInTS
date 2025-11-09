function b = MintsUniformBSpline(n,dtk,t,varargin)
%
% b = MintsUniformBSpline(n,dtk,t,[integral flag])

% INPUT:
% n - order of spline
% dtk - knot separation
% t - vector of sample location
%
% OPTIONS:
% if 1 is specified in the integral flag option, returns the
% integral of b-splines
%
% OUTPUT:
% b - vector of size(t) BSpline basis
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


IFlag = 0;
p = n;
if length(varargin)>=1
  IFlag = varargin{1};
  p = n+1;
end

x = (t./dtk+n+1);
b = zeros(size(x));
for k=0:n+1
  up = (x - k - (n+1)/2).^p;
  b = b + Binomial(n+1,k).*((-1)^k).*...
      up.*((x - k - (n+1)/2)>=0);
end
if IFlag==1
  b = b.*dtk./(n+1);
end
b = b./factorial(n);

return

function [c] = Binomial(p,k)

if ~isempty(find(p-k<0))
   keyboard
end
c = factorial(p)/(factorial(p-k)*factorial(k));

return

%
% update history:
%
% EA Hetland - Caltech 2007
% ehetland@alum.mit.edu
% modified:
% EA Hetland - 02 Jul 2010 (Univ. of Michigan)
%          added option to return the integral of the b-splines
%
% EA Hetland - 16 Aug 2010
%          renamed and included Binomial function
%