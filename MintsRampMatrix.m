function [G] = MintsRampMatrix(X,Y,varargin)
%
% [G] = MintsRampMatrix(X,Y)
%    or
% [G] = MintsRampMatrix(X,Y,N)
%
% default N=4
%
% INPUT:
% X = column vector of x positions of the pixels
% Y = column vector of y positions of the pixels
%
% OPTIONAL INPUT:
% N = number of terms in the ramp (default is 4)
%    ramp = a + b*x + c*y + d*x*y
%           1   2     3     4
%
% OUTPUT:
% G = design matrix of the N-order ramp
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

fcname = 'MintsRampMatrix';

if length(varargin)>=1
  N = varargin{1};
else
  N = 4;
end

% check to make sure positions are vectors
if max(size(X))<length(X(:))
  X = X(:);
end
if max(size(Y))<length(Y(:))
  Y = Y(:);
end
% check to make sure they are column vectors
if size(X,1)==1
  X = X';
end
if size(Y,1)==1
  Y = Y';
end


G = zeros(prod(size(X)),N);
G(:,1) = 1;
if N>=2
  G(:,2) = X;
end
if N>=3
  G(:,3) = Y;
end
if N>=4
  G(:,4) = X.*Y;
end

return

%
% update history:
%
% EA Hetland - Caltech 2007
% ehetland@alum.mit.edu
% modified:
% EA Hetland - 06 Apr 2009 (Univ. of Michigan)
%              renamed to MintsRampMatrix and updated header
%              help
%
%
