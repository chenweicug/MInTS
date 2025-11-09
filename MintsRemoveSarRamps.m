function [varargout] = MintsRemoveSarRamps(X,Y,IGrams,varargin)
%
% IGrams = MintsRemoveSarRamps(X,Y,IGrams,[RAMP])
%
% Removes the best fitting ramp from each interferogram. The (x,y)
% positions and the sar observations are normalized prior to
% inverting for the best fit ramp, and ramp coefficients are
% re-normalized before they are returned.
%
% INPUT:
%     X  = 1D or 2D array of x positions of igrams
%     Y  = 1D or 2D array of x positions of igrams
% IGrams = arrays of igrams, 
%          size(IGram) = [Nx Ny M] or [N M] or [M N], for M igrams
% RAMP   = PN or ramp coefficients, size(RM) = [PN M]
%          if RAMP not defined, estimates the RAMP from the igrams
%          [default RAMP=4]
%
% OUTPUT:
% IGrams = IGrams with ramp removed
%
% See also:
% MintsEstimateSarRamps, MintsRampMatrix
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

fcname = 'MintsRemoveSarRamps';
D1 = size(IGrams,1);
D2 = size(IGrams,2);
D3 = size(IGrams,3);
if length(varargin)>=1
  RAMP = varargin{1};
else
  RAMP = 4;
end
if prod(size(RAMP))==1
  PN = RAMP;
  RAMP = MintsEstimateSarRamps(X,Y,IGrams,PN);
end

PN = size(RAMP,1);
x = X(:);
y = Y(:);
for n=1:size(IGrams,3)
  fprintf(1,['%s: removing a %d-order ramp from igram %d of %d\n'],...
	       fcname,PN,n,size(IGrams,3));
  rho = reshape(IGrams(:,:,n),[prod(size(IGrams(:,:,n))) 1]);
  G = MintsRampMatrix(X,Y,PN);
  IGrams(:,:,n) = reshape(rho - G*RAMP(:,n),[D1 D2]);
end

varargout{1} = IGrams;
varargout{2} = RAMP;

return

%
% update history
%
% EA Hetland - Caltech 12 Jun 2008
% ehetland@alum.mit.edu
% modified:
% EA Hetland - 06 Apr 2009 (Univ. of Michigan)
%              renamed to MintsRemoveSarRamps.m, and updated to
%              operate on a stack of igrams, instead of removing
%              the same ramp from all
% EA Hetland - 05 Jun 2009
%              cleaned script
%