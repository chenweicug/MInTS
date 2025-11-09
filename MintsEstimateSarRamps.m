function [varargout] = MintsEstimateSarRamps(x,y,IGram,varargin)
%
% RAMP = MintsEstimateSarRamps(X,Y,IGram,[PN])
%
% Independently estimates the best fitting ramp from each
% interferogram. The (x,y) positions and the sar observations are
% normalized prior to inverting for the best fit ramp, and ramp
% coefficients are re-normalized before they are returned.
%
% INPUT:
%     X = 1D or 2D array of x positions of igrams
%     Y = 1D or 2D array of x positions of igrams
% IGram = arrays of igrams, 
%         size(IGram) = [Nx Ny M] or [N M] or [M N], for M igrams
%    PN = number of coefficients in the ramp (default PN=4)
%
% OUTPUT:
% RAMP = ramp coefficients size(RM) = [PN M]
%
% See also:
% MintsRampMatrix
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


fcname = 'MintsEstimateSarRamps';
if length(varargin)>=1
  PN = varargin{1};
else
  PN = 4;
end

% normalize X and Y
if size(x,1)~=prod(size(x));
  X = x(:);
end
if size(x,1)==1
  X = x';
end
if size(y,1)~=prod(size(y));
  Y = y(:);
end
if size(y,1)==1
  Y = y';
end

if size(IGram,3)==1
  if size(IGram,1)==length(X)
    drho = IGram;
  elseif size(IGram,2)==length(X)
    drho = IGram';
  elseif prod(size(IGram))~=length(X)
    fprintf(1,['%s: the number of pixels in the IGrams does not equal the' ...
	       ' number of specified positions\n'],fcname);
    varargout{1:nargout} = 0;
    return
  else
    drho = IGram(:);
  end
else
  drho = reshape(IGram,[prod(size(IGram(:,:,1))) size(IGram,3)]);
end

for n=1:size(drho,2)
  fprintf(1,['%s: estimating a %d-order ramp from igram %d of %d\n'],...
	  fcname,PN,n,size(drho,2));
  pos = find(~isnan(drho(:,n)));
  d = drho(pos,n);
  x = X(pos);
  y = Y(pos);

  MuX = mean(x);
  SiX = std(x);
  MuY = mean(y);
  SiY = std(y);
  x = (x-MuX)./SiX;
  y = (y-MuY)./SiY;

  MuR = mean(d);
  SiR = std(d);
  d = (d-repmat(MuR,[size(d,1) 1]))./...
      repmat(SiR,[size(d,1) 1]);

  G = MintsRampMatrix(x,y,PN);
  Gg = inv(G'*G)*G';

  m = Gg*d;
  % re-normalize ramp coefficients
  if PN==1
    m(1) = m(1)*SiR+MuR;
  elseif PN==2
    m(2) = m(2)*SiR/SiX;
    m(1) = m(1)*SiR - m(2)*MuX + MuR;
  elseif PN==3
    m(3) = m(3)*SiR/SiY;
    m(2) = m(2)*SiR/SiX;
    m(1) = m(1)*SiR - m(2)*MuX - m(3)*MuY + MuR;
  elseif PN==4
    m(4) = m(4)*SiR/(SiX*SiY);
    m(3) = m(3)*SiR/SiY - MuX*m(4);
    m(2) = m(2)*SiR/SiX - MuY*m(4);
    m(1) = m(1)*SiR - m(2)*MuX - m(3)*MuY - m(4)*MuX*MuY + MuR;
  else
    fprintf(1,['%s: have not accounted for re-normalization' ...
	       ' for ramps >4 parameters\n'],fcname);
  end    
  RM(:,n) = m;
end

varargout{1} = RM;

return

%
% update history
%
% EA Hetland - Caltech 2007
% ehetland@alum.mit.edu
% modified:
% EA Hetland - 12 Jun 2008
%    renamed to EstimateSarRamps, and updated scipt to only do the
%    estimation
% EA Hetland - 05 Apr 2009 (Univ. of Michigan)
%    renamed to MintsEstimateSarRamps, cleaned help header, updated
%    to estimate ramps on a stack of igrams
% EA Hetland - 05 Jun 2009
%    cleaned script
%
