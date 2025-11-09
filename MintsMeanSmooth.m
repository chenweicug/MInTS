function [varargout] = MintsMeanSmooth(A,rfact,varargin)
%
% B = MintsMeanSmooth(A,rfact,OPTIONS)
% [B,N] = MintsMeanSmooth(A,rfact,OPTIONS)
% 
% Smooths the image (matrix) A, using a moving window.
%
% INPUT:
% A     = matrix
% rfact = number of neighboring pixels to smooth over, window that
%         is smoothed over is (2*rfact+1)x(2*rfact+1) 
%
% OPTIONS:
%
% 'perc',val = only computes the mean at a given pixel location if the
%              window contains at least val% non-NaN values,
%              otherwise sets the value to NaN
%              [default val=0.5]
% 'median'   = uses the median of values in each window
%
% OUTPUT:
% B = smoothed matrix
% N = number of non-NAN pixels used to determine smooth the value
%     at each point
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


% set defaults, and parse the optional input
perc = 0.5;
usemean = 1;
for k=1:length(varargin)
  if ischar(varargin{k})
    switch lower(varargin{k})
     case 'perc'
      perc = varargin{k+1};
     case 'median'
      usemean = 0;
    end
  end
end
      
B = zeros(size(A));
N = zeros(size(A));

for i=1:size(A,1)
  rows = [max([i-rfact 1]):1:min([i+rfact size(A,1)])];
  for j=1:size(A,2)
    cols = [max([j-rfact 1]):1:min([j+rfact size(A,2)])];
    C = A(rows,cols);
    if ((sum(~isnan(C(:)))/numel(C)) < perc)|isnan(A(i,j))
      B(i,j) = NaN;
      N(i,j) = NaN;
    else
      if usemean==1
	B(i,j) = mean(C(~isnan(C)));
      else
	B(i,j) = median(C(~isnan(C)));
      end
      N(i,j) = sum(~isnan(C(:)));
    end
  end
end

varargout{1} = B;
if nargout>=2
  varargout{2} = N;
end

return

%
% update history
%
% EA Hetland - Caltech 2007
% ehetland@alum.mit.edu
% modified:
% EA Hetland - 10 Jun 2008
%       optimized for speed, 18% faster, will not get it faster
%       unless the loops are removed, added the N return, changed
%       rfact to a integer
% EA Hetland - 17 Nov 2010
%       added option for a median smoothing
%