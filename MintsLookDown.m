function [varargout] = MintsLookDown(varargin);
%
% A = MintsLookDown(A,[r1,r2,method,perc])
% [X,Y,A] = MintsLookDown(X,Y,A,[r1,r2,method,perc]);
%
% INPUT:
% A   = igram
% X,Y = arrays of pixel positions
% r1  = the lookdown factor in first dimension of A [default r1=2]
% r2  = the lookdown factor in the second dimension of A 
%      the default is r1
% method = either 'mean' or 'median', the default is 'mean'
% perc = percentage of pixels in each chip that is required to be
% real, otherwise value of chip is NaN [defauly = 0.5];
%
% OUTPUT:
% X,Y,A = looked down arrays
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

fcname = 'MintsLookDown';

r = [2 0];
method = 'mean';
perc = 0.5;
varpos = [2:5];
if length(varargin)>1
  if length(varargin{2})==1
    A{1} = varargin{1};
  else
    A{1} = varargin{1};
    A{2} = varargin{2};
    A{3} = varargin{3};
    varpos = varpos+2;
  end
else
  A{1} = varargin{1};
end

if length(varargin)>=varpos(1)
  r(1) = varargin{varpos(1)};
end
if length(varargin)>=varpos(2)
  r(2) = varargin{varpos(2)};
else
  r(2) = r(1);
end
if length(varargin)>=varpos(3)
  method = varargin{varpos(3)};
end
if length(varargin)>=varpos(4)
  perc = varargin{varpos(4)};
end

if sum(r)==2
  fprintf(1,'%s: not looking down arrays\n',fcname);
  varargout = A;
  return
end

for k=1:length(A)
  varargout{k} = zeros(...
      floor(size(A{1},1)/r(1)),floor(size(A{1},2)/r(2)));
end

fprintf(1,...
	'%s: looking down %d arrays by the %s of %dx%d pixel chip\n',...
	fcname,length(A),method,r(1),r(2));

switch method
 case 'mean'
  UseMean = 1;
 case 'median'
  UseMean = 0;
 otherwise
  fprintf(1,'%s: did not understand method = %s\n',fcname,method);
  return
end

for i=1:floor(size(A{1},1)/r(1))
  rows = [1+(i-1)*r(1):i*r(1)];
  for j=1:floor(size(A{1},2)/r(2))
    cols = [1+(j-1)*r(2):j*r(2)];
    for k=1:length(A)
      C = A{k}(rows,cols);
      if (sum(~isnan(C(:)))/numel(C)) < perc
	varargout{k}(i,j) = NaN;
      else
	if UseMean
	  varargout{k}(i,j) = mean(C(~isnan(C)));
	else
	  varargout{k}(i,j) = median(C(~isnan(C)));
	end
      end
    end
  end
end


return

%
% update history:
%
% EA Hetland - Caltech 2007
% ehetland@alum.mit.edu
% modified:
% EA Hetland - 12 June 2008
%              overloaded to handle three arrays, and changed
%              averaging to not propagate NaN's
% EA Hetland - 17 Feb 2009 (Univ. of Michigan)
%              renamed and cleaned up
%
