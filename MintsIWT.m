function [varargout] = MintsIWT(WaveCoef,varargin)
%
% [image] = MintsIWT(WaveCoefs,OPTIONS)
%
% Determines a spatial map based on wavelet coefficients.
%
% INPUT:
% WaveCoeff    = cell array of all the wavelet coefficients of the
%                igrams, WaveCoeff{k,l}(:,:,n) are the scale k,
%                band l, coefficients for igrams n; k=mstg+1 is
%                only defined for l=1
%
% OPTIONS:
% 'meyer'      = uses the 2D Meyer wavelets (default)
% 'farras'     = uses the 2D Farras wavelets
% 'scale',scl  = uses scl=[min max] as the range of scales to use
%                in the reconstruction of the images, if scl=[max],
%                then uses min=1
%                [default mscl=3]
% 'window',win = win=[min_row max_row min_col max_col] the window of
%                the original igram geometry.
% 'index',ind = only transforms the images in vector ind
%
% OUTPUT:
% image        = 3D array of all the reconstructed images
%
% See also:
% MintsTimeString
% Requires:
% WaveLab850 toolbox from Stanford, http://www-stat.stanford.edu/~wavelab/
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

fcname = 'MintsIWT';

meyerdeg = 1;
meyerl = 3;
ninds = [1:size(WaveCoef{1,1},3)];
%%DIM = size(WaveCoef{1,1}(:,:,1)).*2;
mscl = size(WaveCoef,1);
scale = [1 mscl];
WIN = [];
type = 'meyer';
for v=1:length(varargin)
  if ischar(varargin{v})
    switch lower(varargin{v})
     case 'farras'
      type = 'farras';
      fprintf(1,['\n%s: NOTE Farras wavelets are not fully ', ...
		 'verified in current version\n\trecommend using', ...
		 ' Meyer wavelets\n'],fcname);      
     case 'meyer'
      type = 'meyer';
     case 'index'
      ninds =varargin{v+1};
      if size(ninds,1)>1
	ninds = ninds';
      end
     case 'window'
      WIN = varargin{v+1};
     case 'scale'
      if length(varargin{v+1})==1
	scale(2) = varargin{v+1};
      else
	scale = varargin{v+1}(1:2);
      end
    end
  end
end

N = length(ninds);

% initialize array of times for checkpoints
% times are [start_of_computation;end_of_IDWT]
timing = zeros(2,6);
timing(1,:) = clock;


if scale(2) > mscl
  scale(2) = mscl;
end

switch type
 case 'farras'
  [af,sf] = farras;
end

fprintf(1,'%s: determining %d images from scale %d to %d of %d scales\n',...
	fcname,N,scale,mscl);

% find the dimensions of the scales
NN = zeros([size(WaveCoef,1) 2]);
for k=1:mscl-1
  if ~isempty(WaveCoef{k,1})
    NN(k,1:2) = size(WaveCoef{k,1}(:,:,1));
  end
end
NN(end,1:2) = size(WaveCoef{end,1}(:,:,1));
for l=1:2
  for k=1:mscl-1
    if NN(k,l)==0
      [dm ps] = max(NN(:,l));
      NN(k,l) = dm*2^(ps-k);
    end
  end
  if NN(end,l)==0
    NN(end,l) = NN(end-1,l);
  end
end
DIM = NN(1,1:2).*2;
if isempty(WIN)
  WIN = [1 DIM(1) 1 DIM(2)];
end

% construct images
%varargout{1} = zeros([DIM N]);
varargout{1} = zeros([diff(WIN(1:2))+1 diff(WIN(3:4))+1 N]);
for n=ninds
  fprintf(1,['%s: constructing image %d\n'],fcname,n);
  coef = cell([1 size(WaveCoef,1)]);
  for k=1:mscl-1
    coef{k} = cell([1 3]);
    for l=1:3
      if k>=scale(1)&k<=scale(2)
	coef{k}{l} = WaveCoef{k,l}(:,:,n);
      else
	coef{k}{l} = zeros(NN(k,1:2));
      end
    end
  end
  if scale(2)==mscl
    coef{end} = WaveCoef{end,1}(:,:,n);
  else
    coef{end} = zeros(NN(end,1:2));
  end
  switch type
   case 'farras'
    fprintf(1,'%s: reconstructng image with inverse Farras wavelet transform\n');
    varargout{1}(:,:,n) = idwt2D(coef,mscl-1,sf);
   case 'meyer'
    %fprintf(1,'%s: re-ordering wavelet coefficients into Meyer format\n',fcname);
    f = InvDealMeyerCoeffs(coef,meyerl);
    fprintf(1,'%s: reconstructng image with inverse Meyer wavelet transform\n',fcname);
    %varargout{1}(:,:,n) = IWT2_YM(f,meyerl,meyerdeg);
    at = IWT2_YM(f,meyerl,meyerdeg);
    varargout{1}(:,:,n) = at(WIN(1):WIN(2),WIN(3):WIN(4));
  end
end

%extract igram geometry
%Tvarargout{1} = varargout{1}(WIN(1):WIN(2),WIN(3):WIN(4),:,:);

timing(2,:) = clock;

fprintf(1,'%s: computation time = %s\n',fcname,MintsTimeString(...
    etime(timing(2,:),timing(1,:))));
	

return

function [ff] = InvDealMeyerCoeffs(f,L);
% puts wavelet coefficients back into the array expected by meyer
% wavelet transform programs
% L = Meyer DWT parameter


if size(f,1)==1|size(f,2)==1
  J = length(f);
else
  J = size(f,1);
end

ff = zeros([1 1].*2.^(J+2));


cnt=0;
for j=log2(size(ff,1))-1:-1:L
  cnt=cnt+1;
 
  %f{cnt} = cell([1 3]);
  
  r1 = [1:2^j];
  r2 = [2^j+1:2^j+2^j];

  if size(f,1)==1|size(f,2)==1
    ff(r1,r2) = f{cnt}{1};
    ff(r2,r2) = f{cnt}{2};
    ff(r2,r1) = f{cnt}{3};
  else
    ff(r1,r2) = f{cnt,1};
    ff(r2,r2) = f{cnt,2};
    ff(r2,r1) = f{cnt,3};
  end
end
if size(f,1)==1|size(f,2)==1
  ff(r1,r1) = f{cnt+1};
else
  ff(r1,r1) = f{cnt+1,1};
end

return

%
% update history:
%
% EA Hetland - Univ. of Michigan 14 Feb 2009
% ehetland@alum.mit.edu
% modified:
% EA Hetland - 05 Jun 2009
%              cleaned script
% EA Hetland - 16 Jun 2009
%            fixed bug to allow empty matricies for unused wavelet
%            coeficients
% EA Hetland - 21 Apr 2010
%            added meyer wavelet transforms
% EA Hetland - 05 Aug 2010
%            incoperated InvDealmeyerCoeffs function into this file
% EA Hetland - 16 Nov 2010
%            replaced cputime with clock for timing functions
% EA Hetland - 03 Jun 2011
%            fixed trivial typo of a line continuation
% EA Hetland - 06 Jun 2011
%            corrected help message
% EA Hetland - 28 Jun 2011
%            added index option
%
