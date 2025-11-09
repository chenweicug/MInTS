function [varargout] = MintsWT(A,varargin)
%
% [WaveCoefs] = MintsWT(A,OPTIONS)
% [WaveCoefs,Window] = MintsWT(A,OPTIONS)
% [WaveCoefs,Window,CoefWeights] = MintsWT(A,OPTIONS)
%
% Determines the wavelet coefficients for a stack of
% interferograms. Set decorrelated pixels to NAN. The discrete wavelet
% transform of each igram is determined after the igrams are
% optionally mirrored to avoid edge effects, expanded to a base-2
% dimensions, and all decorrelated holes are filled. Coefficient
% weights are based on the amount of decorrelated regions supporting
% each coefficient. For interferograms with a relatively small area of
% decorrelation, the computation may be significantly faster if you
% use the compliment of the mask (CompMask=1), which is the
% default. For images with relatively large decorrelated areas, the
% computation is faster using the mask directlt.
%
% INPUT:
% A = array of igrams, size(A) = [num_rows num_cols num_igrams]
%
% OPTIONS: (keywords are not case sensitive!)
% 'farras'     = uses the 2D Farras wavelets (must install Selesnick's
%                                             wavelet toolbox)
%               **farras wavelets are no longer supported, and
%               should be used with extreme caution**
% 'meyer'      = uses the 2D Meyer wavelets (default; must install Wavelab850)
% meyer specific options:
% 'scale',scl  = either [minscale maxscale] or [maxscale] uses maxscale
%                as the maximum wavelet scale, and minscale as the minumum
%                [default scl=[1 4]]
% 'MeyerDeg'   = the degree option of the meyer wavelet (see help
%                for FWT2_YM)
%                [default MeyerDeg=1]
%                
% farras specific options:
% 'scale',scl  = either [minscale maxscale] or [maxscale] uses maxscale
%                as the maximum wavelet scale, and minscale as the minumum
%                [default scl=[1 4]]
% 'stage',stg  = either [minstage maxstage] or [maxstage] uses maxstage
%                as the maximum wavelet transform stage, and minstage
%                as the minumum, (maxscale = maxstage+1);
%                specification of stages over-rides MaxScale above
%                [default stg=[1 3]]
%
% 'mirror',mf  = mirrors mf percent of the igrams to avoid edge effects
%                [default mf=0.5]
% 'PreFill',lf = if lf=1 fills the NAN's prior to mirroring
%                if lf=0 fills the NAN's after mirroring
%                [default lf=1] if lf not given, assumes lf=1
% 'CompMask',lf= if lf=1 computes the weights of the wavelet
%                coefficients using the compliment of the mask
%                if lf=0 uses the normal order mask
%                [default lf=1] if lf not given, assumes lf=1
% 'ramp',PN    = removes a PN component ramp from igrams prior to
%                mirroring, filling NAN's, and the DWT, PN=0 does
%                not remove ramps
%                [default PN=0]
% 'dwtm',M     = uses images of size M to generate the masters for
%                the DWT responses, if DWTM=0 uses the mirrored &
%                expanded image size (recommended for Meyer WT)
%                [default DWTM=0]
% 'LibDir',dir = reads the master DWT responses from dir
%                [default ~/matlab/SarSeris/DWTImpulseResp]
%
% OUTPUT:
% WaveCoeff    = cell array of all the wavelet coefficients of the
%                igrams, WaveCoeff{k,l}(:,:,n) are the scale k,
%                band l, coefficients for igrams n; k=mstg+1 is
%                only defined for l=1
% Window       = [min_row max_row min_col max_col] the window of
%                the original igram geometry.
% CoefWeights  = cell array of all the weights of the wavelet
%                coefficients based in the filled in decorrelated
%                regions of each igram
%
% See also:
% MintsRemoveSarRamps, MintsTimeString
% Requires:
% inpaint_nans from John D'Errico, http://www.mathworks.com/matlabcentral/fileexchange/4551
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


fcname = 'MintsWT';
Mirror = 0.5;
ComplMask = 1;
meyerdeg = 1;
meyerl = 3;
minstg = 1;
maxstg = 0;
minscl = 1;
maxscl = 4;
type = 'meyer';
DWTdir = '~/matlab/SarSeries/DWTImpulseResp';
DWTM = 0;
PN = 0;
PreFill = 1;
for v=1:length(varargin)
  if ischar(varargin{v})
    switch lower(varargin{v})
     case 'farras'
      type = 'farras';
      fprintf(1,['\n%s: NOTE Farras wavelets are not fully ', ...
		 'verified in current version\n\trecommend using' ...
		 ' Meyer wavelets\n'],fcname);                   
     case 'meyer'
      type = 'meyer';
     case 'meyerdeg'
      meyerdeg = varargin{v+1};
     case 'scale'
      maxscl = varargin{v+1};
      if length(maxscl)==2
	minscl = maxscl(1);
	maxscl = maxscl(2);
      end
     case 'invertmask'
      ComplMask = 1;
      if ~ischar(varargin{v+1})
	ComplMask = varargin{v+1};
      end
     case 'stage'
      maxstg = varargin{v+1};
      if length(maxstg)==2
	minstg = maxstg(1);
	maxstg = maxstg(2);
      end
     case 'mirror'
      Mirror = varargin{v+1};
     case 'prefill'
      PreFill = 1;
      if length(varargin)>v
	if ~ischar(varargin{v+1})
	  PreFill = varargin{v+1};
	end
      end
     case 'ramp'
      PN = varargin{v+1};
     case 'libdir'
      DWTdir = varargin{v+1};
     case {'dwtm','master'}
      DWTM = varargin{v+1};
    end
  end
end


switch type
 case 'farras'
  if maxstg==0
    maxstg = maxscl-1;
  else
    maxscl = maxstg+1;
  end
  [af,sf] = farras;
 case 'meyer'
  % maxstg is the L parameter in the Meyer wavelet transform
  maxstg = 3;
  maxscl = 0;
end

% initialize array of times for checkpoints
% times are
% [start_of_computation;end_of_DWT;end_of_IR;end_of_weight]
timing = zeros(4,6);
timing(1,:) = clock;

DIM = size(A(:,:,1));
B = ones(size(A));
fprintf(1,'%s: transforming %d images using %s wavelets\n',fcname,size(A,3),upper(type));
fprintf(1,['%s: image footprint is %dx%d pixles\n'],fcname,DIM);

% remove ramps from the igrams
if PN>0
  fprintf(1,['%s: removing %d coefficient ramp' ...
	     ' from images\n'],fcname,PN);
  [X,Y] = meshgrid([1:size(A,2)],[1:size(A,1)]);
  A = MintsRemoveSarRamps(X,Y,A,PN);
end

switch type
 case 'farras'
  if sum(DIM./4<[2 2].^(maxscl-1))>0
    fprintf(1,['%s: maximum scale (%d pixels) is greater' ...
	       ' than 1/4 the image dimension\n'],fcname,2^(maxscl-1));
  else
    fprintf(1,['%s: maximum wavelet scale is %d pixels\n'],...
	    fcname,2.^(maxscl-1));
  end
end

%^^^^
% fill the igram holes (decorrelated regions) - DWT requires the
% image is defined everywhere
for n=1:size(A,3)
  a = A(:,:,n);
  b = B(:,:,n);
  nanpos = find(isnan(a(:)));
  b(nanpos) = 0;
  B(:,:,n) = b;
  if PreFill
    if ~isempty(nanpos)
      fprintf(1,['%s: filling NaNs in %0.2f percent of image %d\n'],...
	      fcname,100*length(nanpos)/numel(a),n);
      a = inpaint_nans(a,4);
      A(:,:,n) = a;
    end
    clear a b;
  end
end
% end of hole filling
%vvvv

%^^^^
% mirror the image to avoid edge effects and expand to base 2 number
% of pixels
if Mirror>0
  Ao = A;
  Window = [1 DIM(1) 1 DIM(2)];
  Grab = Window;
  fprintf(1,'%s: mirroring image to avoid edge effects\n', ...
	  fcname);
  OFF = ceil(DIM.*Mirror);
  A = [flipdim(A,2) A flipdim(A,2)];
  A = [flipdim(A,1);A;flipdim(A,1)];
  B = [flipdim(B,2) B flipdim(B,2)];
  B = [flipdim(B,1);B;flipdim(B,1)];
  Window = Window + [DIM(1) DIM(1) DIM(2) DIM(2)] + 1;
  Grab = Grab + [DIM(1) DIM(1) DIM(2) DIM(2)] + [-OFF(1) +OFF(1) -OFF(2) +OFF(2)];

  WIN = 2.^ceil(log2([Grab(2)-Grab(1) Grab(4)-Grab(3)]));
  if strcmp(lower(type),'meyer')
    WIN = [1 1].*(max(WIN));
  end
  ADD = WIN - [Grab(2)-Grab(1) Grab(4)-Grab(3)];
  n1 = round(ADD(1)/2);
  n2 = round(ADD(2)/2);
  Grab = Grab - [n1 n1 n2 n2] + [1 ADD(1) 1 ADD(2)];
  while Grab(1)<1
    Grab(1:2) = Grab(1:2) + size(A,1).*[1 1];
    Window(1:2) = Window(1:2) + size(A,1).*[1 1];
    A = [flipdim(A,1);A];
    B = [flipdim(B,1);B];
  end
  while Grab(2)>size(A,1)
    A = [A;flipdim(A,1)];
    B = [B;flipdim(B,1)];
  end
  while Grab(3)<1
    Grab(3:4) = Grab(3:4) + size(A,2).*[1 1];
    Window(3:4) = Window(3:4) + size(A,2).*[1 1];
    A = [flipdim(A,2),A];
    B = [flipdim(B,2),B];
  end
  while Grab(4)>size(A,2)
    A = [flipdim(A,2),A];
    B = [flipdim(B,2),B];
  end
  A = A(Grab(1):Grab(2),Grab(3):Grab(4),:);
  B = B(Grab(1):Grab(2),Grab(3):Grab(4),:);
  Window = Window - [Grab(1) Grab(1) Grab(3) Grab(3)];
  
  fprintf(1,['%s: mirrored image is expanded to %dx%d pixels\n'], ...
	  fcname,size(A(:,:,1)));
else
  % expand image size to power of 2
  MIR = ceil((2.^ceil(log2(DIM)) - DIM)./DIM);
  if strcmp(lower(type),'meyer')
    MIR = [1 1].*(max(MIR));
  end
  if sum(MIR)>0
    fprintf(1,['%s: expanding image'...
		  ' to base 2 size by mirroring %d/%d times\n'],...
	    fcname,MIR(2),MIR(1));
    for im=1:MIR(1)
      A = [A flipdim(A,2)];
      B = [B flipdim(B,2)];
    end
    for im=1:MIR(2)
      A = [A;flipdim(A,1)];
      B = [B;flipdim(B,1)];
    end
    WIN = 2.^floor(log2(size(A(:,:,1))));
    if strcmp(lower(type),'meyer')
      WIN = [1 1].*(max(WIN));
    end
    A = A(1:WIN(1),1:WIN(2),:);
    B = B(1:WIN(1),1:WIN(2),:);
    fprintf(1,['%s: expanded image is %dx%d pixels\n'], ...
	    fcname,size(A(:,:,1)));
    Window = [1 DIM(1) 1 DIM(2)];
  end
end

% skip this part, the mirrored parts of the igrams need to be
% considered, otherwise part of the total signal is ignored
if 1==0
  % set B outside of the window to be 0 (i.e., parts of the image
  % that are mirrored igram)
  ps = [[1:Window(1)-2],[Window(2)+1:size(A,1)]];
  B(ps,:,:) = 0;
  ps = [[1:Window(3)-2],[Window(4)+1:size(A,2)]];
  B(:,ps,:) = 0;
end
% end of part to skip

% use the compliment of the mask
if ComplMask==1
  fprintf(1,['%s: calculating the wavelet coefficient weights ',...
	     'using the compliment of the igram mask\n'],fcname);
  B = 1-B;
end

% ensure max(B)=1 and min(B)=0
ps = find(B<0.0);
B(ps) = 0.0;
ps = find(B>1.0);
B(ps) = 1.0;

M = size(A(:,:,1));
MM = min(M);
if DWTM==0
  DWTM = max(M);
end
if MM<2.^maxscl
  fprintf(1,['%s: size of maximum scale %d is greater than '...
	     'image geometry\n'],fcname,maxscl);
  for v=1:nargout
    varargout{v} = [];
  end
  return
end

% determine the window of the original image
if nargout>=2
  varargout{2} = Window;
end

% end of mirror/expand
%vvvv

%^^^^
% wavelet transform the images
for n=1:size(A,3)
  a = A(:,:,n);
  b = B(:,:,n);
  if ~PreFill
    fprintf(1,['%s: filling NaN in %0.2f percent of the image\n'],...
	    fcname,100*length(find(isnan(a(:))))/numel(a));
    a = inpaint_nans(a,4);
    A(:,:,n) = a;
  end
  fprintf(1,['%s: determing %s wavelet coefficients' ...
	     ' for image %d\n'],fcname,upper(type),n);
  switch type
   case 'farras'
    f = dwt2D(a,maxstg,af);
   case 'meyer'
    ff = FWT2_YM(a,meyerdeg,maxstg);
    f = DealMeyerCoeffs(ff,maxstg);
    % set maxscl to correspond to the wavelet structures
    % in the farras decomposition
    maxscl = length(f);
  end
    
  if n==1
    % initialze the wavelet coefficient cell on first pass (n=1)
    varargout{1} = cell([maxscl 3]);	
    for k=1:maxscl-1
      for l=1:3
	varargout{1}{k,l} = zeros([size(f{k}{1}) size(A,3)]);
      end
    end
    varargout{1}{maxscl,1} = zeros([size(f{end}) size(A,3)]);
    if nargout>=3
      varargout{3} = varargout{1};
    end
  end

  % deal the wavelet coefficients into the output cell
  for k=minscl:maxscl-1
    for l=1:3 
      varargout{1}{k,l}(:,:,n) = f{k}{l};
    end
  end
  varargout{1}{maxscl,1}(:,:,n) = f{maxscl};
  clear a b;
end
% end of wavelet transform
%vvvv

timing(2,:) = clock;

%^^^^
% determine the impulse DWT response if needed, uses pre-computed
% if saved in the toolbox folder
W = [];
fid = fopen(sprintf('%s/%s_%d_%d_sparse.mat',DWTdir,type,DWTM,maxstg));
if nargout>=3
  if fid==-1 
    irstring = 'construct';
    fprintf(1,['%s: determining the DWT response for individual' ...
	       ' pixels in a %dx%d master image\n'],fcname,DWTM,DWTM);
    W = cell([maxscl 3]);
    for k = 1:maxscl
      scl = min([k maxscl-1]);
      fprintf(1,'\tworking on response at scale %d of %d stages in the DWT - ',...
	      k,maxscl);
      fprintf(1,'of %d rows, working on row',2^scl);
      if k<=maxscl-1
	for l=1:3
	  W{k,l} = zeros([[1 1].*DWTM./2^scl 2^(2*scl)]);
	end
      else
	W{k,1} = zeros([[1 1].*DWTM./2^scl 2^(2*scl)]);
      end
      wind=0;
      for i=1:2^scl
	fprintf(1,' %d',i);
	for j=1:2^scl
	  wind = wind+1;
	  d = zeros(DWTM);
	  d(floor(DWTM/2)+i,floor(DWTM/2)+j) = 1;
	  switch type
	   case 'farras'
	    F = dwt2D(d,maxstg,af);
	   case 'meyer'
	    ff = FWT2_YM(d,meyerdeg,maxstg);
	    F = DealMeyerCoeffs(ff,maxstg);
	  end
	  if k<=maxscl-1
	    for l=1:3
	      W{k,l}(:,:,wind) = F{k}{l};
	    end
	  else
	    W{k,1}(:,:,wind) = F{k};
	  end
	end
      end
      fprintf(1,'\n');
    end
    for k=1:maxscl-1
      for l=1:3
	W{k,l} = W{k,l}([1:M(1)/2^k]+DWTM/2^(k+1)-M(1)/2^(k+1),...
			[1:M(2)/2^k]+DWTM/2^(k+1)-M(2)/2^(k+1),:);
      end
    end
    W{k+1,1} = W{k+1,1}([1:M(1)/2^k]+DWTM/2^(k+1)-M(1)/2^(k+1),...
			[1:M(2)/2^k]+DWTM/2^(k+1)-M(2)/2^(k+1),:);
    % construct the normalization of the impulse response (i.e.,
    % the total impulse response)
    fprintf(1,'%s: calculating the sum of the total response\n',fcname); 
    % re-order into a sparse representation
    [spU,spV] = MintsWTSparseRep(W);
    clear W;
    DN = MintsWTWeightCalc(ones(size(B(:,:,n))),M,minscl,maxscl,spU,spV);
    for k=1:maxscl
      if k<=maxscl-1
	lend = 3;
      else
	lend = 1;
      end
      for l=1:lend
	Dnorm{k,l} = DN{k,l}(1,1);
      end
    end	
    save(sprintf('%s/%s_%d_%d_sparse.mat',DWTdir,type,DWTM,maxstg),'spU','spV','Dnorm');
  else
    irstring = 'load';
    fprintf(1,['%s: using DWT response previously' ...
	       ' computed for a %dx%d master image\n'],fcname,DWTM,DWTM);
    fclose(fid);
    load(sprintf('%s/%s_%d_%d_sparse.mat',DWTdir,type,DWTM,maxstg));
  end
end
% end of computation of the impulse responses
%vvv

timing(3,:) = clock;

%^^^^
% calculate the weight on the wavelet coefficients
for n=1:size(A,3)
  if nargout>=3
    fprintf(1,['%s: constructing coefficient weights based on' ...
	       ' decorrelated pixels in image %d (this is slow' ...
	       ' for low scales and large images!)\n'],fcname,n);
    Dpart = MintsWTWeightCalc(B(:,:,n),M,minscl,maxscl,spU,spV);
    for k=minscl:maxscl
      if k<=maxscl-1
	lend = 3;
      else
	lend = 1;
      end
      for l=1:lend
	varargout{3}{k,l}(:,:,n) = Dpart{k,l}./Dnorm{k,l};
	if ComplMask==1
	  varargout{3}{k,l}(:,:,n) = 1-varargout{3}{k,l}(:,:,n);
	end
      end
    end
    clear Dnorm Dpart;
  elseif nargout>=3
    varargout{3}{maxscl+1,l}(:,:,n) = 1;
  end
end
% end of coeff weight calculation
%vvvv

timing(4,:) = clock;


fprintf(1,'%s: done\n',fcname);
fprintf(1,'\t%s to wavelet transform the images\n',...
	MintsTimeString(etime(timing(2,:),timing(1,:))));
if nargout>=3
  fprintf(1,'\t%s to %s the impulse response\n',...
	  MintsTimeString(etime(timing(3,:),timing(2,:))),irstring);
  fprintf(1,'\t%s to calculate the coefficient weights\n',...
	  MintsTimeString(etime(timing(4,:),timing(3,:))));
  fprintf(1,'\t%s total time\n',...
	  MintsTimeString(etime(timing(4,:),timing(1,:))));
end

return

function D = MintsWTWeightCalc(b,M,minscl,maxscl,spU,spV)
% modified to include bug fix & optimization from PA Agram
[ii,jj]=find(b);
D = cell([maxscl 3]);
for k=minscl:maxscl
  fprintf(1,'\tworking on coefficient weights at scale %d of %d - completed',...
	  k,maxscl);
  scl = min([k maxscl-1]);
  Num = 2^scl;
  Numrow = M(1)/Num;
  Numcol = M(2)/Num;
  if k<=maxscl-1
    lend = 3;
  else
    lend = 1;
  end
  for l=1:lend
    D{k,l} = zeros(M./Num);
  end
  % create a structure element for each resolution
  sBase = 4;
  sFact = 1.5;
  sj = Num*sBase;
  s_e = strel('square',sj*sFact);
  % dilate the mask
  Bdil = imdilate(b,s_e);
  % erode the mask
  Bero = imerode(b,s_e);
  Btst = Bdil + Bero;
  clear Bero Bdil
  % loop over the wavelet coefficients
  ttloop = Numcol;
  ctloop = 0;
  rploop = 0.05;
  iind = (floor((ii-1)/Num)*Num)+1;
  jind = (floor((jj-1)/Num)*Num)+1;
  for j=1:Num:M(2)
    oind(2) = (j - floor(M(2)/2) - 1);
    oind(2) = oind(2)/Num;
    if sum(sum(Btst(:,j:j+Num-1))~=zeros(1,Num))~=0
      oldc = [1:Numcol] - oind(2);
      oldc = oldc(find(oldc>0&oldc<Numcol+1));
      oldccomp = ones(1,Numcol);
      oldccomp(oldc) = 0;
      oldccomp = find(oldccomp);
      newc = [1:Numcol] + oind(2);
      newc = newc(find(newc>0&newc<Numcol+1));
      newccomp = ones(1,Numcol);
      newccomp(newc) = 0;
      newccomp = find(newccomp);
      for i=1:Num:M(1)
	oind(1) = (i - floor(M(2)/2) - 1);
	oind(1) = oind(1)/Num;
	if sum(sum(Btst(i:i+Num-1,j:j+Num-1)))~=0
	  oldr = [1:Numrow] - oind(1);
	  oldr = oldr(find(oldr>0&oldr<Numrow+1));
          oldrcomp = ones(1,Numrow);
          oldrcomp(oldr)=0;
          oldrcomp=find(oldrcomp);
	  newr = [1:Numrow] + oind(1);
	  newr = newr(find(newr>0&newr<Numrow+1));
          newrcomp = ones(1,Numrow);
          newrcomp(newr)=0;
          newrcomp = find(newrcomp);
          sumpos = find((iind==i) & (jind==j));
	  if ~isempty(sumpos)
            sumpos = ii(sumpos)+(jj(sumpos)-1)*M(2);
	    [it,jt] = ind2sub(M,sumpos);
	    sumpos = sub2ind([1 1].*Num,...
			     mod(it-1,Num)+1,...
			     mod(jt-1,Num)+1);
	  end
          G=zeros(size(D{k,l}));
	  for l=1:lend
	    if ~isempty(sumpos)
              G([newr,newrcomp],[newc,newccomp]) = spU{k,l}([oldr,oldrcomp],sumpos)*eye(length(sumpos))* spV{k,l}([oldc,oldccomp],sumpos)';
              D{k,l} = D{k,l} + G; 
	    end
	  end
	end
      end
    end
    ctloop=ctloop+1;
    if ctloop/ttloop>rploop
      rploop = rploop + 0.05;
      fprintf(1,' %d%s',round(ctloop/ttloop*100),char(37));
    end
  end
  fprintf(1,'\n');
end
return




function [f] = DealMeyerCoeffs(ff,L);
% extracts the meyer wavelet coefficients and puts them in a cell
% format used throughout Mints
% L = Meyer DWT parameter

J = size(ff,1);

f = cell([1 log2(size(ff,1))-2]);

cnt=0;
for j=log2(size(ff,1))-1:-1:L
  cnt=cnt+1;
  f{cnt} = cell([1 3]);
  
  r1 = [1:2^j];
  r2 = [2^j+1:2^j+2^j];
  
  f{cnt}{1} = ff(r1,r2);
  f{cnt}{2} = ff(r2,r2);
  f{cnt}{3} = ff(r2,r1);
  
end
f{cnt+1} = ff(r1,r1);

return

function [spU,spV] = MintsWTSparseRep(W);
% re-orders W into a sparse representation
% code from P.S. Agram (piyush@gps.caltech.edu)

M=size(W);
spU = cell(M);
spV = cell(M);
count = 1;
res =zeros(M(1)*M(2)*3,1);
rat = res;
for m=1:M(1),
    if(m==M(1))
         maxn=1;
    else
         maxn = M(2);
    end;
    for n=1:maxn,
         H = abs(W{m,n});
         dim = size(H,3);  
         len = size(H,1);
         spU{m,n}=zeros(len,dim);
         spV{m,n}=zeros(len,dim);
         for k=1:dim,
             [u,s,v]=svd(H(:,:,k));
             len = size(u,1);
             spU{m,n}(:,k) = u(:,1);
             spV{m,n}(:,k) = s(1,1)*v(:,1);
             res(count) = mean(mean(abs(H(:,:,k)-u(:,1)*s(1,1)*v(:,1)')));
             rat(count) = abs(s(2,2)/s(1,1));
             count= count+1;
         end;
    end;
end;

return

%
% update history:
%
% EA Hetland - Caltech Jun 2008
% ehetland@alum.mit.edu
% modified:
% EA Hetland - 04 Dec 2008 (Univ. of Michigan)
%          cleaned up the image mirroring and indexing, added otion
%          to pre-fill NaN's, to save on computation
% EA Hetland - 14 Feb 2009
%          implimented determing wavelet coefficient weights based
%          on filled in NAN's
% EA Hetlabd - 16 Feb 2009
%          changed Mirror option for better memory
% EA Hetland - 19 Feb 2009
%          fixed bug in how the weighting coefficients were
%          calculated that did not allow for non-square matricies
% EA Hetland - 20 Feb 2009
%          updated to new method to determine the wavelet
%          coefficient weights, that accounts for location
%          invariances at each scale, and avoids looping over all
%          pixels in the images
% EA Hetland - 20 Feb 2009
%          changed to read DWT impulse response from disk if they
%          have previously been computed
% CJ DiCaprio - 24 Mar 2009
%          changed order of output variables to get Window first before
%          Coeficient weights
% EA Hetland - 06 Apr 2009
%          fixed bugs in the expansion of un-mirrored igrams to
%          base 2 dimensions, added stage option, and made MaxScale
%          the actual maximum scale in the WT, cleaned up the
%          warning/error messages that are based on igram geometry
% EA Hetland - 06 Apr 2009
%          fixed index bug in wind, non critical
% EA Hetland - 05 Jun 2009
%          cleaned script up
% EA Hetland - 19 Apr 2010
%          added Meyer wavelets, re-ordered computation, and ...
% EA Hetland - 02 Jul 2010
%          fixed bug in the wavelet weight calculation, adapted it
%          to Meyer wavelets, and implimented new method to speed
%          up the calculation
% EA Hetland - 12 Jul 2010
%          fixed bug in wavelet weight calculation
%          affecting memory in computations on large igrams 
% EA Hetland - 14 Jul 2010
%          added in tracking for computation time, cleaned progress
%          tracking notices
% EA Hetland - 29 Jul 2010
%          seperated the mask from the holefilling step, in order to
%          fix bug that was not updating the mask correctly if the
%          PreFill flag was set to 0.
% EA Hetland - 05 Aug 2010
%          incoperated DealmeyerCoeffs function into this file
% EA Hetland - 20 Oct 2010
%          modifications based on debugging by PS Agram (piyush@gps.caltech.edu)
%          + added the sparse representation of the impulse response
%          (MintsWTSparseRep function below) code from PS Agram
%          + folded in a new version of MintsWTWeightCalc from PS Agram:
%             + fixed a bug in MintsWTWeightCalc were some index
%             vectors were incorrectly flipped
%             + modified to use sparse representation of the
%             impulse response
%             + optimized the set operations, and removed calls to
%             intersect and setdiff          
%          + added in options to use the compliment of the mask
%          when computing the wavelet coefficient weights
%          + added in check to make sure mask values are between
%          zero and one
% EA Hetland - 16 Nov 2010
%          replaced cputime with clock for timing functions
% EA Hetland - 18 Jul 2011
%          rewrote windowing section to fix bug in double mirroring
%