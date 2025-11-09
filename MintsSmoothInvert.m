function [varargout] = MintsSmoothInvert(varargin)
%
% m = MintsSmoothInvert(G,y,FCN,OPTIONS)
% m = MintsSmoothInvert(Gg,yg,Gc,yc,lambda)
% [m,Gg,yg,Gc,yc] = SmoothInvert(...)
%
% inverts G m = y subject to damping/smoothing constraints F_i m = f_i
% i.e. minimizes:
%        ||G m - y||^2 + Sum_i lambda_i ||F_i m - f_i||^2
%
% => m = [(G^t C^{-1} G) + Sum_i lambda_i (F_i^t F_i)]^{-1}*
%                      [G^t C^{-1} y + Sum_i lambda_i F_i^t f_i]
%
%         C is the data covariance operator, default is C = I
% optional output refer to:
% -> m = [Gg + Sum_i lambda_i Gc_i]*[yg + Sum_i lambda_i yc_i]
%
% If you want to re-estimate [m] using new weighting parameters, the
% following can be used to avoid having to redo all of the matrix
% algebra, and assumes that these were output using a previous run of
% SmoothInvert, for instance, . All of these options must be present,
% and all the smoothing options are ignored. There is no checking to
% make sure all the matricies are the right size.
%
% INPUT:
% G    = design matrix
% y    = data vector
% FCN  = structure describing the decomposition of the time
%        function (output of TSeriesMatrix.m)
%
% OPTIONS: non-case sensitive, all can be bundled in one cell
% 'silent'       = to run without printing messages
% 'dataweight',W = uses W as a weighting matrix
%            - types of smoothing -
% 'minimum-function',lambda,[f] = I m = f
% 'flat-function',lambda,[f] = D m = f
%                            D is difference matrix, there must be at
%                            least two parameters in the function for
%                            this constraint to be applied
% 'smooth-function',lambda,t,FCNvar = Del^2 func = 0
%                            Del^2 is the laplacian over the time
%                            vector t, and FCNvar is the cell
%                            containing the fuction options associated
%                            with the function when G was constructed
%                            (FCNvar should be in format as output by
%                            TSeriesMatrix.m)
%
% OUTPUT:
% m  = model vector
% Gg = (G^t C^{-1} G)^{-1}
% yg = G^t C^{-1} y
% Gc(:,:,i) = (F_i^t F_i)^{-1} : 
% yc(:,:,i) = F_i^t f_i
%
% See also:
% MintsGetFunctionPosition
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


fcname = 'MintsSmoothInvert';

for i=1:nargout
  varargout{i} = [];
end


if length(varargin)<3
  fprintf(1,'%s: ERROR input not correct (debug 1)\n',fcname);
  return;
end
  
% unwind varargin in case of multiple options passed as a single cell
% in varargin{4}
if length(varargin)==4
  if iscell(varargin{4})
    for k=length(varargin{4}):-1:1
      varargin{4+k-1} = varargin{4}{k};
    end
  end
end

% set defaults, and parse the optional input
Message = 1;
for k=1:length(varargin)
  if ischar(varargin{k})
    if strcmp(lower(varargin{k}),'message')
      Message = 1;
    end
    if strcmp(lower(varargin{k}),'silent')
      Message = 0;
    end
  end
end

if ~isstruct(varargin{3})
  if length(varargin)>=5
    Gg  = varargin{1};
    yg  = varargin{2};
    Gc  = varargin{3};
    yc  = varargin{4};
    lam = varargin{5};
    if Message
      fprintf(1,['%s: using preformed matricies & lamdba = %0.4eq\n'],...
	      fcname,lam);
    end
  else
    fprintf(1,'%s: ERROR input not correct (debug 2)\n',fcname);
    return
  end
else
  if Message
    printf(1,'%s: parsing input:\n',fcname);
  end
  G = varargin{1};
  y = varargin{2};
  FCN = varargin{3};
  W = eye(length(y));
  cnt=0;
  for k=4:length(varargin)
    if ischar(varargin{k})
      if strcmp(lower(varargin{k}),'silent')|...
	    strcmp(lower(varargin{k}),'message')
      elseif strcmp(lower(varargin{k}),'dataweight')
	if (size(varargin{k+1},1)==size(varargin{k+1},2))&...
	      (size(varargin{k+1},1)==length(y))
	  W = varargin{k+1};
	else
	  fprintf(1,['%s: ERROR the data weighting matrix must be square '...
		     'and of size length(y)\n'],fcanme);
	  return
	end
      else
	func = lower(varargin{k}(strfind(varargin{k},'-')+1:end));
	fpos = MintsGetFunctionPosition(func,FCN,fcname);
	thisr = FCN.ModelInds{fpos};
	if isempty(fpos)
	  return
	end
	if thisr>0
	  cnt=cnt+1;
	  r{cnt} = thisr;
	  F{cnt} = zeros(size(G,2));
	  f{cnt} = zeros(size(F{cnt},2),1);
	  lam(cnt) = varargin{k+1};
	  if ~isempty(strfind(lower(varargin{k}),'minimum'))&r{cnt}>0
	    if Message
	      fprintf(1,['%s: applying a minimum solution',...
			 ' constraint to the %s component (lambda=', ...
			 ' %0.4e)\n'],fcname,func,lam(cnt));
	    end
	    F{cnt}(r{cnt},r{cnt}) = eye(length(r{cnt}));
	    if length(varargin)>=k+2;
	      f{cnt}(r{cnt}) = GetOptionalVector(...
		  varargin{k+2},length(r{cnt}));
	    end
	  elseif ~isempty(strfind(lower(varargin{k}),'flat'))&r{cnt}>0
	    if Message
	      fprintf(1,['%s: applying a solution flatness', ...
			 ' constraint to the %s component (lambda=%0.4e)\n'],...
		      fcname,func,lam(cnt));
	    end
	    if length(r{cnt})>1
	      I=eye(length(r{cnt})-1,length(r{cnt}));
	      I=-I+circshift(I',1)';
	      F{cnt} = PadFMatrixCols(I,r{cnt},size(G,2));
	      f{cnt} = GetOptionalVector(varargin{k+2},size(I,1));
	    else
	      fprintf(1,...
		      '%s: can not apply flatness constraints to "%s"\n',...
		      fcname,func);
	    end
	  elseif ~isempty(strfind(lower(varargin{k}),'first'))&r{cnt}>0
	    if ~ischar(varargin{k+2})
	      if Message
		fprintf(1,['%s: applying a minimum first derivative' ...
			   ' constraint to the %s component over' ...
			   ' interval [%0.2f-%0.2f] (lambda=%0.4e)\n'],...
			fcname,func,varargin{k+2}(1),...
			varargin{k+2}(end),lam(cnt));
	      end
	      F{cnt} = PadFMatrixCols(...
		  gradient(TSeriesFcn(varargin{k+2},func,FCN.FcnVars{fpos})'),...
		  r{cnt},size(G,2));
	      f{cnt} = [];
	    else
	      if Message
		fprintf(1,...
		    '%s: options for %s are wrong, check the help\n',...
			fcname,varargin{k});
	      end
	    end
	  elseif ~isempty(strfind(lower(varargin{k}),'smooth'))&r{cnt}>0
	    if ~ischar(varargin{k+2})
	      if Message
		fprintf(1,['%s: applying a minimum curvature' ...
			   ' constraint to the %s component over' ...
			   ' interval [%0.2f-%0.2f] (lambda=%0.4e)\n'],...
			     fcanem,func,varargin{k+2}(1),...
			     varargin{k+2}(end),lam(cnt));
	      end
	      F{cnt} = PadFMatrixCols(...
		  del2(TSeriesFcn(varargin{k+2},func,FCN.FcnVars{fpos})'),...
		  r{cnt},size(G,2));
	      f{cnt} = [];
	    else
	      if Message
		fprintf(1,...
			'%s: options for %s are wrong, check the help\n',...
			fcname,varargin{k}));
	      end
	    end
	  end
	end
      end
    end
  end
  Gt = G'*W;
  Gg = Gt*G;
  Gc = zeros([size(Gg) length(lam)]);
  yg = Gt*y;
  yc = zeros([size(yg) length(lam)]);
  for k=1:length(F);
    Gc(:,:,k) = F{k}'*F{k};
    if ~isempty(f{k})
      yc(:,:,k) = F{k}'*f{k};
    end
  end
end

if nargout>=2
  varargout{2} = Gg;
end
if nargout>=3
  varargout{3} = yg;
end
if nargout>=4
  varargout{4} = Gc;
end
if nargout>=5
  varargout{5} = yc;
end

Gc = sum(Gc.*repmat(...
    reshape(lam,[1 1 length(lam)]),[size(Gg) 1]),3);
yc = sum(yc.*repmat(...
    reshape(lam,[1 1 length(lam)]),[size(yg) 1]),3);
m = inv(Gg+Gc)*(yg+yc);

varargout{1} = m;


return

function v = GetOptionalVector(option,len)

v = zeros(len,1);
if ~ischar(option)
  if size(option)==[len 1]
    v = option;
  elseif size(option)==[1 len]
    v = option';
  else
    disp(['SmoothInvert: optional vector parameter is wrong size'...
	 ' using default zero vector']);
  end
end
return

function I = PadFMatrixCols(I,r,N)
% pads matrix I with 0 columns to make it a total of N columns, and
% places I in columns given by r
if r(1)>1
  I = [zeros(size(I,1),r(1)-1) I];
end
if r(end)<N
  I = [I zeros(size(I,1),N-r(end))];
end
return

%
% update history
%
% EA Hetland - Caltech Aug 2007
% ehetland@alum.mit.edu
% modified:
% EA Hetland - 18 Nov 2010, Univ. Michigan
%             cleaned up and standardized code
%