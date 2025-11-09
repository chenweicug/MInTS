function [u] = MintsTSeriesFcn(varargin)
%
% U = MintsTSeriesFcn(T,'function',function_options)
%
% INPUT
% T = vector of times to evaulate the function at
%
% The following functions are defined (case insensitive):
%
% Linear functions in time of slope Mags:
%      MintsTSeriesFcn(T,'linear',[Mags]) 
%           U = Mag*T
%
% Step functions centered at Time of magnitude Mag:
%      MintsTSeriesFcn(T,'step',[Times],[Mags])
%           U = 0, T<Time; U = Mag, T>=Time
%
% 1-sided log functions, with time constant Taus, starting at Times,
% with mangitude Mags:
%      MintsTSeriesFcn(T,'log',[Taus],[Times],[Mags])
%           U = Mag*log(1+T/Tau) * TSeriesFcn(T,'step',Time)
%
%
% Step functions centered at Time of magnitude Mag:
%      MintsTSeriesFcn(T,'exp',[Taus],[Times],[Mags]) : 
%           U = -Mag*exp(-T/Tau) * TSeriesFcn(T,'step',Time)
%
% B-Splines of order Order, with uniformally spaced knots, and knot
% coefficients Mags:
%      MintsTSeriesFcn(T,'uspline',Order,Knots,[Mags])
%           U = Sum_i Mag_i * B_Order(t-Knot_i),
%   DelKnot = Knot spacing, the total support of U is [a,b]
%                      a = Knots(1)-DelKnot*(Order+1)/2
%                      b = Knots(end)+DelKnot*(Order+1)/2
%
% B-Splines of order Order, with non-uniformally spaced knots, and
% coefficients Mags:
%      MintsTSeriesFcn(T,'bspline',Order,Knots,[Mags])
%           U = Sum_i Mag_i * B_Order(t-Knot_i),
%   the total support of U is [Knots(1),Knots(end)]
%   note that since the splines are not knot centered, 
%   length(Mag)==length(Knots)-Order-1, and U tapers to zero over
%   the first and last Order intervals
%
% Integral of the B-Splines, with same function options as above:
%      MintsTSeriesFcn(T,'usplinei',...)
%
% Sinusoidal functions of period Periods and magnitude Mags. Each
% period in Periods vector returns two columns/rows in U, one for the
% sin component, and one for the cos component:
%      MintsTSeriesFcn(T,'seasonal',[Periods],[Mags])
%           U = Mag*sin(T*Period), Mag*cos(T*Period)
%
% [option] are optional, defualt is 1 for Mags/Tau, and 0 for Times
% any argument can be a vector, in which case U is a matrix, where
% each column/row (depending on orientation of t) corresponds to each
% element of the [option] vector
%
% OUTPUT:
% U = vector of size(T) corresponding to the time function
%
% See also:
% MintsUniformBSpline, MintsBSpline
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


global fcname;
fcname = 'MintsTSeriesFcn';

t = varargin{1};
u = zeros(size(t));
Transpose = 0;

% expecting t as column vector
if size(t,1)==1
  t = t';
  Transpose = 1;
end

if length(varargin)==1
  varargin{2} = 'linear';
end
if length(varargin)>=3
  if iscell(varargin{3})
    for k=length(varargin{3}):-1:1
      varargin{3+k-1} = varargin{3}{k};
    end
  end
end

switch lower(varargin{2})
 case 'linear'
  fprintf(1,'%s: adding linear rate\n',fcname);
  u = t;
  if length(varargin)>2
    if size(varargin{3},1)~=1
      varargin{3} = varargin{3}';
    end
    u = repmat(u,[1 length(varargin{3})]).*...
	repmat(varargin{3},[length(u) 1]);
  end
 case 'step'
  [Times,Mags] = GetOptionalInput(varargin{3:end});
  [Times,Mags] = CheckTime(Times,Mags);
  [u] = HeavysideFcn(t,Times,Mags);
  %u = (repmat(t,[1 length(Times)])...
  %    >=...
  %    repmat(Times,[length(t) 1])).*...
  %   repmat(Mags,[length(t) 1]);
  for ti=1:length(Times)
    fprintf(1,'%s: adding time step at time %f\n',fcname, ...
	    Times(ti));
  end
 case 'log'
  [Taus,Times,Mags] = GetOptionalInput(varargin{3:end});
  [Taus,Times]   = CheckTau(Taus,Times);
  warning('off','MATLAB:log:logOfZero');
  u = log(1+...
	  (repmat(t,[1 length(Taus)])-...
	   repmat(Times,[length(t) 1]))./...
	  repmat(Taus,[length(t) 1])).*...
      HeavysideFcn(t,Times,Mags);
  %MintsTSeriesFcn(t,'step',Times,Mags);
  u(isnan(u)) = 0.0;
  warning('on','MATLAB:log:logOfZero');
  for ti=1:length(Times)
    fprintf(1,'%s: adding logarithmic function at time %f\n',fcname, ...
	    Times(ti));
  end
 case 'exp'
  [Taus,Times,Mags] = GetOptionalInput(varargin{3:end});
  [Taus,Times]   = CheckTau(Taus,Times);
  u = (-1.*exp(-1.*...
	       (repmat(t,[1 length(Taus)])-...
		repmat(Times,[length(t) 1]))./...
	       repmat(Taus,[length(t) 1]))+1.0).*...
      HeavysideFcn(t,Times,Mags);
  %MintsTSeriesFcn(t,'step',Times,Mags);	
  for ti=1:length(Times)
    fprintf(1,'%s: adding exponential function at time %f\n',fcname, ...
	    Times(ti));
  end
 case 'uspline'
  [Order,Knots,Mags] = GetOptionalInput(varargin{3:end});
  Times = t(1).*ones(size(Knots));
  [Times,Mags] = CheckTime(Times,Mags);
  if CheckKnots(Order,Knots,'uniform')
    u = zeros([length(t) length(Knots)]);
    for k=1:length(Knots);
      u(:,k) = MintsUniformBSpline(Order,diff(Knots(1:2)),t-Knots(k));
    end
    u = u.*HeavysideFcn(t,Times,Mags);%MintsTSeriesFcn(t,'step',Times,Mags);
  end
  fprintf(1,'%s: adding %d uniform B-Splines\n',fcname, ...
	  length(Times));
 case 'uispline'
  [Order,Knots,Mags] = GetOptionalInput(varargin{3:end});
  Times = t(1).*ones(size(Knots));
  [Times,Mags] = CheckTime(Times,Mags);
  if CheckKnots(Order,Knots,'uniform')
    u = zeros([length(t) length(Knots)]);
    for k=1:length(Knots);
      u(:,k) = MintsUniformBSpline(Order,diff(Knots(1:2)),t-Knots(k),1);
    end
    u = u.*HeavysideFcn(t,Times,Mags);%MintsTSeriesFcn(t,'step',Times,Mags);
  end
  fprintf(1,'%s: adding %d integrals of uniform B-Splines\n',fcname, ...
	  length(Times));
 case 'bspline'
  [Order,Knots,Mags] = GetOptionalInput(varargin{3:end});
  Times = t(1).*ones(1,length(Knots)-Order-1);
  [Times,Mags] = CheckTime(Times,Mags);
  if CheckKnots(Order+1,Knots,'nonuniform')
    u = MintsBSpline(Order,Knots,t).*...
	HeavysideFcn(t,Times,Mags);
    %MintsTSeriesFcn(t,'step',Times,Mags);
  end 
  fprintf(1,'%s: adding %d non-uniformly spaced B-Splines\n',fcname, ...
	  length(Times));
 case 'seasonal'
  [Periods,Mags] = GetOptionalInput(varargin{3:end});
  if isempty(Periods)
    Periods = 1;
  end
  u = sin(...
      (repmat(t,[1 2*length(Periods)])./...
       repmat(reshape(...
	   repmat(Periods',[1 2])',...
	   [1 2*length(Periods)]),[length(t) 1])).*(2*pi) + ...
      repmat([0 pi/2],[length(t) length(Periods)])).*...
      HeavysideFcn(t,t(1).*ones(1,length(Periods)*2),Mags);
  %MintsTSeriesFcn(t,'step',t(1).*ones(1,length(Periods)*2),Mags);
end

if Transpose
  u = u';
end

return

function [u] = HeavysideFcn(t,Times,Mags)
  u = (repmat(t,[1 length(Times)])...
       >=...
       repmat(Times,[length(t) 1])).*...
      repmat(Mags,[length(t) 1]);
return  
  

function [Taus,Times] = CheckTau(Taus,Times)
global fcname
if isempty(Taus)
  Taus = 1.0;
end
if size(Taus,1)~=1
  Taus = Taus';
end
if ((length(Times)==1)&(length(Taus)>1))
  Times = Times(1).*ones(size(Taus));
elseif (length(Times)~=length(Taus))
  Times = zeros(size(Taus));
end
return

function [Times,Mags] = CheckTime(Times,Mags)
global fcname
if isempty(Times)
  Times = [0];
end
if isempty(Mags)
  Mags = ones(size(Times));
end
if (length(Mags)<length(Times))&(length(Mags)==1)
  Mags = Mags(1).*ones(size(Times));
elseif (length(Mags)<length(Times))&(length(Mags)>1);
  Mags = [Mags Mags(end).*ones(1,length(Times)-length(Mags))];
else
  Mags = Mags(1:length(Times));
end
return

function [l] = CheckKnots(Order,Knots,flag)
global fcname
l = 1;
if isempty(Order)
  Order = 0;
end
if (length(Knots)-1)<Order
  fprintf(1,['%s: BSpline - '...
	     '(n+1)=%d > %d=number of knots\n'],...
	  fcname,Order+1,length(Knots));
  l = 0;
end
if strcmp(flag,'uniform')&(std(diff(Knots))>1e4*eps)
    fprintf(['%s: Uniform BSpline - '...
	     'knots must be uniformly spaced\n'],fcname);
  l = 0;
end
if strcmp(flag,'nonuniform')&...
      (length(unique(Knots))<length(Knots))
  fprintf(1,['%s: Non-Uniform BSpline - '...
		'knots must not be replicated\n'],fcname);
  l = 0;
end
return

function [varargout] = GetOptionalInput(varargin)
global fcname
varargout = cell(1,nargout);
maxl = 0;
if length(varargin)>nargout
  varargin = varargin{1:nargout};
end
for k=1:length(varargin)
  if size(varargin{k},1)==1
    varargout{k} = varargin{k};
  elseif size(varargin{k},2)==1
    varargout{k} = varargin{k}';
  else
    varargout{k} = reshape(varargin{k}',...
			   [1 prod(size(varargin{k}))]);
  end

end
return

%
% update history
%
% EA Hetland - Caltech 2007
% ehetland@alum.mit.edu
% modified:
% EA Hetland - 18 Jun 2009 (Univ. of Michigan)
%          cleaned up script, added integral of B-Splines (uniform
%          and non-uniform)
% EA Hetland - 02 Jul 1010
%          fixed Bspline integrals, removed non-uniform b-splines
%          (formulation was wrong), fixed bug in CheckKnots,
%          cleaned script, created seperate heavyside function to
%          avoid having to recall MintsTSeriesFcn
%