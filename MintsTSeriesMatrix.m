function [varargout] = MintsTSeriesMatrix(T,varargin)
%
% [G,FcnStruct] = MintsTSeriesMatrix(T,varargin)
%
% INPUT:
% T = a vector of times
%
% OPTIONS: 
%    all options can be packed into a cell, strings are case
%    insensitive, can append "no" to any of these to not add them, if
%    not specified, the default does not add them, and can specify
%    'addfunc' (in general '*func*') to add the function
%
% 'Linear'
% 'Step',StepTimes
% 'Log',LogTaus,LogTimes
% 'Exp',ExpTaus,ExpTimes
% 'Seasonal',SeasonalPeriods
% 'USpline',Order,Knots 
%     "       "  ,KnotSpacing
% 'BSpline',Order,Knots
% 'UISpline',Order,KnotSpacing
%      "       "  ,KnotSpacing
%
% see help for MintsTSeriesFcn.m for more information on the
% function arguments
%
% OUTPUT
% G = design matrix
%
% FcnStruct = structure describing the time series decomposition
%      FcnStruct.FcnNames  : cell of time functions
%      FcnStruct.FcnVars   : cell of the varaible list for each
%                          : FCN, used by TSeriesFcn.m
%      FcnStruct.ModelInds : cell of column positions of G associated
%                            with the time function FCN{k}
%
% See also:
% MintsTseriesFcn
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


fcname = 'MintsTSeriesMatrix';

if length(varargin)==1
  if iscell(varargin{1})
    for k=length(varargin{1}):-1:1
      varargin{1+k-1} = varargin{1}{k};
    end
  end
end

% these function names need to be defined in TSeriesFCN.m
FCN{1} = 'linear';
FCNvar{1} = { 1 };
FCN{2} = 'step';
FCNvar{2} = { 1 , 1 };
FCN{3} = 'log';
FCNvar{3} = { 1 , 1 };
FCN{4} = 'exp';
FCNvar{4} = { 1 , 1 };
FCN{5} = 'seasonal';
FCNvar{5} = { 1 };
FCN{6} = 'uspline';
FCNvar{6} = { 0 , 0};
FCN{7} = 'bspline';
FCNvar{7} = { 0 , 0};
FCN{8} = 'uispline';
FCNvar{8} = { 0 , 0};
% this is turned off since it was incorrect
%FCN{9} = 'bipline';
%FCNvar{9} = { 0 , 0};

Find = zeros(1,length(FCN));

for k=1:length(varargin)
  if ischar(varargin{k})
    for j=1:length(FCN)
      if ~isempty(strfind(lower(varargin{k}),FCN{j})) & ...
            isempty(strfind(lower(varargin{k}),'no'))
        Find(j) = 1;
        if (j==2)|(j==5)
          FCNvar{j} = { varargin{k+1} };
        elseif (j>=3)&(j<=4)
          FCNvar{j} = { varargin{k+1:k+2} };
        elseif (j==6)|(j==8)
	  if length(varargin{k+1})==1
	    if mod(T(end)-T(1),varargin{k+2})
	      newspc = (T(end)-T(1))/round((T(end)-T(1))/varargin{k+2});
	      fprintf(1,['%s: uniform B-Spline knot spacing is ',...
                         'not correct, adjusting from %f to %f\n'],...
		      fcname,varargin{k+2},newspc);
	      varargin{k+2} = newspc;
	    end
	    FCNvar{j} = { varargin{k+1} , ...
			  [T(1)-varargin{k+2} ,...
			   [T(1):varargin{k+2}:T(end)] , ...
			   T(end)+varargin{k+2}] };
	    %FCNvar{j} = { varargin{k+1} , ...
		%	  [T(1):varargin{k+2}:T(end)] };
	   
	  else
	    FCNvar{j} = { varargin{k+1} , ...
			  varargin{k+2}};
	  end
        elseif (j==7)|(j==9)
          FCNvar{j} = { varargin{k+1} , ...
                        varargin{k+2}};
        end
      end
    end
  end
end

nunks = [1 length(FCNvar{2}{1}) ...
         length(FCNvar{3}{1}) ...
         length(FCNvar{4}{1}) ...
         2*length(FCNvar{5}{1}) ...
         length(FCNvar{6}{2}) ...
         length(FCNvar{7}{2}) ...
         length(FCNvar{8}{2})];% ... length(FCNvar{9}{2})];


nunk = Find*nunks';
fprintf(1,'%s: creating %d model dimension design matrix\n',fcname,nunk);

% assume T is a column vector
if size(T,1)==1
  T = T';
end

G = zeros(length(T),nunk);
Gcol = cell(1,length(FCN));

stcol = 0;
for k=1:length(FCN)
  if Find(k)
    Gcol{k} = stcol+[1:nunks(k)];
    G(:,Gcol{k}) = MintsTSeriesFcn(T,FCN{k},FCNvar{k});
    stcol = Gcol{k}(end);
  end
end

varargout{1} = G;
if nargout>=2
  D.FcnNames = FCN;
  D.FcnVars = FCNvar;
  D.ModelInds = Gcol;
  varargout{2} = D;
end

return

%
% update history:
%
% EA Hetland - Caltech 2007
% ehetland@alum.mit.edu
% modified:
% EA Hetland - 18 Jun 2009 (Univ. of Michigan)
%          cleaned up script, added integral of B-Splines (uniform
%          and non-uniform), added option to pass either knot
%          spacing or knots for uniform splines, and added check of
%          the knot spacing for uniform splines
% EA Hetland - 02 Jul 2010
%          removed the integral of non-uniform b-spines, cleaned script
%