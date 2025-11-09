function ind = MintsGetFcnPosition(func,FCN)
% 
% ind = MintsGetFcnPosition(func,FCN)
%
% INPUT:
% func = function name
% FCN  = function structure used by MintsTSeriesFcn.m
%
% OUTPUT:
% ind = index of that function in FCN
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


fcname = 'MintsGetFcnPosition';

p = strfind(FCN.FcnNames,func);
nofunc = 1;
for j=1:length(p);
  if ~isempty(p{j})
    nofunc = 0;
    ind = j;
  end
end
if nofunc
  ind=[];
  if length(varargin)>=1
    disp(sprintf(...
	'%s: can not find "%s" in FCN structure',fcanme,func));
  end
end

return

%
% update history:
%
% EA Hetland - Caltech 2007
% ehetland@alum.mit.edu
% modified:
% EA Hetland - 17 Aug 2010 (Univ. of Michigan)
%
