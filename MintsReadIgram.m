function [varargout] = MintsReadIgram(file)
%
% [phs] = MintsReadIgram(file_name)
% [X,Y,phs] = MintsReadIgram(file_name)
% [X,Y,phs,mag] = MintsReadIgram(file_name)
%
% Reads a ROI_PAC formatted file, either *.int, *.unw, *.hgt
% the *.rsc file must also be in the same location
%
% INPUT:
% file_name = name of the file to read, including the directory path
%
% OUTPUT:
% phs       = array of phase
% X,Y       = location arrays
% mag       = array of magnitude
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

fcname = 'MintsReadIgram';

fprintf(1,'%s: reading file %s\n',fcname,file);

fid = fopen(sprintf('%s',file));
if fid<0
  fprintf(1,'%s: file %s can not be found\n',fcname,file);
  for v=1:nargout
    varargout{v} = [];
  end
  return
end

rscfile = sprintf('%s.rsc',file);
[w,l,xmin,ymin,xmax,ymax,xo,xs,yo,ys] = GetValuesFromRscFile(rscfile, ...
						  'WIDTH', ...
						  'FILE_LENGTH', ...
						  'XMIN','YMIN', ...
						  'XMAX','YMAX', ...
						  'X_FIRST','X_STEP', ...
						  'Y_FIRST','Y_STEP');

A=fread(fid,[w*2,l],'float32');
fclose(fid);


phs = A([w+1:1:w*2],[1:1:size(A,2)])';
phs(find(phs==0))=NaN;

mag = A([1:1:w],[1:1:size(A,2)])';
mag(find(mag==0))=NaN;


if nargout>=3
  if isempty(xo)
    yvec = [0:1:size(phs,1)-1];
    xvec = [0:1:size(phs,2)-1];
  else
    yvec = [yo:ys:yo+ys*(l-1)];
    xvec = [xo:xs:xo+xs*(w-1)];
  end
  [X,Y]=meshgrid(xvec,yvec);
end

if nargout==1
  varargout{1} = phs;
elseif nargout>=3
  varargout{1} = X;
  varargout{2} = Y;
  varargout{3} = phs;
  if nargout==4
    varargout{4} = mag;
  end
end

return

%
% update history
%
% EA Hetland - Univ. of Michigan Jun 2009
% ehetland@alum.mit.edu
% modified:
%
