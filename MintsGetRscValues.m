function  [varargout] = MintsGetRscValues(filename,varargin)
%
% [values] = MintsGetRscValues(file_name,OPTIONS)
%
% Reads the specified values from a ROI_PAC *.rsc file, values are
% returned in order of requested, OPTIONS is a list of rsc key-words
%
% INPUT:
%
% file_name    = name of the file to read, including the directory path
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


if nargout~=length(varargin)
  fprintf(['MintsGetRscValue: expects 1 output variable for every' ...
	   ' request\n']);
  return
end

if isempty(strfind(filename,'.rsc'))
  filename = sprintf('%s.rsc',filename);
end
[namestring,valuestring] = textread(filename,'%s%s');


xo = NaN;
yo = NaN;
for k=1:length(namestring)
  for j=1:length(varargin)
    %switch namestring{k}
    %case varargin{j}
    if strcmp(namestring{k},varargin{j})
      %if strcmp(namestring{k},'Y_STEP')
	%keyboard
      %end
      if isempty(str2num(valuestring{k}));
	varargout{j} = valuestring{k};
      else
	if isempty(regexp(valuestring{k},'\d+-\d+'))
	  varargout{j} = str2num(valuestring{k});
	else
	  valuestring{k}(max(find(valuestring{k}=='-'))) = ' ';
	  varargout{j} = str2num(valuestring{k});
	end
      end
    end
  end
end
if length(varargout)<nargout
  for k=length(varargout)+1:nargout
    varargout{k} = [];
  end
end

return

%
% update history:
%
% EA Hetland - Univ. of Michigan Jun 2009
% ehetland@alum.mit.edu
% modified:
%
