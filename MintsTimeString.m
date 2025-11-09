function [timestr] = MintsTimeString(timenum)
%
% [string] = MintsTimeString(time)
%
% forms a text string of the form "xx units' where xx=the time and
% units = seconds, minutes, hours, days
%
% INPUT:
% time = time in seconds
%
% OUTPUT:
% string       = time string w/ units
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


if timenum<(2*60)
  timestr = 'seconds';
elseif timenum<(2*60*60)
  timenum = timenum/60;
  timestr = 'minutes';
elseif timenum<(2*24*60*60)
  timenum = timenum/3600;
  timestr = 'hours';
else
  timenum = timenum/86400;
  timestr = 'days';
end

timestr = sprintf('%0.3f %s',timenum,timestr);

return

%
% update history:
%
% EA Hetland - Univ. Michigan Jan 2010
% ehetland@alum.mit.edu
%
