function [varargout] = MintsPlotMatrix(A,varargin)
%
% [plot_handle] = MintsPlotMatrix(A,options)
%
% INPUT:
% A = 2D matrix to plot
%
% options: 
% 'nozero' : does not plot a zeros value
% 'nobar' : does not plat a color bar
% 'nobalance' : does not balance the color scale
% 'blockview' : sets all values to 1, & 'nobar'
% 'noreverse' : use normal y-axis
% 'clabel',clabel : uses clabel to label the colorbar
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


cbar = 1;
ceven = 1;
block = 0;
norev = 0;
clabel = [];
for k=1:length(varargin)
  switch lower(varargin{k})
   case 'nozero'
    A(find(A==0)) = NaN;
   case {'nobar','nocolorbar'}
    cbar = 0;
   case 'nobalance'
    ceven = 0;
   case 'blockview'
    block = 1;
   case 'noreverse'
    norev = 1;
   case 'clabel'
    clabel = varargin{k+1};
  end
end

A(isinf(A)) = NaN;

if block
  A(find(A~=0)) = 1;
  cbar = 0;
  ceven = 1;
end

h = pcolor([[A,zeros(size(A,1),1)];zeros(1,size(A,2)+1)]);
if max(size(A))>50
  shading flat
end
if ceven
  if max(abs(A(:)))==0
    caxis([-1 +1]);
  else
    caxis([-max(max(abs(A))) max(max(abs(A)))]);
  end
end
if cbar
  ch = colorbar;
  if ~isempty(clabel)
    set(get(ch,'YLabel'),'String',clabel);
  end
end

if norev
  set(gca,'XDir','normal','YDir','normal');
else
  set(gca,'XDir','normal','YDir','reverse');
end
set(gca,'XTickLabelMode','manual','YTickLabelMode','manual');

xl = get(gca,'XTickLabel');
yl = get(gca,'YTickLabel');

yl(find(mod(str2num(yl),1)>0),:)=' ';
xl(find(mod(str2num(xl),1)>0),:)=' ';

set(gca,'XTick',get(gca,'XTick')+0.5,'XTickLabel',xl);
set(gca,'YTick',get(gca,'YTick')+0.5,'YTickLabel',yl);

if nargout==1
  vargout{1} = h;
end

return

%
% update history:
%
% EA Hetland - Univ. Michigan 2009
% ehetland@alum.mit.edu
%
