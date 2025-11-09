function MintsPlotWaveCeoffs(F,varargin)
%
% MintsPlotWaveCoeffs(coefs,OPTIONS)
%
% Plot the wavelet coefficients, in MInTS format for wavelet
% coefficients storgage.
%
% INPUT:
% coefs = wavelet coefficients
%
% OPTIONS:
% 'figure',N      = plots in figure N
% 'image',n       = plots only the n'th image in the coefficient
%                   arrays [default n=1]
% 'minscale',scl1 = plots coefficients starting at scale=scl
% 'maxscale',scl2 = plots coefficients up to scale=scl, note that
%                   scl=MaxScale+1 plots the final scale
%
% See also:
% MintsPlotMatrix
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



fcname = 'MintsPlotWaveCoeffs';
minscl = 1;
maxscl = size(F,1);
figwin = 1;
n = 1;
for v=1:length(varargin)
  if ischar(varargin{v})
    switch lower(varargin{v})
     case 'image'
      n = varargin{v+1};
     case 'figure'
      figwin = varargin{v+1};
     case 'minscale'
      minscl = varargin{v+1};
     case 'maxscale'
      maxscl = varargin{v+1};
    end
  end
end
    
for i=1:size(F,1)
  for j=1:size(F,2);
    if ~isempty(F{i,j})
      F{i,j} = F{i,j}(:,:,n);
    end
  end
end

mhand = FigureLayup(figwin,'orient','land','rows',1,'cols',1,...
		    'xdel',0.0,'ydel',0.0,...
		    'x0',0.1,'y0',0.1,'xf',0.1,'yf',0.1);
fhand(3) = gcf;
ps = get(gcf,'Position');
set(gcf,'Position',[1500 800 ps(3:4)]);
pos = [1 2 4];
cnt = 0;
for k=minscl:maxscl-1
  %if k-K(1)==2
  %  cnt=cnt+1;
  %  mhand = FigureLayup(cnt,'orient','land','rows',1,'cols',1,...
  %		'xdel',0.0,'ydel',0.0,...
  %		'x0',0.1,'y0',0.1,'xf',0.1,'yf',0.1);
  %  fhand(3+cnt) = gcf;
  %  ps = get(gcf,'Position');
  %  set(gcf,'Position',[1500+50*cnt 800-50*cnt ps(3:4)]);
  %end
  
  thand = SubdivideAxes(mhand,'rows',2,'cols',2,...
			'xdel',0,'ydel',0);
  for i=1:2
    for j=1:2
      set(thand(i,j),'XTick',[],'YTick',[])
    end
  end
  for l=1:3
    axes(thand(pos(l)))
    if max(abs(F{k,l}(:)))>0
      G = F{k,l}./max(abs(F{k,l}(:)));
    else
      G = F{k,l};
    end
    MintsPlotMatrix(G,'nobar','balance');
    axis tight
    shading flat
  end
  mhand = thand(1,2);
end
if maxscl==size(F,1)
  axes(mhand);
  if max(abs(F{end}(:)))>0
    G = F{end,1}./max(abs(F{end,1}(:)));
  else
    G = F{end,1};
  end
  MintsPlotMatrix(G,'nobar','balance');
  axis tight
  shading flat;
end

return

%
% update history:
%
% EA Hetland - Univ. Michigan Jan 2009
% ehetland@alum.mit.edu
% modified:
% EA Hetland - 06 Apr 2009
%              fixed index bug, added min/max scale, added help header
% EA Hetland - 05 May 2010
%              added option 'image'
%
