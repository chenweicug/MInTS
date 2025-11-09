function [varargout] = MintsConnectMatrix(dates,varargin)
%
%  A = MintsConnectMatrix(dates,[len])
% [A,cycle] = MintsConnectMatrix(dates,[len])
% [A,cycle,start] = MintsConnectMatrix(dates,[len])
% [A,cycle,start,B] = MintsConnectMatrix(dates,[len])
%
% Creates the connectivity matrix (i.e., the transpose of the directed
% incidence matrix) corresponding to a igram stack, and finds all the
% simple cycles of the graph (i.e., all the closed loops, minimum
% cycle length=3). Dates should be with the greater date in the first
% column.
%
% INPUT:
% dates = array of dates of igrams, each line corresponds to a
%         pair, and each column a date
%
% OPTIONS:
% len   = finds cycles of maximum length len
%         [default len = 5];
%
% OUTPUT:
% A     = connectivity matrix
% cycle = dates of all simple cycles (i.e., closed loops)
% start = earliest date of each sub-graph
% B     = connectivity matrix of all simple cycles (i.e.,
%         closed loops)
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


global cycle cyclecnt;

fcname = 'MintsConnectMatrix';
N = 5;
if length(varargin)>=1
  N = varargin{1};
end

M = size(dates,1);
T = unique(dates(:));

pairs = zeros([M 2]);
for m=1:M
  pairs(m,1:2) = [find(dates(m,1)==T) find(dates(m,2)==T)];
end
pairs = fliplr(pairs);


% create the directed connectivity matrix
fprintf(1,'%s: constructing connectivity matrix\n',fcname);
A=zeros(size(pairs,1),max(max(pairs)));
for k=1:size(pairs,1);
  A(k,pairs(k,1)) = -1;
  A(k,pairs(k,2))=+1;
end

varargout{1} = A;

if nargout>=2
  % find all the simple cycles (this is essentially a hack, but it
  % seems to work)
  fprintf(1,'%s: finding all simple cycles up to length %d\n',...
	  fcname,N);
  cycle = [];
  cyclecnt = 0;
  last = 0;
  for n=3:N
    for r = 1:size(A,1)
      clear lp
      lp(1) = find(A(r,:)==-1);
      lp(2) = find(A(r,:)==+1);
      clear pos
      pos = find(A(:,lp(2))==-1);
      TraceConnections(A,lp,pos,2,n);
    end
    numb = length(cycle);
    fprintf(1,'%s: found %d simple cycles of length %d\n',...
	    fcname,numb-last,n);
    last = numb;
  end
  for k=1:length(cycle)
    cycle{k} = T(cycle{k});
  end
  varargout{2} = cycle;
else
  return
end
if nargout>=3
  % another hack to find the earliest time in all subgraphs
  fprintf(1,'%s: finding the starting date to all subgraphs\n',...
	  fcname);
  visited = [];
  notvisited = 1;
  cnt=0;
  while ~isempty(notvisited)
    cnt=cnt+1;
    start(cnt) = notvisited(1);
    visited = [visited start(cnt)];
    row = find(abs(A(:,start(cnt)))==1)';
    visited = sort(FindAllConnected(A,row,visited));
    notvisited = setdiff([1:size(A,2)],visited);
  end
  fprintf(1,'%s: found %d subgraphs\n',...
	  fcname,length(start));
  varargout{3} = start;
else
  return
end
if nargout>=4
  fprintf(1,'%s: creating constraint equations for all simple cycles\n',...
	  fcname);
  B = zeros([length(cycle) size(dates,1)]);
  for k=1:length(cycle)
    pos = zeros([1 length(cycle{k})]);
    for i=1:length(cycle{k})-1
      pos(i) = find(dates(:,1)==cycle{k}(i+1)&...
		    dates(:,2)==cycle{k}(i));
    end
    pos(end) = find(dates(:,1)==cycle{k}(end)&...
		    dates(:,2)==cycle{k}(1));
    B(k,pos(1:end-1)) = -1;
    B(k,pos(end)) = +1;
  end
  varargout{4} = B;
end

return

function visited = FindAllConnected(A,row,visited)

cnt=0;
for r=row
  cnt = cnt+1;
  p = setdiff(find(abs(A(r,:))),visited);
  if ~isempty(p)
    pos(cnt) = p;
    visited = [visited pos(cnt)];
    newrow = setdiff(find(abs(A(:,pos(cnt)))==1)',r);
    visited = FindAllConnected(A,newrow,visited);
  end
end

return

function TraceConnections(A,lp,pos,i,n)

global cycle cyclecnt;

if i==n
  for k=1:length(pos)
    if find(A(pos(k),:)==-1)==lp(1)
      % completed a cycle
      cyclecnt=cyclecnt+1;
	cycle{cyclecnt} = lp;
    end
  end
else
  for k=1:length(pos)
    lp(i+1) = find(A(pos(k),:)==+1);
    if length(lp)==n
      pos2 = setdiff(find(A(:,lp(i+1))==+1),pos(k));
    else
      pos2 = setdiff(find(A(:,lp(i+1))==-1),pos(k));
    end
    TraceConnections(A,lp,pos2,i+1,n);
  end
end

return

%
% update history:
%
% EA Hetland - Caltech Aug 2007
% ehetland@alum.mit.edu
% modified:
% EA Hetland  - 09 Jun 2009 (Univ. of Michigan)
%             cleaned script, included construction of pairs array
%             into script, removed maximum index option, added search
%             for simple cycles
% EA Hetland - 11 Jun 2009
%             added search for disconected subgraphs
%
