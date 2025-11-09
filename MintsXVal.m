function [varargout] = MintsXVal(varargin)
%
% m = MintsXVal(G,y,OPTIONS)
% [m,penalty] = MintsXVal(G,y,OPTIONS)
% [m,penalty,svdtruncation] = MintsXVal(G,y,OPTIONS)
% 
% Performs a cross-validation on the penalty parameters of a smoothed
% inversion of Gm=y for underdetermined models. The smoothing can also
% include a smoothing based on the model resolution (the so-called
% "shape smoothing"), so that model coeficients not supported as
% supported by observations are damped more heavily than those with
% better support. The model resultion used for the shape smoothing is
% based on an SVD inversion of G, and the degree of resolution based
% smoothing is effectively controlled by the truncation of the SVD
% pseudoinverse (an integer value, 1 uses only the first singular
% vector, size(G,2) uses all singular vectors). Both the penalty and
% the shape-smoothing parameters can be jointly cross validated, but
% note that the objective cross validation curve does not always have
% a clear minumum for variations in the shape-smoothing
% parameter. If the testing set is empty, then inverts all data
% using only the first penalty parameter and shape smoothing
% parameter in the supplied range.
%
%
% INPUT:
% G = design matrix
% y = vector of observations
%
% OPTIONS: non-case sensitive, all can be bundled in one cell
% 'silent'    = runs without giving output
% 'verbose'   = displays information of use
% 'datacovariance',C = uses the data covariance in matrix C
% 'dataweight',W     = uses the data weight given by matrix X
% 'regularization',F = uses the regularization matrix F (see below)
% 'smoothing',sn = sn=1 damps the model parameters [default]
%                  sn=2 minimizes the gradient of the model parameters
%                  sn=3 minimizes the lapacian of the model parameters
% 'PlotAll',i = if i=1 plots the cross-validation curve and the
%               standard l-curve in a new window, only works for
%               1-D cross validation on the penalty parameter
%               [default is i=0]
% 'noshapesmoothing' = does not add regularization based on the
%               model resolution matrix [default it to add shape
%               smoothing]
% 'shapepower',pwr   = uses shape regularization at (1-MR_ii)^pwr,
%               and effectively controlls the dynamic range of the
%               additional shape smoothing matrix
%               [default pwr=1/2]
% 'SubSets',xval = uses the training/testing subsets in the
%               structure xval for the cross validation
%               xval.train{i} = data indecies to invert in the i'th
%               iteration of the cross-validation
%               xval.test{i} = data indecies to test against in the
%               i'th iteration [default is to do leave-one-out
%               cross validation] If SubSets.test=[], SubSets.train
%               is ignored and inverts all data using the first
%               value in penaltyrange and shaperange.
% 'folds',k   = does k-fold cross validation, ignored if SubSets
%               is passed [default k=length(y), i.e., leave-one-out]
% 'penaltyrange',[minpwr maxpwr num] = searches from 10^minpwr to
%               10^maxpwr in num steps in the cross-validation for
%               the penaly parameter [default =[-5 +5 41]
% 'shaperange',[minp maxp num] = searches the svd trunctaion from
%               minp to maxp, this only affects the construction of
%               the additional shape smoothing matrix and is only
%               used when shape smoothing is turned on, fixed bug
%               to determine the final model using the most optimum
%               shape smoothing matrix
%
% OUTPUT:
% m  = the estimated model
% penalty = the optimum penalty parameter found in the
%           cross-validation
% svdtruncation = the optimum svd trunction in the formation of the
%                 model resolution shape smoothing matrix
% Requires:
% Regularization Tools (regu toolbox; tested with version 4.1),
% from P.C. Hansen http://www2.imm.dtu.dk/~pch/
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
fcname = 'MintsXVal';

for i=1:nargout
  varargout{i} = [];
end

% unwind varargin in case of multiple options passed as a single cell
% in varargin{3}
if length(varargin)==3
  if iscell(varargin{3})
    for k=length(varargin{3}):-1:1
      varargin{3+k-1} = varargin{3}{k};
    end
  end
end

% G is the design matrix
G = varargin{1};
% y is the data vector
y = varargin{2};
% C is the data covariance matrix
C = [];
% W is the data weight matrix
W = eye(length(y));
% F is the smoothing matrix
F = [];

DEBUG_MODE = 0;
DEBUG_MOVIE = 0;

% set defaults, and parse the optional input
Message = 0;
ShapeSmoothing = 1;
rpow = 1/2;
SubSets = [];
XValFolds = length(y);
lamrange = [-8 8 41];
prange = [1 size(G,2) size(G,2)];
Smoothing = 1;
PlotAll = 0;
% PlotAll=1 will plot the ocv and l-curve, you need to hit enter to
% move to the next coefficient...
for k=1:length(varargin)
  if ischar(varargin{k})
    switch lower(varargin{k})
     case 'silent'
      Message=0;
     case 'verbose'
      Message=1;
     case 'plotall'
      PlotAll = varargin{k+1};
     case 'noshapesmoothing'
      ShapeSmoothing = 0;
     case 'regularization'
      F = varargin{k+1};
     case 'dataweight'
      W = varargin{k+1};
     case 'datacovariance'
      C = varargin{k+1};
     case 'shapepower'
      rpow = varargin{k+1};
     case 'subsets'
      SubSets = varargin{k+1};
     case {'folds','groups'}
      XValFolds = varargin{k+1};
     case {'penaltyrange'}
      lamrange = varargin{k+1};
     case {'shaperange','shapesvdrange'}
      prange = varargin{k+1};
     case 'smoothing'
      if ischar(varargin{k+1})
	switch lower(varargin{k+1})
	 case 'damped'
	  Smoothing = 1;
	 case {'gradient','flat'}
	  Smoothing = 2;
	 case {'laplacian','curvature'}
	  Smoothing = 3;
	end
      else
	Smoothing = varargin{k+1};
      end
    end
  end
end

% set up the vector of penalty parameters to try
lamvec = logspace(lamrange(1),lamrange(2),lamrange(3));
if ShapeSmoothing==1
  % set up the vector of svd truncation values to try (used for
  % constructing the initial model resolution & the shape smoothing)
  pvec = round(linspace(prange(1),prange(2),prange(3)));
else
  pvec = 1;
end

if Smoothing==1
  if Message==1
    fprintf(1,'%s: damping the solution\n',fcname);
  end
  F = eye(size(G,2));
elseif Smoothing==2
  % difference smoothing
  if Message==1
    fprintf(1,['%s: regularizing the solution with the gradient of' ...
	       ' the model parameters\n'],fcname);
  end
  for r=1:size(F,2)-1
    F(r,r+1) = -1;
  end
elseif Smoothing==3
  % curvature
  if Message==1
    fprintf(1,['%s: regularizing the solution with the laplacian of' ...
	       ' the model parameters\n'],fcname);
  end
  for r=2:size(F,2)-1
    F(r,r-1) = -1;
    F(r,r+1) = -1;
    F(r,r) = 2;
  end
  F(1,1:2) = [1 -1];
  F(end,end-1:end) = [-1 1];
  F = -F;
else
  fprintf(1,'%s: smoothing choice is not correct, see help\n', ...
	  fcname);
  return
end

if isempty(SubSets)
  % create the training and testing subsets
  if Message==1
    fprintf(1,'%s: using %d-fold cross-validation\n',XValFolds);
  end
  pos = randperm(size(G,1));
  % break up into XValFolds groups (XValFolds = length(y) is
  % leave-one-out cross-validation
  breakpos = round(size(dates,1)/XValFolds);
  subsamples{1} = [1:breakpos];
  for k=2:XValFolds-1
    subsamples{k} = subsamples{k-1}(end)+[1:breakpos];
  end
  subsamples{ngroup} = [subsamples{XValFolds-1}(end)+1:length(T)];
  for k=1:XValFolds
    SubSets.test{k} = sort(pos(subsamples{k}));
    SubSets.train{k} = setdiff([1:size(dates,1)]',SubSets.test{k});
  end
end

if ~isempty(C)
  % find the inverse of the covariance matrix 
  % NOTE: the inverse does not exist for closed Bperp graphs, and
  % so take the psuedoinverse if condition is bad
  [Ci] = MatrixInvert(C);
  % find the composite weighting matrix and its square-root
  Z = W'*Ci*W;
  Zsq = MatrixRoot(Z,'data weight');
else
  Zsq = W;
end
% weight the design matrix and the data vector
Gw = Zsq*G;
dw = Zsq*y;

if DEBUG_MODE==1
  % mcnt is a counter for creating movies of the model at each
  % lambda/p choice, it is used only in the debug mode
  mcnt=0;
end

% calculate the svd of the design matrix for the resolution based smoothing
[U,S,Vr] = svd(G);
clear U S;
if ~isempty(SubSets.test)
  % array to store the sum-of-squares for each test set
  rss = zeros(length(pvec),length(lamvec),length(SubSets.test));
  % start loop over shape smoothing parameters
  for pp=1:length(pvec)
    if ShapeSmoothing==1
      % data resolution using the first pvec(pp) singular vectors
      R = Vr(:,1:pvec(pp))*Vr(:,1:pvec(pp))';
      S = diag(abs(1-diag(R)).^rpow)*F;
      clear R;
    else
      S = F;
    end
    % start loop over training subsets
    for jj=1:length(SubSets.test)
      % the generalized svd for each subset using the given
      % smoothing matrix
      [U,sm,X,V,W] = cgsvd(Gw(SubSets.train{jj},:),S);
      clear V W;
      % start loop over penalty parameters
      for kk=1:length(lamvec)
	% save the model parameters determined from each trial set only
	% for each lambda/p choice
	% do the tikhonov regularization using the cgsvd from above
	[mp,rt,et] = tikhonov(U,sm,X,dw(SubSets.train{jj}),lamvec(kk));
	rss(pp,kk,jj) = sum((dw(SubSets.test{jj}) - ...
			     Gw(SubSets.test{jj},:)*mp).^2);
	%
	%%%% block only used for testing - start
	if DEBUG_MODE==1
	  % set j==the subset number to check
	  if jj==1
	    mcnt=mcnt+1
	    figure(1)
	    clf
	    plot([1:size(G,1)],dw,'ko');
	    hold on
	    plot(SubSets.train{jj},dw(SubSets.train{jj}),'b*','MarkerSize',16)
	    plot(SubSets.test{jj},dw(SubSets.test{jj}),'rp','MarkerSize',16)
	    plot([1:size(G,1)],Gw*mp,'g-')
	    hold off
	    legend('all data','train','test','prediction');
	    
	    title(sprintf('p=%d, lamdba = %0.2e',pvec(pp),lamvec(kk)));
	    if jj==1
	      sax = axis;
	    else
	      axis(sax);
	    end
	    % either generate a movie, or just pause the animation
	    if DEBUG_MOVIE==1
	      AA(mcnt)=getframe(gcf);
	    else
	      pause(0.5);
	    end
	  end
	end
	%%%% block only used for testing - end
	%
      end
      % end loop over penalty parameters
    end
    % end loop over training subsets
  end
  % end loop over shape smoothing parameters
  if DEBUG_MOVIE==1
    movie2avi(AA,'hold_ten_out_xval.avi','fps',5);
    keyboard
  end
  
  % determine the x-val function, find the minimum, and determine the
  % best model using the best penalty parameter
  rss_vec = sum(rss.^2,3);
  [Fmin,pos] = min(rss_vec(:));
  [pCV,lamCV] = ind2sub(size(rss_vec),pos);
  P = pvec(pCV);
  LAM = lamvec(lamCV);
else
  P = pvec(1);
  LAM = lamvec(1);
end

% re-invert with the best P/LAM
if ShapeSmoothing==1
  % data resolution using the first pvec(pp) singular vectors
  R = Vr(:,1:P)*Vr(:,1:P)';
  S = diag(abs(1-diag(R)).^rpow)*F;
  clear R;
else
  S = F;
end
[U,sm,X,V,W] = cgsvd(Gw,S);
if PlotAll==1
  f_vec = zeros(length(lamvec),2);
  for kk=1:length(lamvec)
    [mp,rho,eta] = tikhonov(U,sm,X,dw,lamvec(kk));
    f_vec(kk,1:2) = [rho,eta];
  end
end
[M,rho,eta] = tikhonov(U,sm,X,dw,LAM);


if PlotAll==1
  figure;
  clf
  subplot(1,2,1)
  loglog(lamvec, rss_vec, 'b.'); grid on; hold on;
  loglog(lamvec(lamCV),rss_vec(lamCV),'rp');hold off
  xlabel(' Regularization parameter, log (\lambda)');
  ylabel(' RSS at datapoints NOT used in computing model');
  subplot(1,2,2)
  plot(f_vec(:,1),f_vec(:,2),'o');hold on;
  plot(f_vec(lamCV,1),f_vec(lamCV,2),'rp','MarkerSize',16);hold off
  xlabel('||y - Gm||_2','FontSize',14)
  ylabel('||Sm||_2','FontSize',14)
  grid on
  fprintf(1,'%s: hit enter to continue\n',fcname);
  pause
end 


varargout{1} = M;
if nargout>=2
  varargout{2} = LAM;
end
if nargout>=3
  varargout{3} = P;
end

return

function [Ci] = MatrixInvert(C,varargin)

if length(varargin)>=1
  wlvl = varargin{1};
else
  wlvl = 1e-10;
end

if rcond(C)>1e4*eps
  % invert C directly 
  Ci = inv(C);
else
  % use the psuedoinverse of C 
  [U,S,V] = svd(C);
  spos = find(abs(S)>wlvl);
  Si = zeros(size(S));
  Si(spos) = 1./S(spos);
  Ci = V*Si*U';
end

return

function [Zsq] = MatrixRoot(Z,flag)
global fcname

[T,p]=chol(Z);
if rcond(Z)>1e4*eps & p==0
  Zsq = chol(Z);
else
  Zsq = sqrtm(Z);
  if max(imag(Zsq))>max(real(Zsq))*1e-6
    fprintf(1,['\n%s: WARNING imaginary values on the ',...
  	       'sqrt of the %s inverse, taking the real part only', ...
		' (second warning will appear if root is inexact)\n'],...
	       fcname,flag);
  end
  Zsq = real(Zsq);
end
if max(Zsq'*Zsq-Z)>1e-8
  fprintf(1,['%s: WARNING the approximated sqrt of the %s ',...
	     'matrix is inexact\n'],fcname,flag);
end

return

%
% update history:
%
% EA Hetland - Caltech 2007
% ehetland@alum.mit.edu
% modified:
% EA Hetland - 03 Jun 2009, Univ. Michigan
% EA Hetland - Apr 2010
%             included leave one out ordinary cross validation routine
%             w/ P. Muse
% EA Hetland - 02 Jul 2010
%             embedded ocv functions written by P. Muse into this
%             file, cleaned up the script
% EA Hetland - 16 Aug 2010
%             almost total revision of the script, implemented
%             k-fold cross validation, and fixed numerous
%             bugs in the script
% EA Hetland - 17 Nov 2010
%             re-orderd loops to remove redundencies and minimize
%             computations, cleaned up comments, and added some
%             clear commands to better manage memory
% EA Hetland - 18 Jul 2011
%             added ability to skip the cross validation loop is
%             SubSets.test is empty, just inverts
%