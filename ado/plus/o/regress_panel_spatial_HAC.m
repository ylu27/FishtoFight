function [output] = regress_panel_spatial_HAC(y, X, loc, year, panelvar, locCutoff, lagCutoff, display)

%
% S. HSIANG
% SMH2137@COLUMBIA.EDU
% 5/10
%
% ----------------------------
%
% output = regress_panel_spatial_HAC(y, X, loc, year, panelvar, locCutoff, lagCutoff, display)
%
% Function calculates non-parametric (GMM) spatial and autocorrelation 
% structure using a panel data set.  Spatial correlation is estimated for all
% observations within a given period.  Autocorrelation is estimated for a
% given individual over multiple periods up to some lag length. Var-Covar
% matrix is robust to Heteroskedasticity.
% 
% Arguments: 
% 
% y: dependent variable (Nx1) vector 
% X: independnet variables (Nxk) vector (INCLUDE constant as column)
% loc: (Nx2) vector of locations (Lat, Lon) for each obs in DEGREES
% year: (Nx1) vector of year indicators
% panelvar: (Nx1) vector of panel variables that are common to an individual
% locCutoff: (scalar) cutoff for nonparametric windows in KILOMETERS
% lagCutoff: (scalar) cutoff for number of periods that are autocorrelated
%               {note, Greene recommends at least T^0.25}
% display: (string) 'yes' or 'no', determines if output is displayed
%
% Implementation:
% 
% The kernal used to weight spatial correlations is a uniform kernal that
% discontinously falls from 1 to zero at length locCutoff in all directions.
% This is the kernal recommented by Conley (2008).
%
% (previous versions of this code *did not* use an isotropic kernal).
%
% Observations of the same individual over multiple periods seperated by
% lag L are weighted by 
%
%       w(L) = 1 - L/(lagCutoff+1)
%
% (previous versions that implemented autocorrelation treated it as a third
% dimension, with weights that were multiplied by spatial weights, leading
% to intertemporal correlations across space)
%
% Notes:
%
% Location arguments LOC should specify lat-lon units in DEGREES, however
% locCutoff should be specified in KILOMETERS. 
%
% LocCutoff must exceed zero. CAREFUL: do not supply
% coordinate locations in modulo(360) if observations straddle the
% zero-meridian or in modulo(180) if they straddle the date-line. 
%
% Distances are computed by approximating the planet's surface as a plane
% around each observation.  This allows for large changes in LAT to be
% present in the dataset (it corrects for changes in the length of
% LON-degrees associated with changes in LAT). However, it does not account
% for the local curvature of the surface around a point, so distances will
% be slightly lower than true geodesics. This should not be a concern so
% long as locCutoff is < O(~2000km), probably.
%
% Each time-series for an individual observation in the panel is treated
% with Heteroskedastic and Autocorrelation Standard Errors. If lagCutoff =
% T (the number of sequential observation periods), than this treatment is
% equivelent to the "cluster" command in Stata.
%
% References:
%
%      TG Conley "GMM Estimation with Cross Sectional Dependence" 
%      Journal of Econometrics, Vol. 92 Issue 1(September 1999) 1-45
%      http://www.elsevier.com/homepage/sae/econworld/econbase/econom/frame.htm
%      
%      and 
%
%      Conley "Spatial Econometrics" New Palgrave Dictionary of Economics,
%      2nd Edition, 2008
%
%      and
%
%      Greene, Econometric Analysis, p. 546
%
% (modified from a script written by Ruben Lebowski and Wolfram Schlenker)


if strcmp(display, 'yes')== false
    if strcmp(display, 'no')== false
        disp('SOL: DISPLAY must be a string "YES" or "NO"')
        output = false;  
        return
    end
end


% first-stage OLS
b = X\y;

% get information on number of variables
K = size(X,2);
spatialDim = size(loc,2);    
% set variance-covariance matrix equal to zeros
XeeX = zeros(K);

% SPATIAL CORRELATION ENTRIES IN VAR-COVAR MATRIX
% ------------------------------------------------------------
% loop over different years
yearUnique = sort(unique(year));
Nyears = max(size(yearUnique));
for ti = 1:Nyears
    rows = (year == yearUnique(ti));
    
    % get subset of variables
    X1 = X(rows,:);
    loc1 = loc(rows,:);
    e1 = y(rows) - X1*b;
    
    % get information on number of observations
    N = size(X1,1);
   
    % loop over all observations
    for i = 1:N
        
        % ----------------------------------------------------------------
        % step a: get non-parametric weight
        % Bartlett window of Newey-West
        
        %This is Euclidean distance IN KILOMETERS
        lon_scale = cosd(loc1(i,1))*111; %*ones(size(loc1,1),1);
        lat_scale = 111;
        
        %
        % 1 deg lat = 111km
        % 1 deg lon = 111 km * cos(lat)
        %
        %distance scales lat and lon degrees differently depending on
        %latitude.  The distance here assumes a distortion of Euclidean
        %space that is locally uniform around the location of 'i'
        distance = ((lat_scale*(loc1(i,1)-loc1(:,1))).^2 + (lon_scale*(loc1(i,2)-loc1(:,2))).^2).^0.5;
        
        
        %this zeros terms beyond locCutoff, and weights nearby observations
        %all equally
        %(this is isotropic, unlike an earlier version of this code [SMH])
        window = (locCutoff > distance);
                

%         %---(USE THIS TO IMPLEMENT A LINEAR KERNEL)-----------------------|
%         %this weights contemporaneous observations as a linear function of|
%         %distance                                                         |
%         weight = 1-distance/locCutoff;%                                   |
%         %                                                                 |
%         %this zeros terms beyond locCutoff                                |
%         window = weight.*(weight > 0);%                                   |
%         %-----------------------------------------------------------------|


        % ----------------------------------------------------------------
        % step b: construct X'e'eX for given observation
        XeeXh = repmat(X1(i,:)',1,N).*repmat(e1(i)*e1'.*window',K,1) * X1;
        
        %XeeX = XeeX + (XeeXh + XeeXh')/2;
        %line above is unneccessary, matrix will already be symmetric,
        %adjusted below. [SMH]
        
        XeeX = XeeX + XeeXh; 
    end
end

invXX = inv(X'*X) * size(X,1);

if strcmp(display,'yes')
    disp('================================================================')

    XeeX_temp = XeeX/size(X,1);

    % now get corrected variance-covariance matrix

    V = invXX * XeeX_temp * invXX / size(X,1);
    disp(' ')
    disp('VARIANCE-COV MATRIX WITH CONTEMPOREANEOUS SPATIAL CORR ONLY:')
    disp(' ')
    disp(V) 

end


% AUTOCORRELATION ENTRIES IN VAR-COVAR MATRIX
% ------------------------------------------------------------
% loop over different observations
panelvarUnique = sort(unique(panelvar));
Npanelvar = max(size(panelvarUnique));
for pi = 1:Npanelvar
    
    rows = (panelvar == panelvarUnique(pi));
    
    % get subset of variables
    X1 = X(rows,:);
    year1 = year(rows,:);
    e1 = y(rows) - X1*b;
    
    % get information on number of observations
    N = size(X1,1);
   
    % loop over all observations
    for i = 1:N
        
        % ----------------------------------------------------------------
        % step a: get non-parametric weight
        
        %this is the weight for Newey-West
        weight = 1-abs(year1(i)-year1)/(lagCutoff+1);
        
        %obs var enough apart in time are prescribed to have no estimated
        %correlation (Greene recomments lagCutoff >= T^0.25 {pg 546})
        window = weight.*(abs(year1(i) - year1) <= lagCutoff);
        
        %this is required so diagonal terms in var-covar matrix are not
        %double counted (since they were counted once above for the spatial
        %correlation estimates:
        window = window.*(year1(i) ~= year1);                   
            
        % ----------------------------------------------------------------
        % step b: construct X'e'eX for given observation
        XeeXh = repmat(X1(i,:)',1,N).*repmat(e1(i)*e1'.*window',K,1) * X1;
        
        XeeX = XeeX + XeeXh; 
    end
end


XeeX = XeeX/size(X,1);

% now get corrected variance-covariance matrix

V = invXX * XeeX * invXX / size(X,1);


b = b';
s = diag(V.^0.5)';


n = size(X,1);
k = size(X,2);
p = (1-tcdf(abs(b./s),n-k))*2;


if strcmp(display,'yes')

    disp('-------------------------------------')
    disp(' ')
    disp('VARIANCE-COV MATRIX WITH CONTEMPOREANEOUS SPATIAL CORR AND LOCATION-SPECIFIC AUTOCORRELATION:')
    disp(' ')
    disp(V) 
    disp('-------------------------------------')
    disp(' ')
    disp('COEFFICIENTS')
    disp(' ')
    disp(b)
    disp('STANDARD ERRORS')
    disp(' ')
    disp(s)
    disp('TWO-SIDED P-VAL')
    disp(' ')
    disp(p)
    disp('================================================================')
    
end


parameters = struct('N',Npanelvar,'T',Nyears,'locCutoff',locCutoff,'lagCutoff',lagCutoff,'RunDate', date);
output = struct('coeff',b,'se',s,'p_val',p,'parameters',parameters);

return

