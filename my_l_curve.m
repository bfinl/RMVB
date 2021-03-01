function [reg_corner,rho,eta,reg_param] = my_l_curve(U, sm, b, disp_opt, varargin)

% Set defaults.
smin_ratio = 16*eps;  % Smallest regularization parameter.

% Initialization.
[m,n] = size(U); [p,ps] = size(sm);
beta = U'*b; beta2 = norms(b).^2 - norms(beta).^2;
if (ps==1)
  s = sm; beta = beta(1:p, :);
else
  s = sm(p:-1:1,1)./sm(p:-1:1,2); beta = beta(p:-1:1, :);
end
xi = beta(1:p, :)./repmat(s, [1, size(b, 2)]);
xi( isinf(xi) ) = 0;



%% Main (Tikhonov Regularization)
npoints  = 200;
reg_param  = zeros(npoints, 1);

lambda_min  = max([s(p),s(1)*smin_ratio]);
lambda_max  = s(1);

reg_param(npoints) = lambda_min;
ratio = (lambda_max/reg_param(npoints))^(1/(npoints-1));
for i=npoints-1:-1:1, reg_param(i) = ratio*reg_param(i+1); end

eta = zeros(npoints,1); rho = eta; s2 = s.^2;
for i=1:npoints
    f = repmat(s2./(s2 + reg_param(i)^2), [1, size(b, 2)]);
    eta(i) = sqrt(mean(norms(f.*xi, 2, 1).^2));
    rho_tmp(i, :) = norms((1-f).*beta(1:p, :), 2, 1).';
end
if (m > n) 
    beta2(beta2 < 0) = 0; 
    rho_tmp = sqrt(rho_tmp.^2 + repmat(beta2, [npoints, 1])); 
end
rho  = sqrt(mean(rho_tmp.^2, 2));


%% Locate the "corner" of the L-curve, if required.
if(isempty(varargin))
    reg_opt  = 'Lin';
else
    reg_opt  = varargin{:};
end
switch reg_opt
    case 'Lin'
        [reg_corner, rho_c, eta_c]  = my_l_corner(exp(rho.^2), exp(eta.^2), reg_param);
    case 'Log'
        [reg_corner, rho_c, eta_c]  = my_l_corner(rho, eta, reg_param);
    case 'Closed-Form'
          opt  = 'Lin';
          g = my_lcfun(reg_param,s,beta,xi, opt);
          [~,gi] = min(g);
          reg_corner = fminbnd('lcfun',...
            reg_param(min(gi+1,length(g))),reg_param(max(gi-1,1)),...
            optimset('Display','off'),s,beta,xi); % Minimizer.
          kappa_max = - my_lcfun(reg_corner,s,beta,xi); % Maximum curvature.
          if (kappa_max < 0)
            lr = length(rho);
            reg_corner = reg_param(lr); rho_c = rho(lr); eta_c = eta(lr);
          else
            f = (s.^2)./(s.^2 + reg_corner^2);
            eta_c = norm(f.*xi);
            rho_c = norm((1-f).*beta(1:length(f)));
            if (m>n), rho_c = sqrt(rho_c^2 + norm(b0)^2); end
          end
end


%% Make plot.
[~, ind] = min(abs(reg_corner - reg_param)); 
if( strcmp(disp_opt, 'Display') )
        figure
        plot(rho.^2, eta.^2, 'bo');
        hold on
        plot(rho(ind).^2, eta(ind).^2, 'r*');
        axis([0 max([rho.^2; eta.^2]) 0 max([rho.^2; eta.^2])])
        title([' \lambda_{corner}: ',num2str(reg_corner^2) ', \lambda_{max}: ' num2str(reg_param(1)^2) ', \lambda_{min}: ' num2str(reg_param(end)^2)]);
        xlabel('residual norm || A x - b ||_2^2')
        ylabel('solution norm || x ||_2^2')
        
        figure
        loglog(rho.^2, eta.^2, 'bo');
        hold on
        loglog(rho(ind).^2, eta(ind).^2, 'r*');
        axis([0 max([rho.^2; eta.^2]) 0 max([rho.^2; eta.^2])])
        title([' \lambda_{corner}: ',num2str(reg_corner^2) ', \lambda_{max}: ' num2str(reg_param(1)^2) ', \lambda_{min}: ' num2str(reg_param(end)^2)]);
        xlabel('residual norm || A x - b ||_2^2')
        ylabel('solution norm || x ||_2^2')
end


%%
% f = s2./(s2 + reg_corner^2);
% J = 0;
% for i = 1 : p
%   J  = J + f(i) * xi(i) * VK(:, i);
% end
  

  
end








function [reg_c, rho_c, eta_c]  = my_l_corner(rho, eta, reg_param)


% Set this logical variable to 1 (true) if the corner algorithm
% should always be used, even if the Spline Toolbox is available.
alwayscorner = 0;

% Set threshold for skipping very small singular values in the
% analysis of a discrete L-curve.
s_thr = eps;  % Neglect singular values less than s_thr.

% Set default parameters for treatment of discrete L-curve.
deg   = 2;  % Degree of local smooting polynomial.
q     = 2;  % Half-width of local smoothing interval.
order = 4;  % Order of fitting 2-D spline curve.


  % Use the adaptive pruning algorithm to find the corner, if the
  % Spline Toolbox is not available.
  if ~exist('splines','dir') || alwayscorner
    %error('The Spline Toolbox in not available so l_corner cannot be used')
    reg_c = corner(rho,eta);
    rho_c = rho(reg_c);
    eta_c = eta(reg_c);
    return
  end

  % Otherwise use local smoothing followed by fitting a 2-D spline curve
  % to the smoothed discrete L-curve. Restrict the analysis of the L-curve
  % according to s_thr.
  if (nargin > 3)
    if (nargin==6)       % In case the bound M is in action.
      s = s(index,:);
    end
    index = find(s > s_thr);
    rho = rho(index); eta = eta(index); reg_param = reg_param(index);
  end

  % Convert to logarithms.
  lr = length(rho);
%   lrho = (rho.^2); leta = (eta.^2); slrho = lrho; sleta = leta;
  lrho = log(rho); leta = log(eta); slrho = lrho; sleta = leta;

  % For all interior points k = q+1:length(rho)-q-1 on the discrete
  % L-curve, perform local smoothing with a polynomial of degree deg
  % to the points k-q:k+q.
  v = (-q:q)'; A = zeros(2*q+1,deg+1); A(:,1) = ones(length(v),1);
  for j = 2:deg+1, A(:,j) = A(:,j-1).*v; end
  for k = q+1:lr-q-1
    cr = A\lrho(k+v); slrho(k) = cr(1);
    ce = A\leta(k+v); sleta(k) = ce(1);
  end

  % Fit a 2-D spline curve to the smoothed discrete L-curve.
  sp = spmak((1:lr+order),[slrho';sleta']);
  pp = ppbrk(sp2pp(sp),[4,lr+1]);

  % Extract abscissa and ordinate splines and differentiate them.
  % Compute as many function values as default in spleval.
  P     = spleval(pp);  dpp   = fnder(pp);
  D     = spleval(dpp); ddpp  = fnder(pp,2);
  DD    = spleval(ddpp);
  ppx   = P(1,:);       ppy   = P(2,:);
  dppx  = D(1,:);       dppy  = D(2,:);
  ddppx = DD(1,:);      ddppy = DD(2,:);

  % Compute the corner of the discretized .spline curve via max. curvature.
  % No need to refine this corner, since the final regularization
  % parameter is discrete anyway.
  % Define curvature = 0 where both dppx and dppy are zero.
  k1    = dppx.*ddppy - ddppx.*dppy;
  k2    = (dppx.^2 + dppy.^2).^(1.5);
  I_nz  = find(k2 ~= 0);
  kappa = zeros(1,length(dppx));
  kappa(I_nz) = -k1(I_nz)./k2(I_nz);
  [kmax,ikmax] = max(kappa);
  x_corner = ppx(ikmax); y_corner = ppy(ikmax);

  % Locate the point on the discrete L-curve which is closest to the
  % corner of the spline curve.  Prefer a point below and to the
  % left of the corner.  If the curvature is negative everywhere,
  % then define the leftmost point of the L-curve as the corner.
  if (kmax < 0)
    reg_c = reg_param(lr); rho_c = rho(lr); eta_c = eta(lr);
  else
    index = find(lrho < x_corner & leta < y_corner);
    if ~isempty(index)
      [~,rpi] = min((lrho(index)-x_corner).^2 + (leta(index)-y_corner).^2);
      rpi = index(rpi);
    else
      [~,rpi] = min((lrho-x_corner).^2 + (leta-y_corner).^2);
    end
    reg_c = reg_param(rpi); rho_c = rho(rpi); eta_c = eta(rpi);
  end
  
end





function g = my_lcfun(lambda,s,beta,xi,opt)

% Auxiliary routine for l_corner; computes the NEGATIVE of the curvature.
% Note: lambda may be a vector.  PCH, DTU Compute, Jan. 31, 2015.

% Initialization.
phi = zeros(size(lambda)); dphi = phi; psi = phi; dpsi = phi;
eta = phi; rho = phi;
if length(beta) > length(s)  % A possible least squares residual.
    LS = true;
    rhoLS2 = beta(end)^2;
    beta = beta(1:end-1);
else
    LS = false;
end

% Compute some intermediate quantities.
for i = 1:length(lambda)
  if (nargin==4)
    f  = (s.^2)./(s.^2 + lambda(i)^2);
  else
    f  = s./(s + lambda(i));
  end
  cf = 1 - f;
  eta(i) = norm(f.*xi);
  rho(i) = norm(cf.*beta);
  f1 = -2*f.*cf/lambda(i);
  f2 = -f1.*(3-4*f)/lambda(i);
  phi(i)  = sum(f.*f1.*abs(xi).^2);
  psi(i)  = sum(cf.*f1.*abs(beta).^2);
  dphi(i) = sum((f1.^2 + f.*f2).*abs(xi).^2);
  dpsi(i) = sum((-f1.^2 + cf.*f2).*abs(beta).^2);
end
if LS  % Take care of a possible least squares residual.
    rho = sqrt(rho.^2 + rhoLS2);
end

% Now compute the first and second derivatives of eta and rho
% with respect to lambda;
deta  =  phi./eta;
drho  = -psi./rho;
ddeta =  dphi./eta - deta.*(deta./eta);
ddrho = -dpsi./rho - drho.*(drho./rho);

% Convert to derivatives of log(eta) and log(rho).
dlogeta  = deta./eta;
dlogrho  = drho./rho;
ddlogeta = ddeta./eta - (dlogeta).^2;
ddlogrho = ddrho./rho - (dlogrho).^2;

if( strcmp(opt, 'Log') )
    Let g = curvature.
    g = - (dlogrho.*ddlogeta - ddlogrho.*dlogeta)./...
          (dlogrho.^2 + dlogeta.^2).^(1.5);
elseif( strcmp(opt, 'Lin') )
    g = - (drho.*ddeta - ddrho.*deta)./...
          (drho.^2 + deta.^2).^(1.5);
end
  
  
end









  
