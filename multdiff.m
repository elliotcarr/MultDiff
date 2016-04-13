function [u,x] = multdiff(m,kappa,l0,lm,l,u0,Lbnd,Rbnd,tspan,interface,varargin)
% MULTDIFF Solves the one-dimensional multilayer diffusion problem using a
%                 Semi-Analytic method.
%
%   MULTDIFF solves the transient diffusion equation in a one-dimensional 
%   composite slab of finite length consisting of multiple layers. The code
%   is applicable to both perfect and imperfect contact at the interfaces 
%   between adjacent layers and either Dirichlet, Neumann or Robin boundary
%   conditions at the ends of the slab.
%  
%   MULTDIFF is an implementation of the semi-analytic method
%   proposed by Carr and Turner based on the Laplace Transform and an 
%   orthogonal eigenfunction expansion.
%
%   Full details can be found in the paper: 
%   E. J. Carr and I. W. Turner.
%
%   Description:
%   -----------------------------------------------------------------------
%   MULTDIFF solves the standard diffusion equation in each layer 
%   (l(i-1) < x < l(i)):
%
%      du_(i)/dt = d/dx * (kappa(i) * du_(i)/dx),   i = 1,...,m,
%   
%   subject to the following initial and external boundary conditions:
%   
%      u_(i)(x,t) = u0(x)                           at t = 0
%      aL * u_(1)(x,t) + bL * du_(1)/dx(x,t) = cL   at x = l0
%      aR * u_(m)(x,t) + bR * du_(m)/dx(x,t) = cR   at x = lm
%
%   where u_(i) is the solution in layer i, kappa(i) is the diffusivity in 
%   layer i (constant) and aL, bL, cL, aR, bR and cR are constants.
%
%   Either perfect or imperfect contact is imposed at the interfaces 
%   between adajacent layers (at x = l(i), i = 1,...,m-1):
%
%    - Perfect contact 
%       u_(i)(x,t) = u_(i+1)(x,t)                           
%       kappa(i) * u_(i)(x,t) = kappa(i+1) * u_(i+1)(x,t)
%
%    - Imperfect contact
%       kappa(i) * du_(i)/dx(x,t) = H(i) * (u_(i+1)(x,t) - u_(i)(x,t))
%       kappa(i+1) * du_(i+1)/dx(x,t) = H(i) * (u_(i+1)(x,t) - u_(i)(x,t))       
%   
%   Usage:
%   -----------------------------------------------------------------------
%   [U,X] = multdiff(m,kappa,l0,lM,l,u0,bcs,tspan,'Perfect')
%   [U,X] = multdiff(m,kappa,l0,lM,l,u0,bcs,tspan,'Perfect',options)
%   [U,X] = multdiff(m,kappa,l0,lM,l,u0,bcs,tspan,'Imperfect',H)
%   [U,X] = multdiff(m,kappa,l0,lM,l,u0,bcs,tspan,'Imperfect',H,options)
%
%   Input Arguments:
%   -----------------------------------------------------------------------
%   m           Number of layers. Must be an integer greater than or equal 
%               to 3.  
%   kappa       A vector of length M containing the diffusivity values 
%               in each layer such that the diffusivity in Layer i is given
%               by kappa(i) (i = 1,...,m).  
%   l0          x coordinate of the left boundary of the slab
%   lm          x coordinate of the right boundary of the slab
%   l           A vector of length M-1 containing the locations of the
%               interfaces between adjacent layers such that the interface 
%               between Layer i and Layer i+1 is located at L(i) 
%               (i = 1,...,m-1).
%   u0          A function handle specifying the initial condition. The 
%               function uint = u0(X) should accept a vector argument x and 
%               return a vector result uint. Use array operators .*, ./ and 
%               .^ in the definition of u0 so that it can be evaluated with
%               a vector argument.
%   Lbnd        A cell array specifying the boundary condition at the x=l0.
%               Lbnd takes the form {type,aL,bL,cL}, where type is one of
%                    'Dirichlet': aL ~= 0 and bL = 0
%                    'Neumann':   aL = 0 and bL ~= 0
%                    'Robin':     aL ~= 0 and bL ~= 0
%   Rbnd        A cell array specifying the boundary condition at the x=lm.
%               Rbnd takes the form {type,aR,bR,cR}, where type is one of
%                    'Dirichlet': aR ~= 0 and bR = 0
%                    'Neumann':   aR = 0 and bR ~= 0
%                    'Robin':     aR ~= 0 and bR ~= 0
%   tspan       A vector specifying the times at which a solution is 
%               requested. To obtain solutions at specific times 
%               t0,t1,...,tf, use TSPAN = [t0,t1,...,tf].
%   interface   Internal boundary conditions at interfaces between adjacent
%               layers. inteface can be either 'Perfect' or 'Imperfect'.
%   H           A vector of length m-1 containing the contact 
%               transfer coeffecients at the interfaces between adjacent 
%               layers such that the coefficient between layer i and 
%               layer i+1 is given by H(i) (i = 1,..,m-1). 
%               * Applicable to imperfect contant only.
%   options     An (optional) set of solver options. Fields in the 
%               structure options are
%                - N    number of eigenvalues to use in expansions
%                       [N = 50 by default]
%                - NX   number of divisions within each slab. U(:,j) gives
%                       the solution at x = l(i-1):(l(i)-l(i-1))/NX:l(i) 
%                       and t = tspan(j).
%                       [NX = 50 by default]
%                - NZ   number of poles in CF method (see cf.m)  
%                       [NZ = 14 by default]
%                - Hp   value of contact transfer coefficient to 
%                       approximate perfect contact condition
%                       [Hp = 1e6 by default]                      
%
%   Output Arugments:
%   -----------------------------------------------------------------------
%   u   Matrix of solution values. u(:,j) gives the solution on the entire
%       slab (l0 <= x <= lm) at time t = tspan(j) and at the grid points 
%       returned in the output vector x.
%   x   Vector of grid points at which solution is given. x is a vector 
%       taking the form x = [x1,x2,...,xm]', where:
%          x1 = l0:(l(1)-l0)/NX:l(1)
%          xi = l(i-1):(l(i)-l(i-1))/NX:l(i), i = 2,..,m-1
%          xm = l(m-1):(lm-l(m-1))/NX:lm
%
%   Example:
%   -----------------------------------------------------------------------
%   u0 = @(x) zeros(size(x));
%   [u,x] = multdiff(3,[1,0.1,1],0.0,1.0,[0.3,0.7],u0,...
%          {'Dirichlet',1.0,0.0,1.0},{'Dirichlet',1.0,0.0,0.5},...
%          [0.02,0.05,0.1,0.2,0.5,1.0],'Perfect');
%

% Default values
AbsTol = 1e-10;
RelTol = 1e-6;

% AbsTol = 1e-14;
% RelTol = 1e-10;

% -------------------------------------------------------------------------
% Check inputs
% -------------------------------------------------------------------------
if nargin < 10
    error('Not enough input arguments.');
elseif nargin == 10
    if strcmp(interface,'Imperfect')
        error('H must be specified for imperfect contact at interfaces.');
    end
    options = struct;    
elseif nargin == 11
    if strcmp(interface,'Perfect')
        options = varargin{1};
    elseif strcmp(interface,'Imperfect')
        H = varargin{1};
        options = struct;
    end
elseif nargin == 12
    if strcmp(interface,'Perfect')
        error('Too many input arguments for interface = ''Perfect''.');
    elseif strcmp(interface,'Imperfect')
        H = varargin{1};
        options = varargin{2};
    end    
else
    error('Too many input arguments.');
end

% Number of layers
if round(m) ~= m || m < 3
    error('m must be an integer greater than or equal to 3.')
end

% Diffusivities
if length(kappa) ~= m || sum(kappa > 0) ~= m
    error('kappa must be a vector of length m with kappa(i)>0.')
end

% Slab left and right boundary
if l0 > lm
    error('l0 must be less than lm.')
end

% Interfaces
if length(l) ~= m-1 || sum(diff(l)>0) ~= m-2
    error('l must be a vector of length m-1 with with increasing values.')
elseif l(1) <= l0 || l(m-1) >= lm
    error('l(1) must be greater than l0 and l(m-1) must be less than lm.');
end

% Initial condition
if ~isa(u0,'function_handle') || nargin(u0) ~= 1
    error('u0 must be a function handle of the form uint = u0(x).');
end

% Boundary conditions
if ~isa(Lbnd,'cell') || length(Lbnd) ~= 4
    error(['Lbnd must be a cell array of length 4.']);
end
if ~isa(Rbnd,'cell') || length(Rbnd) ~= 4
    error(['Rbnd must be a cell array of length 4.']);
end

% Time vector
tlength = length(tspan);
if sum(tspan > 0) ~= tlength
    error('tspan must have entries that are greater than or equal to 0.')
end

% Internal boundary conditions at interfaces
if strcmp(interface,'Perfect') || strcmp(interface,'Imperfect')
else
    error('interface must be either ''Perfect'' or ''Imperfect''.')
end

% Check options structure
if ~isa(options,'struct')
    error('options must be a structure.')
end
Names = {'N','NX','NZ','Hp'};
fn = fieldnames(options);
for i = 1:length(fn)
    j = strcmp(fn(i),Names);
    if sum(j) == 0
        error('Invalid option ''%s''.',fn{i});
    end
end
% Number of eigenvalues to use in expansions
if isfield(options,'N')
    N = options.N;
    if round(N) ~= N && N < 1
        error('options.N must be an integer greater than or equal to 1.')
    end
else
    N = 50; % Default
end
% Number of divisions within each slab
if isfield(options,'NX')
    NX = options.NX; 
    if round(NX) ~= NX && NX < 1
        error('options.NX must be an integer greater than or equal to 1.')
    end
else
    NX = 50; % Default
end
% Number of poles in CF method (see cf.m)  
if isfield(options,'NZ'), 
    NZ = options.NZ;
    if round(NZ) ~= NZ && NX < 1
        error('options.NZ must be an integer greater than or equal to 1.')
    end
else
    NZ = 14; % Default
end
% Value of contact transfer coefficient to approximate perfect contact 
% condition
if isfield(options,'Hp')
    if strcmp(interface,'Perfect')
        Hp = options.Hp;
        if Hp < 0
            error('options.Hp must be greater than or equal to 0.')
        end
    else
        warning('options.Hp is specified but not used.')
    end
else
    Hp = 1.0e6; % Default
end
if strcmp(interface,'Perfect')
    H = Hp*ones(m-1,m);
end

% Get boundary condition constants
Ltype = Lbnd{1};
Rtype = Rbnd{1};
aL    = Lbnd{2};
bL    = Lbnd{3};
cL    = Lbnd{4};
aR    = Rbnd{2};
bR    = Rbnd{3};
cR    = Rbnd{4};

% Check boundary conditions are implemented correctly
if sum(strcmp(Ltype,{'Dirichlet','Neumann','Robin'})) == 0
    error(['Boundary condition at left boundary must be one of either',...
        ' ''Dirichlet'', ''Neumann'' or ''Robin''.']);
end
if sum(strcmp(Rtype,{'Dirichlet','Neumann','Robin'})) == 0
    error(['Boundary condition at right boundary must be one of either',...
        ' ''Dirichlet'', ''Neumann'' or ''Robin''.']);
end
if aL == 0 && strcmp(Ltype,'Dirichlet')
    error('Dirichlet condition at left boundary cannot have aL = 0.');
end
if bL == 0 && strcmp(Ltype,'Neumann')
    error('Neumann condition at left boundary cannot have bL = 0.');
end
if aR == 0 && strcmp(Rtype,'Dirichlet')
    error('Dirichlet condition at right boundary cannot have aR = 0.');
end
if bR == 0 && strcmp(Rtype,'Neumann')
    error('Neumann condtion at right boundary cannot have bR = 0.');
end
if strcmp(Ltype,'Dirichlet') && bL ~= 0
    error('Dirichlet condition at left boundary cannot have bL = 0.');
end
if strcmp(Rtype,'Dirichlet') && bR ~= 0
    error('Dirichlet condition at right boundary cannot have bR = 0.');
end
if strcmp(Ltype,'Neumann') && aL ~= 0
    error('Neumann condition at left boundary cannot have aL = 0.');
end
if strcmp(Rtype,'Neumann') && aR ~= 0
    error('Neumann condition at right boundary cannot have aR = 0.');
end
if (aL == 0 || bL == 0) && strcmp(Ltype,'Robin')
    error('Robin condition at left boundary cannot have aL = 0 or bL = 0');
end
if (aR == 0 || bR == 0) && strcmp(Rtype,'Robin')
    error('Robin condition at left boundary cannot have aR = 0 or bR = 0');
end
if strcmp(Ltype,'Robin') && aL/bL < 0
    warning('Robin condition at left boundary must have aL/bL > 0.');
end
if strcmp(Rtype,'Robin') && aR/bR < 0
    warning('Robin condition at right boundary must have aR/bR > 0.');
end
if aL == 0 && bL == 0
    error('Boundary condition is incorrect at left boundary (aL = bL = 0).')
end
if aR == 0 && bR == 0
    error('Boundary condition is incorrect at left boundary (aR = bR = 0).')
end

if strcmp(Ltype,'Neumann') && strcmp(Rtype,'Neumann') && ...
        (kappa(1)*bR*cL + kappa(m)*bL*cR) ~= 0
    error(['If Neumann boundary conditions are applied at both ends then ', ...
        'kappa(1)*bR*cL + kappa(m)*bL*cR must be equal to zero.'])
end
    

% -------------------------------------------------------------------------
% Compute function w(x) that satisfies non-homogeneous BCs
%
% Slab 1:
% w1(x) = w(1) + w(2)*x            for Dirichlet, Robin
% w1(x) = w(1)*x + w(2)*x^2        for Neumann
%
% Slab i = 2,..,m-1:
% wi(x) = w(2*i-1) + w(2*i)*x      for Dirichlet, Robin, Neumann
%
% Slab m:
% wm(x) = w(2*m-1) + w(2*m)*x      for Dirichlet, Robin
% wm(x) = w(2*m-1)*x + w(2*m)*x^2  for Neumann
%
% So we want to solve the linear system of equations:
%
%  kappa(1)*w_(1)'(l(1)) = H(1)*(w_(2)(l(1)) - w_(1)(l(1))
%  kappa(1)*w_(1)'(l(1)) = kappa(2)w_(2)'(l(1))
%                     .
%                     .
%                     .
%  kappa(m-1)*w_(m-1)'(l(m-1)) = H(m-1)*(w_(m)'(l(m-1))-w_(m-1)'(l(m-1)))
%  kappa(m-1)*w_(m-1)'(l(m-1)) = kappa(m) * w_(m)'(l(m-1))
%  aL * w_(1)(l0) + bL * w_(1)'(l0) = cL
%  aR * w_(m)(lm) + bR * w_(m)'(lm) = cR
% -------------------------------------------------------------------------

A = zeros(2*m,2*m);
b = zeros(2*m,1);

% if strcmp(Ltype,'Neumann') && strcmp(Rtype,'Neumann')
%     
%     % Interface conditions
%     for i = 1:m-1
%         A(2*i-1,2*i-1) = l(i)+1/H(i)*kappa(i);
%         A(2*i-1,2*i)   = 2*kappa(i)*l(i)/H(i)+l(i)^2;
%         A(2*i-1,2*i+1) = -l(i);
%         A(2*i-1,2*i+2) = -l(i)^2;
%         b(2*i-1) = 0;
%         
%         A(2*i,2*i-1) = kappa(i);
%         A(2*i,2*i) = 2*kappa(i)*l(i);
%         A(2*i,2*i+1) = -kappa(i+1);
%         A(2*i,2*i+2) = -2*kappa(i+1)*l(i);
%         b(2*i) = 0;
%     end
%     
%     % Right boundary condition
%     A(2*m-1,2*m)   = 2*bR*lm;
%     A(2*m-1,2*m-1) = bR;
%     b(2*m-1)       = cR;
%     
%     % Left boundary condition
%     A(2*m,1) = bL;
%     b(2*m)   = cL;
    

% else
    
    % Interface conditions
    for i = 1:m-1
        A(2*i-1,2*i-1) = 1.0;
        A(2*i-1,2*i)   = (kappa(i)/H(i)) + l(i);
        A(2*i-1,2*i+1) = -1.0;
        A(2*i-1,2*i+2) = -l(i);
        b(2*i-1)       = 0.0;
        
        A(2*i,2*i)   = kappa(i);
        A(2*i,2*i+2) = -kappa(i+1);
        b(2*i)       = 0.0;
    end
    
    % Right boundary condition
    A(2*m-1,2*m)   = bR+aR*lm;
    A(2*m-1,2*m-1) = aR;
    b(2*m-1)       = cR;
    
    % Left boundary condition
    A(2*m,1) = aL;
    A(2*m,2) = -bL+aL*l0;
    b(2*m)   = cL;
    
% end

if strcmp(Ltype,'Neumann') && strcmp(Rtype,'Neumann')
    i = 1;
    A(2*m+1,2*i-1) = l(1) - l0;
    A(2*m+1,2*i)   = (l(1)^2-l0^2)/2;
    b(2*m+1)       = integral(u0,l0,l(1),'AbsTol',AbsTol,'RelTol',RelTol);
    for i = 2:m-1
        A(2*m+1,2*i-1) = l(i)-l(i-1);
        A(2*m+1,2*i)   = (l(i)^2-l(i-1)^2)/2;
        b(2*m+1)       = b(2*m+1) + integral(u0,l(i-1),l(i),'AbsTol',AbsTol,'RelTol',RelTol);
    end
    i = m;
    A(2*m+1,2*i-1) = lm-l(m-1);
    A(2*m+1,2*i)   = (lm^2-l(m-1)^2)/2;
    b(2*m+1)       = b(2*m+1) + integral(u0,l(m-1),lm,'AbsTol',AbsTol,'RelTol',RelTol);
end


w = A\b;
% pause;
% w = 0.5*ones(size(w))
% pause;

% -------------------------------------------------------------------------
% Eigenvalues
% -------------------------------------------------------------------------
eigs = zeros(N,m);

% Slab 1 (First slab)
switch Ltype
    case 'Dirichlet'
        eigs(:,1) = (2*[0:N-1]+1)'*pi/(2*(l(1)-l0));
    case 'Neumann'
        eigs(:,1) = [0:N-1]'*pi/(l(1)-l0);
    case 'Robin'
        eigs(:,1) = eigvals(l(1),l0,aL/bL,N);
end

% Slab 2,..., m-1 (Middle slabs)
for i = 2:m-1
    eigs(:,i) = [0:N-1]'*pi/(l(i)-l(i-1));
end

% Slab m (End slab)
switch Rtype
    case 'Dirichlet'
        eigs(:,m) = (2*[0:N-1]+1)'*pi/(2*(lm-l(m-1)));
    case 'Neumann'
        eigs(:,m) = [0:N-1]'*pi/(lm-l(m-1));
    case 'Robin'
        eigs(:,m) = eigvals(lm,l(m-1),aR/bR,N);
end

% -------------------------------------------------------------------------
% Eigenfunction normalisation constants
% -------------------------------------------------------------------------
eigs_norm = zeros(N,m);

% Slab 1 (First slab)
for n = 1:N
    
    if strcmp(Ltype,'Dirichlet')
        eigs_norm(n,1) = sqrt(2/(l(1)-l0));
    elseif strcmp(Ltype,'Neumann')
        if n == 1
            eigs_norm(n,1) = sqrt(1/(l(1)-l0));
        else
            eigs_norm(n,1) = sqrt(2/(l(1)-l0));
        end
    elseif strcmp(Ltype,'Robin')
        lambda = eigs(n,1);
        eigs_norm(n,1) = sqrt((2*(1+(bL^2/aL^2)*lambda^2))/...
            ((bL/aL)+(l(1)-l0)*(1+(bL^2/aL^2)*lambda^2)));
    end
    
end

% Slab 2,..., m-1 (Middle slabs)
for i = 2:m-1
    eigs_norm(1,i) = sqrt(1/(l(i)-l(i-1)));
    for n = 2:N
        eigs_norm(n,i) = sqrt(2/(l(i)-l(i-1)));
    end
end

% Slab m (End slab)
for n = 1:N
    
    if strcmp(Rtype,'Dirichlet')
        eigs_norm(n,m) = sqrt(2/(lm-l(m-1)));
    elseif strcmp(Rtype,'Neumann')
        if n == 1
            eigs_norm(n,m) = sqrt(1/(lm-l(m-1)));
        else
            eigs_norm(n,m) = sqrt(2/(lm-l(m-1)));
        end
    elseif strcmp(Rtype,'Robin')
        lambda = eigs(n,m);
        eigs_norm(n,m) = sqrt((2*(1+(bR^2/aR^2)*lambda^2))/...
            ((bR/aR)+(lm-l(m-1))*(1+(bR^2/aR^2)*lambda^2)));
    end
    
end

% eigs_norm

% -------------------------------------------------------------------------
% Grid spacing within each slab
% -------------------------------------------------------------------------
xgrid = zeros(NX+1,m);

% Slab 1 (First slab)
xgrid(:,1) = l0:(l(1)-l0)/NX:l(1);

% Slabs 2,...,m
for i = 2:m-1
    xgrid(:,i) = l(i-1):(l(i)-l(i-1))/NX:l(i);
end

% Slab m (First slab)
xgrid(:,m) = l(m-1):(lm-l(m-1))/NX:lm;

% -------------------------------------------------------------------------
% Initial conditions - expand in terms of eigenfunctions
% -------------------------------------------------------------------------
c = zeros(N,m); % Coefficients

% Slab 1 (First slab)
for n = 1:N
    lambda = eigs(n,1);
    prod = @(x) (u0(x)-wfunc(1,Ltype,Rtype,x,m,w)) .* ...
        eigfunc(lambda,1,Ltype,Rtype,x,m,l0,lm,l);
    c(n,1) = eigs_norm(n,1)*integral(prod,l0,l(1),'AbsTol',AbsTol,'RelTol',RelTol);
end

% Slabs 2,...,m
for i = 2:m-1
    for n = 1:N
        lambda = eigs(n,i);
        prod = @(x) (u0(x)-wfunc(i,Ltype,Rtype,x,m,w)) .* ...
            eigfunc(lambda,i,Ltype,Rtype,x,m,l0,lm,l);
        c(n,i) = eigs_norm(n,i)*integral(prod,l(i-1),l(i),'AbsTol',AbsTol,'RelTol',RelTol);
    end
end

% Slab m
for n = 1:N
    lambda = eigs(n,m);
    prod = @(x) (u0(x)-wfunc(m,Ltype,Rtype,x,m,w)) .* ...
        eigfunc(lambda,m,Ltype,Rtype,x,m,l0,lm,l);
    c(n,m) = eigs_norm(n,m)*integral(prod,l(m-1),lm,'AbsTol',AbsTol,'RelTol',RelTol);
end
% wfunc(1,Ltype,Rtype,x,m,w)
% wfunc(m,Ltype,Rtype,x,m,w)
% c

% Plot initial condition
u = zeros(NX+1,m);
for n = 1:N
    for i = 1:m
        lambda = eigs(n,i);
        u(:,i) = u(:,i)+c(n,i)*...
            eigs_norm(n,i)*...
            eigfunc(lambda,i,Ltype,Rtype,xgrid(:,i),m,l0,lm,l);
    end
end

usoln = reshape(u,(NX+1)*m,1);
x = reshape(xgrid,(NX+1)*m,1);

% % % figure;
% % % % for i = 1:m
% % % %     uint(:,i) = u0(xgrid(:,i))-wfunc(i,Ltype,Rtype,xgrid(:,i),m,w);
% % % % end
% % % % for i = 1:m
% % % %     u(:,i) = u(:,i) + wfunc(i,Ltype,Rtype,xgrid(:,i),m,w);
% % % % end
% % % usoln = reshape(u,(NX+1)*m,1);
% % % plot(x,usoln,'r')
% % % % hold on
% % % %plot(x,reshape(uint,(NX+1)*m,1),'b')
% % % drawnow
% % % hold off
% % % % pause;

% eigs
% pause;

% norm(u(end,1)-u(1,2),inf)
% pause;

% Get weights and poles for use in inverse transform
[zk,ck] = cf(NZ);

usoln = zeros((NX+1)*m,tlength);

% Time loop
for j = 1:tlength
    
    t = tspan(j);
    
    % Solution (at given time)
    u = zeros(NX+1,m);
    
    for i = 1:m
        u(:,i) = wfunc(i,Ltype,Rtype,xgrid(:,i),m,w);
    end
    
%     if strcmp(Ltype,'Neumann') && strcmp(Rtype,'Neumann')
%         u(:,1) = u(:,1)+2*kappa(1)*w(2)*sqrt(l(1)-l0)*t;
%         for i = 2:m-1
%             u(:,i) = u(:,i)+2*kappa(i)*w(2*i)*sqrt(l(i)-l(i-1))*t;
%         end
%         u(:,m) = u(:,m)+2*kappa(m)*w(2*m)*sqrt(lm-l(m-1))*t;
%     end
    
    % Compute inverse Laplace transform of interface functions
    vbar = zeros(m-1,NZ/2);
    
    for k = 1:NZ/2
        
        A    = zeros(m-1,m-1);
        vr   = zeros(m-1,1);
        b    = zeros(m-1,1);
        
        poles = 2*k-1;
        s = zk(poles)/t;
        
%         if strcmp(Ltype,'Neumann') && strcmp(Rtype,'Neumann')
%             b(1) =  2*kappa(2)*w(4)*sqrt(l(2)-l(1))/s^2 ...
%                 -2*kappa(1)*w(2)*sqrt(l(1)-l0)/s^2;
%             for i = 2:m-2
%                 b(i) = 2*kappa(i+1)*w(2*(i+1))*sqrt(l(i+1)-l(i))/s^2 ...
%                    -2*kappa(i)*w(2*i)*sqrt(l(i)-l(i-1))/s^2;
%             end
%             b(m-1) = 2*kappa(m)*w(2*m)*sqrt(lm-l(m-1))/s^2 ...
%                -2*kappa(m-1)*w(2*(m-1))*sqrt(l(m-1)-l(m-2))/s^2;
%         end
        
        for n = 1:N
            
            % Interface 1 (between layers 1 and 2)
            lambda = eigs(n,1);
            phin_r = eigfunc(lambda,1,Ltype,Rtype,l(1),m,l0,lm,l);
            s1     = s+kappa(1)*lambda^2;
            A(1,1) = A(1,1)+eigs_norm(n,1)^2*phin_r*phin_r/s1;
            b(1)   = b(1)-c(n,1)*eigs_norm(n,1)*phin_r/s1;
            
            lambda = eigs(n,2);
            phin_l = eigfunc(lambda,2,Ltype,Rtype,l(1),m,l0,lm,l);
            phin_r = eigfunc(lambda,2,Ltype,Rtype,l(2),m,l0,lm,l);
            s2     = s+kappa(2)*lambda^2;
            A(1,1) = A(1,1)+eigs_norm(n,2)^2*phin_l*phin_l/s2;
            A(1,2) = A(1,2)-eigs_norm(n,2)^2*phin_l*phin_r/s2;
            b(1)   = b(1)+c(n,2)*eigs_norm(n,2)*phin_l/s2;
            
            % Middle interfaces
            for i = 2:m-2
                lambda   = eigs(n,i);
                phin_l   = cos(lambda*(l(i)-l(i-1)));
                phin_r   = 1.0;
                s1       = s+kappa(i)*lambda^2;
                A(i,i-1) = A(i,i-1)-eigs_norm(n,i)^2*phin_l*phin_r/s1;
                A(i,i)   = A(i,i)+eigs_norm(n,i)^2*phin_r*phin_r/s1;
                b(i)     = b(i)-c(n,i)*eigs_norm(n,i)*phin_r/s1;
                
                lambda   = eigs(n,i+1);
                phin_l   = cos(lambda*(l(i+1)-l(i)));
                phin_r   = 1.0;
                s2       = s+kappa(i+1)*lambda^2;
                A(i,i)   = A(i,i)+eigs_norm(n,i+1)^2*phin_l*phin_l/s2;
                A(i,i+1) = A(i,i+1)-eigs_norm(n,i+1)^2*phin_r*phin_l/s2;
                b(i)     = b(i)+c(n,i+1)*eigs_norm(n,i+1)*phin_l/s2;
            end
            
            % Interface m (between layers m and m-1)
            lambda     = eigs(n,m-1);
            phin_l     = eigfunc(lambda,m-1,Ltype,Rtype,l(m-2),m,l0,lm,l);
            phin_r     = eigfunc(lambda,m-1,Ltype,Rtype,l(m-1),m,l0,lm,l);
            s2         = s+kappa(m-1)*lambda^2;
            A(m-1,m-2) = A(m-1,m-2)-eigs_norm(n,m-1)^2*phin_l*phin_r/s2;
            A(m-1,m-1) = A(m-1,m-1)+eigs_norm(n,m-1)^2*phin_r*phin_r/s2;
            b(m-1)     = b(m-1)-c(n,m-1)*eigs_norm(n,m-1)*phin_r/s2;
            
            lambda     = eigs(n,m);
            phin_l     = eigfunc(lambda,m,Ltype,Rtype,l(m-1),m,l0,lm,l);
            s1         = s+kappa(m)*lambda^2;
            A(m-1,m-1) = A(m-1,m-1)+eigs_norm(n,m)^2*phin_l*phin_l/s1;
            b(m-1)     = b(m-1)+c(n,m)*eigs_norm(n,m)*phin_l/s1;
            
        end
        
        for i = 1:m-1
            A(i,i) = A(i,i)+(1.0/H(i));
        end
        
        % Laplace transform of v evaluated at zk(poles)/t
        vbar(:,k) = A\b;
%         norm(b - A*vbar(:,k),'inf');
%         pause;
        %vbar(:,k) = zeros(size(b));
        
    end
    
    % Form the sums corresponding to interfaces
    for n = 1:N
        
        % Compute inverse Laplace transform of v(1)/(s+kappa*lambda^2)
        lambda = eigs(n,1);
        vr(1) = 0;
        for k = 1:NZ/2
            poles = 2*k-1;
            s = zk(poles)/t;
            vr(1) = vr(1)-...
                ck(poles)*vbar(1,k)/(t*(s+kappa(1)*lambda^2));
        end
        vr(1) = 2*real(vr(1));
        u(:,1)  = u(:,1)+eigs_norm(n,1)^2*vr(1)*...
            eigfunc(lambda,1,Ltype,Rtype,l(1),m,l0,lm,l)*...
            eigfunc(lambda,1,Ltype,Rtype,xgrid(:,1),m,l0,lm,l);
%         n
%         abs(eigs_norm(n,1)^2*vr(1)*...
%             eigfunc(lambda,1,Ltype,Rtype,l(1),m,l0,lm,l)*...
%             eigfunc(lambda,1,Ltype,Rtype,xgrid(:,1),m,l0,lm,l))
%         pause;
        
        for i = 2:m-1
            lambda = eigs(n,i);
            vr(i-1) = 0;
            for k = 1:NZ/2
                poles = 2*k-1;
                s = zk(poles)/t;
                vr(i-1) = vr(i-1)-...
                    ck(poles)*vbar(i-1,k)/(t*(s+kappa(i)*lambda^2));
            end
            vr(i-1) = 2*real(vr(i-1));
            u(:,i) = u(:,i)-eigs_norm(n,i)^2*vr(i-1)*...
                eigfunc(lambda,i,Ltype,Rtype,l(i-1),m,l0,lm,l)*...
                eigfunc(lambda,i,Ltype,Rtype,xgrid(:,i),m,l0,lm,l);
            vr(i) = 0;
            for k = 1:NZ/2
                poles = 2*k-1;
                s = zk(poles)/t;
                vr(i) = vr(i)-ck(poles)*...
                    vbar(i,k)/(t*(s+kappa(i)*lambda^2));
            end
            vr(i) = 2*real(vr(i));
            u(:,i) = u(:,i)+eigs_norm(n,i)^2*vr(i)*...
                eigfunc(lambda,i,Ltype,Rtype,l(i),m,l0,lm,l)*...
                eigfunc(lambda,i,Ltype,Rtype,xgrid(:,i),m,l0,lm,l);
        end
        
        % Compute inverse Laplace transform of v(m-1)/(s+kappa*lambda^2)
        lambda = eigs(n,m);
        vr(m-1) = 0;
        for k = 1:NZ/2
            poles = 2*k-1;
            s = zk(poles)/t;
            vr(m-1) = vr(m-1)-ck(poles)*...
                vbar(m-1,k)/(t*(s+kappa(m)*lambda^2));
        end
        vr(m-1) = 2*real(vr(m-1));
        u(:,m) = u(:,m)-eigs_norm(n,m)^2*vr(m-1)*...
            eigfunc(lambda,m,Ltype,Rtype,l(m-1),m,l0,lm,l)*...
            eigfunc(lambda,m,Ltype,Rtype,xgrid(:,m),m,l0,lm,l);        
    end
    
%     u = zeros(size(u));
    % Form the sums corresponding to initial condition
    for n = 1:N        
        % Slabs 1,..,m
        for i = 1:m
            lambda = eigs(n,i);
            u(:,i) = u(:,i)+exp(-t*kappa(i)*lambda^2)*c(n,i)*...
                eigs_norm(n,i)*...
                eigfunc(lambda,i,Ltype,Rtype,xgrid(:,i),m,l0,lm,l);
%             if i == 1
%                 n
%                 exp(-t*kappa(i)*lambda^2)*c(n,i)*...
%                 eigs_norm(n,i)*...
%                 eigfunc(lambda,i,Ltype,Rtype,xgrid(:,i),m,l0,lm,l)
%                 pause;
%             end
        end        
    end

    usoln(:,j) = reshape(u,(NX+1)*m,1);
    
end

u = usoln;
x = reshape(xgrid,(NX+1)*m,1);