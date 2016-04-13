function [u,x] = multdiff_fvm(m,kappa,l0,lm,l,u0,Lbnd,Rbnd,tspan,interface,dt,varargin)
% MULTDIFF_FVM Solves the one-dimensional multilayer diffusion problem
%              using the Finite Volume Method.
%
%   MULTDIFF_FVM solves the transient diffusion equation in a 
%   one-dimensional composite slab of finite length consisting of multiple 
%   layers. The code is applicable to both perfect and imperfect contact at
%   the interfaces between adjacent layers and either Dirichlet, Neumann or
%   Robin boundary conditions at the ends of the slab.
%  
%   MULTDIFF_FVM uses the vertex-centered Finite Volume Method together 
%   with a backward Euler discretisation in time using a fixed time 
%   stepsize. 
%
%   Full details can be found in the paper: 
%   E. J. Carr and I. W. Turner.
%
%   Description:
%   -----------------------------------------------------------------------
%   MULTDIFF_FVM solves the standard diffusion equation in each layer 
%   (l(i-1) < x < l(i)):
%
%      du_(i)/dt = d/dx * (kappa(i) * du_(i)/dx),   i = 1,...,m,
%   
%   subject to the following initial and external boundary conditions:
%   
%      u_(i)(x,t) = u0(x)                           at t = 0
%      aL * u_(1)(x,t) - bL * du_(1)/dx(x,t) = cL   at x = l0
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
%       kappa(i) * du_(i)/dx(x,t) = kappa(i+1) * du_(i+1)/dx(x,t)
%
%    - Imperfect contact
%       kappa(i) * du_(i)/dx(x,t) = H(i) * (u_(i+1)(x,t) - u_(i)(x,t))
%       kappa(i+1) * du_(i+1)/dx(x,t) = H(i) * (u_(i+1)(x,t) - u_(i)(x,t))       
%   
%   Usage:
%   -----------------------------------------------------------------------
%   [U,X] = multdiff_fvm(m,kappa,l0,lM,l,u0,Lbnd,Rbnd,tspan,'Perfect',dt)
%   [U,X] = multdiff_fvm(m,kappa,l0,lM,l,u0,Lbnd,Rbnd,tspan,'Perfect',dt,options)
%   [U,X] = multdiff_fvm(m,kappa,l0,lM,l,u0,Lbnd,Rbnd,tspan,'Membrane',dt,gamma)
%   [U,X] = multdiff_fvm(m,kappa,l0,lM,l,u0,Lbnd,Rbnd,tspan,'Membrane',dt,gamma,options)
%   [U,X] = multdiff_fvm(m,kappa,l0,lM,l,u0,Lbnd,Rbnd,tspan,'Imperfect',dt,H)
%   [U,X] = multdiff_fvm(m,kappa,l0,lM,l,u0,Lbnd,Rbnd,tspan,'Imperfect',dt,H,options)
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
%               t0,t1,...,tf, use TSPAN = [t0,t1,...,tf]. t0,t1,...,tf must
%               be multiples of dt.
%   interface   Internal boundary conditions at interfaces between adjacent
%               layers. inteface can be either 'Perfect' or 'Imperfect'.
%   dt          Time stepsize to use in the backward Euler discretisation.
%   H           A vector of length m-1 containing the contact 
%               transfer coeffecients at the interfaces between adjacent 
%               layers such that the coefficient between layer i and 
%               layer i+1 is given by H(i) (i = 1,..,m-1). 
%               * Applicable to imperfect contant only.
%   options     An (optional) set of solver options. Fields in the 
%               structure options are
%                - NX   number of divisions within each slab. U(:,j) gives
%                       the solution at x = l(i-1):(l(i)-l(i-1))/NX:l(i) 
%                       and t = tspan(j).
%                       [NX = 50 by default]
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
%   bcs.Ltype = 'Dirichlet'; bcs.aL = 1.0; bcs.bL = 0.0; bcs.cL = 1.0;
%   bcs.Rtype = 'Dirichlet'; bcs.aR = 1.0; bcs.bR = 0.0; bcs.cR = 0.5;
%   [u,x] = multdiff_fvm(3,[1,0.1,1],0.0,1.0,[0.3,0.7],u0,bcs,...
%           [0.02,0.05,0.1,0.2,0.5,1.0],'Perfect',0.01);
%
% -------------------------------------------------------------------------
% Check inputs
% -------------------------------------------------------------------------
if nargin < 11
    error('Not enough input arguments.');
elseif nargin == 11
    if strcmp(interface,'Imperfect')
        error('H must be specified for imperfect contact at interfaces.');
    end
    options = struct;    
elseif nargin == 12
    if strcmp(interface,'Perfect')
        options = varargin{1};
    elseif strcmp(interface,'Imperfect')
        H = varargin{1};
        options = struct;
    elseif strcmp(interface,'Membrane')
        gamma = varargin{1};
        options = struct;
    end
elseif nargin == 13
    if strcmp(interface,'Perfect')
        error('Too many input arguments for interface = ''Perfect''.');
    elseif strcmp(interface,'Imperfect')
        H = varargin{1};
        options = varargin{2};
    elseif strcmp(interface,'Membrane')
        gamma = varargin{1};
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
    %warning('u0 must be a function handle of the form uint = u0(x).');
end

% Boundary conditions
if ~isa(Lbnd,'cell') || length(Lbnd) ~= 4
    error('Lbnd must be a cell array of length 4.');
end
if ~isa(Rbnd,'cell') || length(Rbnd) ~= 4
    error('Rbnd must be a cell array of length 4.');
end

% Time vector
tlength = length(tspan);
if sum(tspan > 0) ~= tlength
    error('tspan must have entries that are greater than or equal to 0.')
end

% Internal boundary conditions at interfaces
if strcmp(interface,'Perfect') || strcmp(interface,'Imperfect') || ...
        strcmp(interface,'Membrane')
else
    error('interface must be either ''Perfect'', ''Imperfect'' and ''Membrane''.')
end

% Check time step
if dt < 0
    error('dt must be positive.')
end
%if sum(mod(abs(tspan/dt),1)) ~= 0
%    error('tspan must have entries that are multiples of dt.');
%end

% Check options structure
if ~isa(options,'struct')
    error('options must be a structure.')
end
Names = {'NX','Hp'};
fn = fieldnames(options);
for i = 1:length(fn)
    j = strcmp(fn(i),Names);
    if sum(j) == 0
        error('Invalid option ''%s''.',fn{i});
    end
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
% Value of contact transfer coefficient to approximate perfect contact 
% condition
if isfield(options,'Hp')
    if strcmp(interface,'Perfect') || strcmp(interface,'Membrane')
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
    H = Hp*ones(m-1,1);
    gamma = ones(m-1,1);
elseif strcmp(interface,'Imperfect')
    gamma = ones(m-1,1);
elseif strcmp(interface,'Membrane')
    H = Hp*ones(m-1,1);
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
if strcmp(Ltype,'Neumann') && aL ~= 0
    error('Neumann condition at left boundary cannot have aL = 0.');
end
if strcmp(Rtype,'Neumann') && aR ~= 0
    error('Neumann condition at right boundary cannot have aR = 0.');
end
% if (aL == 0 || bL == 0) && strcmp(Ltype,'Robin')
%     error('Robin condition at left boundary cannot have aL = 0 or bL = 0');
% end
% if (aR == 0 || bR == 0) && strcmp(Rtype,'Robin')
%     error('Robin condition at left boundary cannot have aR = 0 or bR = 0');
% end
if strcmp(Ltype,'Robin') && aL/bL < 0
%     aL,bL
%     error('Robin condition at left boundary must have aL/bL > 0.');
end
if strcmp(Rtype,'Robin') && aR/bR < 0
%     error('Robin condition at right boundary must have aR/bR > 0.');
end
if aL == 0 && bL == 0
    error('Boundary condition is incorrect at left boundary (aL = bL = 0).');
end
if aR == 0 && bR == 0
    error('Boundary condition is incorrect at left boundary (aR = bR = 0).');
end

% -------------------------------------------------------------------------
% FVM geometric properties
% -------------------------------------------------------------------------
xgrid = zeros(NX+1,m);

% Layer 1
xgrid(:,1) = l0:(l(1)-l0)/NX:l(1);

% Layer 2,...,m-1
for i = 2:m-1
    xgrid(:,i) = l(i-1):(l(i)-l(i-1))/NX:l(i);
end
% Layer m
xgrid(:,m) = l(m-1):(lm-l(m-1))/NX:lm;

% Node spacing
dx = diff(xgrid);

% Control volume widths
cv = zeros(size(xgrid));

for i = 1:m
    cv(1,i) = dx(1,i)/2;
    for j = 2:NX
        cv(j,i) = dx(j-1,i)/2 + dx(j,i)/2;
    end
    cv(NX+1,i) = dx(NX,i)/2;
end

% -------------------------------------------------------------------------
% Form FVM matrix
% -------------------------------------------------------------------------
N = (NX+1)*m; % No of unknowns
A = zeros(N,N);
b = zeros(N,1);

% Left boundary
A(1,1) = bL - (dt*kappa(1)/cv(1,1))*(-aL - bL/dx(1,1));
A(1,2) = -dt*kappa(1)*bL/(cv(1,1)*dx(1,1));

% Right boundary
A(N,N)   = bR + (dt*kappa(m)/cv(NX+1,m))*(aR + bR/dx(NX,m));
A(N,N-1) = -dt*kappa(m)*bR/(cv(NX+1,m)*dx(NX,m));

% Internal nodes within each layer
for j = 1:m
    for i = 2:NX
        r        = (j-1)*(NX+1)+i;
        A(r,r-1) = -dt*kappa(j)/(cv(i,j)*dx(i-1,j));
        A(r,r)   = 1.0 + dt*kappa(j)/cv(i,j)*(1/dx(i-1,j) + 1/dx(i,j));
        A(r,r+1) = -dt*kappa(j)/(cv(i,j)*dx(i,j));
    end
end

% Interface node (left end of layer)
for j = 2:m
    i = 1;
    r = (j-1)*(NX+1)+i;
    if strcmp(interface,'Perfect')
        cv_area  = cv(i,j)+cv(NX+1,j-1);
        A(r,r-2) = -(dt/cv_area)*kappa(j-1)/dx(NX,j-1);
        A(r,r)   = 1.0 + (dt/cv_area)*(kappa(j)/dx(i,j) + kappa(j-1)/dx(NX,j-1));
        A(r,r+1) = -(dt/cv_area)*kappa(j)/dx(i,j);
    elseif strcmp(interface,'Imperfect')
        A(r,r-1) = -dt/cv(i,j);
        A(r,r)   = (1.0/H(j-1)) + dt/cv(i,j)*(kappa(j)/(H(j-1)*dx(i,j)) + 1.0);
        A(r,r+1) = -(dt/cv(i,j))*kappa(j)/(H(j-1)*dx(i+1,j));
    end
end

% Interface node (right end of layer)
for j = 1:m-1
    i = NX+1;
    r = (j-1)*(NX+1)+i;
    if strcmp(interface,'Perfect')
%         cv_area  = cv(i,j)+cv(1,j+1);
%         A(r,r-1) = -(dt/cv_area)*kappa(j)/dx(i-1,j);
%         A(r,r)   = 1.0 + (dt/cv_area)*(kappa(j)/dx(i-1,j) + kappa(j+1)/dx(1,j+1));
%         A(r,r+1) = -(dt/cv_area)*kappa(j+1)/dx(1,j+1); 
        A(r,r)    = 1.0;
        A(r,r+1)  = -1.0;
    elseif  strcmp(interface,'Imperfect')
        A(r,r-1) = -(dt/cv(i,j))*kappa(j)/(H(j)*dx(i-1,j));
        A(r,r)   = (1.0/H(j)) + dt/cv(i,j)*(kappa(j)/(H(j)*dx(i-1,j)) + 1.0);
        A(r,r+1) = -dt/cv(i,j);
    end
end

% A
% pause;

% if strcmp(interface,'Perfect')
%     A([1:m-1]*(NX+1),:) = [];
%     %A(:,[1:m-1]*(NX+1)) = [];
% end

% A

% pause;

A = sparse(A);

% -------------------------------------------------------------------------
% Time loop
% -------------------------------------------------------------------------
t = 0;
k = 1;
tend = tspan(end);

% Initial condition
if isa(u0,'function_handle')
    u = u0(xgrid);
elseif sum(size(u0) == [NX+1,m]) == 2
    u = u0;
else
    error('u0 not implemented correctly');
end

usoln = zeros((NX+1)*m,length(tspan));
solntemp = zeros((NX+1)*m,1);
xgrid = reshape(xgrid,(NX+1)*m,1);

while t <= tend
    
    t = t + dt;
    
    % Left boundary
    b(1) = bL*u(1,1)+(dt*kappa(1)*cL/cv(1,1));
    
    % Right boundary
    b(N) = bR*u(NX+1,m)+(dt*kappa(m)*cR/cv(NX+1,m));
    
    % Middles nodes within each layer
    for j = 1:m
       for i = 2:NX
           r = (j-1)*(NX+1)+i;
           b(r) = u(i,j);
       end
    end

    % Middles nodes within each layer
%     for j = 1:m
% %        for i = 2:NX
%            i = 2:NX;
%            r = (j-1)*(NX+1)+i;
%            b(r) = u(i,j);
% %        end
%     end


    % Interface node (left)
    for j = 2:m
        i = 1;
        r = (j-1)*(NX+1)+i;
        if strcmp(interface,'Perfect')
            b(r) = u(i,j);
        elseif strcmp(interface,'Imperfect')
            b(r) = u(i,j)/H(j-1);
        end
            
    end
    
    % Interface node (right)
    for j = 1:m-1
        i = NX+1;
        r = (j-1)*(NX+1)+i;
        if strcmp(interface,'Perfect')
            b(r) = 0.0;
        elseif strcmp(interface,'Imperfect')
            b(r) = u(i,j)/H(j);
        end
    end

%     if strcmp(interface,'Perfect')
%         b([1:m-1]*(NX+1)) = [];
%     end
%     b
%     b = max(b,0);
%     pause;
%     b
% kappa
% [full(A),b]
    soln = A\b;
%     pause;
    
%     if strcmp(interface,'Perfect')
%         for j = 1:m-1
%             for i = 1:NX+1
%                 if i == NX+1
%                     solntemp((j-1)*(NX+1)+i) = soln((j-1)*NX+NX);
%                 else
%                     solntemp((j-1)*(NX+1)+i) = soln((j-1)*NX+i);
%                 end
%             end
%         end
%         j = m;
%         for i = 1:NX+1
%             solntemp((j-1)*(NX+1)+i) = soln((j-1)*NX+i);
%         end
%         soln = solntemp;
%     end
    u = reshape(soln,NX+1,m);
    %pause;
    
    % Store solution at requested values in tspan
    if abs(tspan(k)-t) < 1e-6*dt
        usoln(:,k) = soln;
        k = k + 1;
    end
    
    if k == length(tspan)+1
        break;
    end
    
end

%eigs = eig(full(A));
%figure;
%plot(real(eigs),imag(eigs),'.')

u = usoln;
x = reshape(xgrid,(NX+1)*m,1);

