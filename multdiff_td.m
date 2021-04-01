function varargout = multdiff_td(m,D,gamma,theta,l0,lm,l,u0,Lbnd,Rbnd,tspan,H,varargin)
% MULTDIFF_TD Solves the one-dimensional multilayer diffusion problem using a
%                 Semi-Analytic method.
%
%   MULTDIFF_TD solves the transient diffusion equation in a one-dimensional 
%   composite slab of finite length consisting of multiple layers. The code
%   is applicable to both perfect and imperfect contact at the interfaces 
%   between adjacent layers and time-dependent boundary conditions of Dirichlet,
%   Neumann or Robin at the ends of the slab.
%  
%   MULTDIFF_TD is an implementation of the semi-analytical method
%   proposed by Carr and March.
%
%   Full details can be found in the paper: 
%   EJ Carr and NG March, Semi-analytical solution of multilayer diffusion problems 
%   with time-varying boundary conditions and general interface conditions,
%   Appl. Math. Comput. 333 (2018), 286-303.

% Default values
AbsTol = 1e-10;
RelTol = 1e-6;

% AbsTol = 1e-14;
% RelTol = 1e-10;

% -------------------------------------------------------------------------
% Check inputs
% -------------------------------------------------------------------------
if nargin < 12
    error('Not enough input arguments.'); 
elseif nargin == 12
    options = struct;
elseif nargin == 13
    options = varargin{1};
else
    error('Too many input arguments.');
end

% Number of layers
if round(m) ~= m || m < 2
    error('m must be an integer greater than or equal to 3. For m = 2 with diffusivity D1 on l0 < x < L and D2 on L < x < lm, use m = 4. For example, set D = [D1,D1,D2,D2] and l = [L/2,L,3*L/4].')
end

% Diffusivities
if length(D) ~= m || sum(D > 0) ~= m
    error('D must be a vector of length m with D(i)>0.')
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
if ~isa(u0,'function_handle') || nargin(u0) ~= 2
    error('u0 must be a function handle of the form uint = u0(x).');
end

% Boundary conditions
if ~isa(Lbnd,'cell') || length(Lbnd) ~= 5
    error(['Lbnd must be a cell array of length 5.']);
end
if ~isa(Rbnd,'cell') || length(Rbnd) ~= 5
    error(['Rbnd must be a cell array of length 5.']);
end

% Time vector
tlength = length(tspan);
% if sum(tspan > 0) ~= tlength
%     error('tspan must have entries that are greater than or equal to 0.')
% end

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

% Get boundary condition constants
Ltype = Lbnd{1};
Rtype = Rbnd{1};
aL    = Lbnd{2};
bL    = Lbnd{3};
cL    = Lbnd{4};
Laplace_cL = Lbnd{5};
aR    = Rbnd{2};
bR    = Rbnd{3};
cR    = Rbnd{4};
Laplace_cR = Rbnd{5};

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

% -------------------------------------------------------------------------
% Compute function w(x) that satisfies non-homogeneous BCs


% -------------------------------------------------------------------------
% Eigenvalues
% -------------------------------------------------------------------------
lambda = zeros(N,m);

% Slab 1 (First slab)
switch Ltype
    case 'Dirichlet'
        lambda(:,1) = (2*[0:N-1]+1)'*pi/(2*(l(1)-l0));
    case 'Neumann'
        lambda(:,1) = [0:N-1]'*pi/(l(1)-l0);
    case 'Robin'
        lambda(:,1) = eigvals(l(1),l0,aL/bL,N);
end

% Slab 2,..., m-1 (Middle slabs)
for i = 2:m-1
    lambda(:,i) = [0:N-1]'*pi/(l(i)-l(i-1));
end

% Slab m (End slab)
switch Rtype
    case 'Dirichlet'
        lambda(:,m) = (2*[0:N-1]+1)'*pi/(2*(lm-l(m-1)));
    case 'Neumann'
        lambda(:,m) = [0:N-1]'*pi/(lm-l(m-1));
    case 'Robin'
        lambda(:,m) = eigvals(lm,l(m-1),aR/bR,N);
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
        lambda_val = lambda(n,1);
        eigs_norm(n,1) = sqrt((2*(1+(bL^2/aL^2)*lambda_val^2))/...
            ((bL/aL)+(l(1)-l0)*(1+(bL^2/aL^2)*lambda_val^2)));
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
        lambda_val = lambda(n,m);
        eigs_norm(n,m) = sqrt((2*(1+(bR^2/aR^2)*lambda_val^2))/...
            ((bR/aR)+(lm-l(m-1))*(1+(bR^2/aR^2)*lambda_val^2)));
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
% Expand relevant terms in terms of eigenfunctions
% -------------------------------------------------------------------------
c = zeros(N,m); % Coefficients
beta1 = zeros(N,m);
beta2 = zeros(N,m);
beta3 = zeros(N,m);
beta4 = zeros(N,m);

% Slab 1 (First slab)
for n = 1:N
    prod = @(x) u0(x,1) .* eigfunc(lambda(n,1),1,Ltype,Rtype,x,m,l0,lm,l);
    c(n,1) = eigs_norm(n,1)*integral(prod,l0,l(1),'AbsTol',AbsTol,'RelTol',RelTol);
    
    prod = @(x) psifunc(1,1,Lbnd,Rbnd,x,m,l0,lm,l,gamma) .* ...
        eigfunc(lambda(n,1),1,Ltype,Rtype,x,m,l0,lm,l);
    beta1(n,1) = eigs_norm(n,1)*integral(prod,l0,l(1),'AbsTol',AbsTol,'RelTol',RelTol);
    
    prod = @(x) psifunc(1,2,Lbnd,Rbnd,x,m,l0,lm,l,gamma) .* ...
        eigfunc(lambda(n,1),1,Ltype,Rtype,x,m,l0,lm,l);
    beta2(n,1) = eigs_norm(n,1)*integral(prod,l0,l(1),'AbsTol',AbsTol,'RelTol',RelTol);
    
    prod = @(x) psifunc_2der(1,1,Lbnd,Rbnd,x,m,l0,lm,l,gamma) .* ...
        eigfunc(lambda(n,1),1,Ltype,Rtype,x,m,l0,lm,l);
    beta3(n,1) = eigs_norm(n,1)*integral(prod,l0,l(1),'AbsTol',AbsTol,'RelTol',RelTol);
    
    prod = @(x) psifunc_2der(1,2,Lbnd,Rbnd,x,m,l0,lm,l,gamma) .* ...
        eigfunc(lambda(n,1),1,Ltype,Rtype,x,m,l0,lm,l);
    beta4(n,1) = eigs_norm(n,1)*integral(prod,l0,l(1),'AbsTol',AbsTol,'RelTol',RelTol);
    
end

% Slabs 2,...,m
for i = 2:m-1
    for n = 1:N
        prod = @(x) u0(x,i) .* eigfunc(lambda(n,i),i,Ltype,Rtype,x,m,l0,lm,l);
        c(n,i) = eigs_norm(n,i)*integral(prod,l(i-1),l(i),'AbsTol',AbsTol,'RelTol',RelTol);
        
        prod = @(x) psifunc(i,1,Lbnd,Rbnd,x,m,l0,lm,l,gamma) .* ...
            eigfunc(lambda(n,i),i,Ltype,Rtype,x,m,l0,lm,l);
        beta1(n,i) = eigs_norm(n,i)*integral(prod,l(i-1),l(i),'AbsTol',AbsTol,'RelTol',RelTol);
        
        prod = @(x) psifunc(i,2,Lbnd,Rbnd,x,m,l0,lm,l,gamma) .* ...
            eigfunc(lambda(n,i),i,Ltype,Rtype,x,m,l0,lm,l);
        beta2(n,i) = eigs_norm(n,i)*integral(prod,l(i-1),l(i),'AbsTol',AbsTol,'RelTol',RelTol);
        
        prod = @(x) psifunc_2der(i,1,Lbnd,Rbnd,x,m,l0,lm,l,gamma) .* ...
            eigfunc(lambda(n,i),i,Ltype,Rtype,x,m,l0,lm,l);
        beta3(n,i) = eigs_norm(n,i)*integral(prod,l(i-1),l(i),'AbsTol',AbsTol,'RelTol',RelTol);
        
        prod = @(x) psifunc_2der(i,2,Lbnd,Rbnd,x,m,l0,lm,l,gamma) .* ...
            eigfunc(lambda(n,i),i,Ltype,Rtype,x,m,l0,lm,l);
        beta4(n,i) = eigs_norm(n,i)*integral(prod,l(i-1),l(i),'AbsTol',AbsTol,'RelTol',RelTol);        
        
    end
end

% Slab m
for n = 1:N
    prod = @(x) u0(x,m) .* eigfunc(lambda(n,m),m,Ltype,Rtype,x,m,l0,lm,l);
    c(n,m) = eigs_norm(n,m)*integral(prod,l(m-1),lm,'AbsTol',AbsTol,'RelTol',RelTol);
    
    prod = @(x) psifunc(m,1,Lbnd,Rbnd,x,m,l0,lm,l,gamma) .* ... 
        eigfunc(lambda(n,m),m,Ltype,Rtype,x,m,l0,lm,l);
    beta1(n,m) = eigs_norm(n,m)*integral(prod,l(m-1),lm,'AbsTol',AbsTol,'RelTol',RelTol);
    
    prod = @(x) psifunc(m,2,Lbnd,Rbnd,x,m,l0,lm,l,gamma) .* ... 
        eigfunc(lambda(n,m),m,Ltype,Rtype,x,m,l0,lm,l);
    beta2(n,m) = eigs_norm(n,m)*integral(prod,l(m-1),lm,'AbsTol',AbsTol,'RelTol',RelTol);
    
    prod = @(x) psifunc_2der(m,1,Lbnd,Rbnd,x,m,l0,lm,l,gamma) .* ... 
        eigfunc(lambda(n,m),m,Ltype,Rtype,x,m,l0,lm,l);
    beta3(n,m) = eigs_norm(n,m)*integral(prod,l(m-1),lm,'AbsTol',AbsTol,'RelTol',RelTol);
    
    prod = @(x) psifunc_2der(m,2,Lbnd,Rbnd,x,m,l0,lm,l,gamma) .* ... 
        eigfunc(lambda(n,m),m,Ltype,Rtype,x,m,l0,lm,l);
    beta4(n,m) = eigs_norm(n,m)*integral(prod,l(m-1),lm,'AbsTol',AbsTol,'RelTol',RelTol);    
    
end

% Plot initial condition
u = zeros(NX+1,m);
for n = 1:N
    for i = 1:m
        u(:,i) = u(:,i)+c(n,i)*...
            eigs_norm(n,i)*...
            eigfunc(lambda(n,i),i,Ltype,Rtype,xgrid(:,i),m,l0,lm,l);
    end
end

% Get weights and poles for use in inverse transform
[zk,ck] = cf(NZ);

% usoln = zeros((NX+1)*m,tlength);
% Solution (at given time)
u = zeros(NX+1,m,tlength);

%figure; 
% Time loop
for j = 1:tlength
    
    t = tspan(j);
    
    % Form the sums corresponding to initial condition
    for n = 1:N        
        % Slabs 1,..,m
        for i = 1:m
            u(:,i,j) = u(:,i,j)+exp(-t*D(i)*lambda(n,i)^2)*c(n,i)*...
                eigs_norm(n,i)*...
                eigfunc(lambda(n,i),i,Ltype,Rtype,xgrid(:,i),m,l0,lm,l);
        end        
    end
    
    % Compute inverse Laplace transform of interface functions
    gbar = zeros(m-1,NZ/2);
    
    for k = 1:NZ/2
        
        A    = zeros(m-1,m-1);
        vr   = zeros(m-1,1);
        b    = zeros(m-1,1);
        
        poles = 2*k-1;
        s = zk(poles)/t;
        
        A(1,1) = psifunc(1,2,Lbnd,Rbnd,l(1),m,l0,lm,l,gamma) ...
            - theta(1)*psifunc(2,1,Lbnd,Rbnd,l(1),m,l0,lm,l,gamma);
        A(1,2) = -theta(1)*psifunc(2,2,Lbnd,Rbnd,l(1),m,l0,lm,l,gamma);
        %Laplace_cL = cL(t)/s;
        %Laplace_cL = integral(@(t) cL(t).*exp(-s*t),0,Inf)
        b(1) = -Laplace_cL(s)*psifunc(1,1,Lbnd,Rbnd,l(1),m,l0,lm,l,gamma);
        
        for i = 2:m-2
            A(i,i-1) = psifunc(i,1,Lbnd,Rbnd,l(i),m,l0,lm,l,gamma);
            A(i,i)   = psifunc(i,2,Lbnd,Rbnd,l(i),m,l0,lm,l,gamma) ...
                - theta(i)*psifunc(i+1,1,Lbnd,Rbnd,l(i),m,l0,lm,l,gamma);
            A(i,i+1) = -theta(i)*psifunc(i+1,2,Lbnd,Rbnd,l(i),m,l0,lm,l,gamma);
        end
        
        if m > 2
            A(m-1,m-2) = psifunc(m-1,1,Lbnd,Rbnd,l(m-1),m,l0,lm,l,gamma);
            A(m-1,m-1) = psifunc(m-1,2,Lbnd,Rbnd,l(m-1),m,l0,lm,l,gamma) ...
                - theta(m-1)*psifunc(m,1,Lbnd,Rbnd,l(m-1),m,l0,lm,l,gamma);
            %Laplace_cR = cR(t)/s;
            %Laplace_cR = integral(@(t) cR(t).*exp(-s*t),0,Inf);
            b(m-1) = Laplace_cR(s)*theta(m-1)*psifunc(m,2,Lbnd,Rbnd,l(m-1),m,l0,lm,l,gamma);
        end
        
        for n = 1:N
            
            % Interface 1 (between layers 1 and 2)
            phi1   = eigs_norm(n,1) ...
                *eigfunc(lambda(n,1),1,Ltype,Rtype,l(1),m,l0,lm,l);
            phi2   = eigs_norm(n,2) ...
                *eigfunc(lambda(n,2),2,Ltype,Rtype,l(1),m,l0,lm,l);
            s1     = s+D(1)*lambda(n,1)^2;
            s2     = s+D(2)*lambda(n,2)^2;
            A(1,1) = A(1,1) ...
                + (D(1)*(beta4(n,1)+lambda(n,1)^2*beta2(n,1))/s1 - beta2(n,1))*phi1 ...
                - theta(1)*(D(2)*(beta3(n,2)+lambda(n,2)^2*beta1(n,2))/s2 - beta1(n,2))*phi2;
            A(1,2) = A(1,2) ...
                - theta(1)*(D(2)*(beta4(n,2)+lambda(n,2)^2*beta2(n,2))/s2-beta2(n,2))*phi2;
            b(1)   = b(1) + theta(1)*c(n,2)*phi2/s2 - c(n,1)*phi1/s1 ...
                - Laplace_cL(s) * (D(1)*(beta3(n,1)+lambda(n,1)^2*beta1(n,1))/s1-beta1(n,1))*phi1;
            
            % Middle interfaces
            for i = 2:m-2
                phi1   = eigs_norm(n,i) ...
                    *eigfunc(lambda(n,i),i,Ltype,Rtype,l(i),m,l0,lm,l);
                phi2   = eigs_norm(n,i+1) ...
                    *eigfunc(lambda(n,i+1),i+1,Ltype,Rtype,l(i),m,l0,lm,l);
                s1     = s+D(i)*lambda(n,i)^2;
                s2     = s+D(i+1)*lambda(n,i+1)^2;
                A(i,i-1) = A(i,i-1) ...
                    - beta1(n,i)*phi1 ...
                    + D(i)*phi1*(beta3(n,i)+lambda(n,i)^2*beta1(n,i))/s1;
                A(i,i) = A(i,i) ...
                    + (D(i)*(beta4(n,i)+lambda(n,i)^2*beta2(n,i))/s1- beta2(n,i))*phi1 ...
                    - theta(i)*(D(i+1)*(beta3(n,i+1)+lambda(n,i+1)^2*beta1(n,i+1))/s2-beta1(n,i+1))*phi2;
                A(i,i+1) = A(i,i+1) ...
                    - theta(i)*(D(i+1)*(beta4(n,i+1)+lambda(n,i+1)^2*beta2(n,i+1))/s2-beta2(n,i+1))*phi2;
                b(i) = b(i) + theta(i)*c(n,i+1)*phi2/s2 - c(n,i)*phi1/s1;
            end
            
            if m > 2
                % Interface m (between layers m and m-1)
                phi1   = eigs_norm(n,m-1) ...
                    *eigfunc(lambda(n,m-1),m-1,Ltype,Rtype,l(m-1),m,l0,lm,l);
                phi2   = eigs_norm(n,m) ...
                    *eigfunc(lambda(n,m),m,Ltype,Rtype,l(m-1),m,l0,lm,l);
                s1     = s+D(m-1)*lambda(n,m-1)^2;
                s2     = s+D(m)*lambda(n,m)^2;
                A(m-1,m-2) = A(m-1,m-2) ...
                    + (D(m-1)*(beta3(n,m-1)+lambda(n,m-1)^2*beta1(n,m-1))/s1-beta1(n,m-1))*phi1;
                A(m-1,m-1) = A(m-1,m-1) ...
                    + (D(m-1)*(beta4(n,m-1)+lambda(n,m-1)^2*beta2(n,m-1))/s1 - beta2(n,m-1))*phi1...
                    - theta(m-1)*(D(m)*(beta3(n,m)+lambda(n,m)^2*beta1(n,m))/s2 -beta1(n,m))*phi2;
                b(m-1) = b(m-1) + theta(m-1)*c(n,m)*phi2/s2 - c(n,m-1)*phi1/s1 ...
                    + Laplace_cR(s) * theta(m-1)*(D(m)*(beta4(n,m)+lambda(n,m)^2*beta2(n,m))/s2-beta2(n,m))*phi2;
            end
        end
        
        for i = 1:m-1
            A(i,i) = A(i,i)+1.0/H(i);
        end
        
        
%         A,b,s
%         pause;
        
        % Laplace transform of v evaluated at zk(poles)/t
        %size(A),size(b)
        gbar(:,k) = A\b;
        
    end
    
    % Form the sums corresponding to interfaces
    % Find interface functions: g_{i}(t)
    g = zeros(m-1,1);
    for k = 1:NZ/2
        poles = 2*k-1;
        g = g-ck(poles)*gbar(:,k)/t;
    end
    g = 2*real(g);
    g = [cL(t); g; cR(t)];
    gint(:,j) = g;
    
    for i = 1:m
        u(:,i,j) = u(:,i,j)...
            + g(i)*psifunc(i,1,Lbnd,Rbnd,xgrid(:,i),m,l0,lm,l,gamma) ...
            + g(i+1)*psifunc(i,2,Lbnd,Rbnd,xgrid(:,i),m,l0,lm,l,gamma);
    end
    
    for n = 1:N
        
        % Compute inverse Laplace transform of v(1)/(s+D*lambda^2)
%         int1 = 1/(D(1)*lambda(n,1)^2)
%         int1 = integral(@(tau) cL(tau).*exp(-(t-tau)*D(1)*lambda(n,1)^2),0,t)
%         pause;

        int1 = 0;
        for k = 1:NZ/2
            poles = 2*k-1;
            s = zk(poles)/t;
            int1 = int1-...
                ck(poles)*Laplace_cL(s)/(t*(s+D(1)*lambda(n,1)^2));
        end
        int1 = 2*real(int1);
%         pause;

        u(:,1,j)  = u(:,1,j)+D(1)*(beta3(n,1)+lambda(n,1)^2*beta1(n,1)) * ...
            int1 * ...
            eigs_norm(n,1)*eigfunc(lambda(n,1),1,Ltype,Rtype,xgrid(:,1),m,l0,lm,l);
        vr(1) = 0;
        for k = 1:NZ/2
            poles = 2*k-1;
            s = zk(poles)/t;
            vr(1) = vr(1)-...
                ck(poles)*gbar(1,k)/(t*(s+D(1)*lambda(n,1)^2));
        end
        vr(1) = 2*real(vr(1));
        u(:,1,j)  = u(:,1,j)+eigs_norm(n,1)*vr(1)*...
            D(1)*(beta4(n,1)+lambda(n,1)^2*beta2(n,1))*...
            eigfunc(lambda(n,1),1,Ltype,Rtype,xgrid(:,1),m,l0,lm,l);
        
        u(:,1,j) = u(:,1,j)...
            -g(1)*beta1(n,1)*eigs_norm(n,1)*eigfunc(lambda(n,1),1,Ltype,Rtype,xgrid(:,1),m,l0,lm,l)...
            -g(2)*beta2(n,1)*eigs_norm(n,1)*eigfunc(lambda(n,1),1,Ltype,Rtype,xgrid(:,1),m,l0,lm,l);
        
        for i = 2:m-1
            vr(i-1) = 0;
            for k = 1:NZ/2
                poles = 2*k-1;
                s = zk(poles)/t;
                vr(i-1) = vr(i-1)-...
                    ck(poles)*gbar(i-1,k)/(t*(s+D(i)*lambda(n,i)^2));
            end
            vr(i-1) = 2*real(vr(i-1));
            u(:,i,j) = u(:,i,j)+eigs_norm(n,i)*vr(i-1)*...
                D(i)*(beta3(n,i)+lambda(n,i)^2*beta1(n,i))*...
                eigfunc(lambda(n,i),i,Ltype,Rtype,xgrid(:,i),m,l0,lm,l);
            vr(i) = 0;
            for k = 1:NZ/2
                poles = 2*k-1;
                s = zk(poles)/t;
                vr(i) = vr(i)-ck(poles)*...
                    gbar(i,k)/(t*(s+D(i)*lambda(n,i)^2));
            end
            vr(i) = 2*real(vr(i));
            u(:,i,j) = u(:,i,j)+eigs_norm(n,i)*vr(i)*...
                D(i)*(beta4(n,i)+lambda(n,i)^2*beta2(n,i))*...
                eigfunc(lambda(n,i),i,Ltype,Rtype,xgrid(:,i),m,l0,lm,l);
            
            u(:,i,j) = u(:,i,j)...
                -g(i)*beta1(n,i)*eigs_norm(n,i)*eigfunc(lambda(n,i),i,Ltype,Rtype,xgrid(:,i),m,l0,lm,l)...
                -g(i+1)*beta2(n,i)*eigs_norm(n,i)*eigfunc(lambda(n,i),i,Ltype,Rtype,xgrid(:,i),m,l0,lm,l);
            
        end
        
        % Compute inverse Laplace transform of v(m-1)/(s+D*lambda^2)
        u(:,m,j)  = u(:,m,j)+D(m)*(beta4(n,m)+lambda(n,m)^2*beta2(n,m)) * ...
            integral(@(tau) cR(tau).*exp(-(t-tau)*D(m)*lambda(n,m)^2),0,t) * ...
            eigs_norm(n,m)*eigfunc(lambda(n,m),m,Ltype,Rtype,xgrid(:,m),m,l0,lm,l);
        vr(m-1) = 0;
        for k = 1:NZ/2
            poles = 2*k-1;
            s = zk(poles)/t;
            vr(m-1) = vr(m-1)-ck(poles)*...
                gbar(m-1,k)/(t*(s+D(m)*lambda(n,m)^2));
        end
        vr(m-1) = 2*real(vr(m-1));
        u(:,m,j) = u(:,m,j)+D(m)*(beta3(n,m)+lambda(n,m)^2*beta1(n,m)) * vr(m-1) * ...
            eigs_norm(n,m)*eigfunc(lambda(n,m),m,Ltype,Rtype,xgrid(:,m),m,l0,lm,l);
        
        u(:,m,j) = u(:,m,j)...
            -g(m)*beta1(n,m)*eigs_norm(n,m)*eigfunc(lambda(n,m),m,Ltype,Rtype,xgrid(:,m),m,l0,lm,l)...
            -g(m+1)*beta2(n,m)*eigs_norm(n,m)*eigfunc(lambda(n,m),m,Ltype,Rtype,xgrid(:,m),m,l0,lm,l);
                    
    end
    
end

% u = usoln;
x = xgrid;

varargout{1} = u;
varargout{2} = x;
if nargout == 3
    varargout{3} = gint(2:m,:); % Return interface fluxes
end
