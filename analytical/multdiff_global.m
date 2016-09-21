function [u,x] = multdiff_global(m,kappa,l0,lm,l,u0,Lbnd,Rbnd,tspan,interface,varargin)

%% Check inputs
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

%% Get boundary condition constants
Ltype = Lbnd{1};
Rtype = Rbnd{1};
aL    = Lbnd{2};
bL    = Lbnd{3};
cL    = Lbnd{4};
aR    = Rbnd{2};
bR    = Rbnd{3};
cR    = Rbnd{4};

N     = options.N; % Number of eigenvalues
NX    = options.NX; % Number of divisions in each slab (to determine x values 
                    % at which the solution is computed)

%% Contact transfer coefficient to approximate perfect contact condition
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
    Hp = Inf; % Default
end
if strcmp(interface,'Perfect')
    H = Hp*ones(m-1,m);
end

%% Compute function w(x) that satisfies non-homogeneous BCs
A = zeros(2*m,2*m);
b = zeros(2*m,1);

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

if strcmp(Ltype,'Neumann') && strcmp(Rtype,'Neumann')
    i = 1;
    A(2*m+1,2*i-1) = l(1) - l0;
    A(2*m+1,2*i)   = (l(1)^2-l0^2)/2;
    b(2*m+1)       = integral(u0,l0,l(1));
    for i = 2:m-1
        A(2*m+1,2*i-1) = l(i)-l(i-1);
        A(2*m+1,2*i)   = (l(i)^2-l(i-1)^2)/2;
        b(2*m+1)       = b(2*m+1) + integral(u0,l(i-1),l(i));
    end
    i = m;
    A(2*m+1,2*i-1) = lm-l(m-1);
    A(2*m+1,2*i)   = (lm^2-l(m-1)^2)/2;
    b(2*m+1)       = b(2*m+1) + integral(u0,l(m-1),lm);
end

w = A\b;

%% Eigenvalues
detA = @(lambda) det(Amat(lambda,m,kappa,l0,lm,l,aL,bL,aR,bR,H));

intx = [0.2,100.0];
g = chebfun(detA,intx,'vectorize');
eigvals = roots(g,'nozerofun','norecursion');

neigvals = length(eigvals);
i = 0;
while neigvals < N
    i = i + 1;
    inter = intx+[0.0,i*10.0];
    g = chebfun(detA,inter,'vectorize'); % Case B
    eigvals = roots(g);
    neigvals = length(eigvals);
end

eigvals = eigvals(1:N);

coeffs = zeros(2*m,length(eigvals));

for i = 1:N
    A = Amat(eigvals(i),m,kappa,l0,lm,l,aL,bL,aR,bR,H);
    [X,D] = eig(A);
    [~,indx] = min(abs(diag(D)));
    coeffs(:,i) = real(X(:,indx)/norm(X(:,indx),2));
end

%% Eigenfunction normalisation constants
Nscal = zeros(N,1);

for n = 1:N
    Nscal(n) = 0.0;
    lambda = eigvals(n);
    for i = 1:m
        if i == 1
            func = @(x) (coeffs(2*i-1,n)*sin(lambda*(x-l0)/sqrt(kappa(1))) + ...
                coeffs(2*i,n)*cos(lambda*(x-l0)/sqrt(kappa(1)))).^2;
            Nscal(n) = Nscal(n) + integral(func,l0,l(1));
        elseif i == m
            func = @(x) (coeffs(2*i-1,n)*sin(lambda*(x-l(m-1))/sqrt(kappa(m))) + ...
                coeffs(2*i,n)*cos(lambda*(x-l(m-1))/sqrt(kappa(m)))).^2;
            Nscal(n) = Nscal(n) + integral(func,l(m-1),lm);
        else
            func = @(x) (coeffs(2*i-1,n)*sin(lambda*(x-l(i-1))/sqrt(kappa(i))) + ...
                coeffs(2*i,n)*cos(lambda*(x-l(i-1))/sqrt(kappa(i)))).^2;
            Nscal(n) = Nscal(n) + integral(func,l(i-1),l(i));
        end
    end
end

%% Grid spacing within each slab
xgrid = zeros(NX+1,m);
xgrid(:,1) = l0:(l(1)-l0)/NX:l(1); % Slab 1
for i = 2:m-1
    xgrid(:,i) = l(i-1):(l(i)-l(i-1))/NX:l(i); % Slabs 2,...,m
end
xgrid(:,m) = l(m-1):(lm-l(m-1))/NX:lm; % Slab m

%% Initial conditions - expand in terms of eigenfunctions
c = zeros(N,1);

for n = 1:N
    lambda = eigvals(n);
    for i = 1:m
        if i == 1
            prodt = @(x) ((u0(x)-(w(2*i-1) + w(2*i)*x)) .* ...
                (coeffs(2*i-1,n)*sin(lambda*(x-l0)/sqrt(kappa(1))) + ...
                coeffs(2*i,n)*cos(lambda*(x-l0)/sqrt(kappa(1)))));
            c(n) = c(n) + integral(prodt,l0,l(1));
        elseif i == m
            prodt = @(x) ((u0(x)-(w(2*i-1) + w(2*i)*x)) .* ...
                (coeffs(2*i-1,n)*sin(lambda*(x-l(m-1))/sqrt(kappa(m))) + ...
                coeffs(2*i,n)*cos(lambda*(x-l(m-1))/sqrt(kappa(m)))));
            c(n) = c(n) + integral(prodt,l(m-1),lm);
        else
            prodt = @(x) ((u0(x)-(w(2*i-1) + w(2*i)*x)) .* ...
                (coeffs(2*i-1,n)*sin(lambda*(x-l(i-1))/sqrt(kappa(i))) + ...
                coeffs(2*i,n)*cos(lambda*(x-l(i-1))/sqrt(kappa(i)))));
            c(n) = c(n) + integral(prodt,l(i-1),l(i));
        end
    end
end

%% Time loop
tlength = length(tspan);
usoln = zeros((NX+1)*m,tlength); % No of unknowns
for cnt = 1:tlength
    
    t = tspan(cnt);
    
    % Solution (at given time)
    u = zeros(NX+1,m);
    
    for i = 1:m
        u(:,i) = w(2*i-1) + w(2*i)*xgrid(:,i);
    end

    % Form the sums corresponding to initial condition
    for n = 1:N
        lambda = eigvals(n);
        expc = exp(-t*lambda^2) * c(n);
        for i = 1:m
            if i == 1
                u(:,i) = u(:,i) + expc ...
                    * (coeffs(2*i-1,n)*sin(lambda*(xgrid(:,i)-l0)/sqrt(kappa(i))) + ...
                    coeffs(2*i,n)*cos(lambda*(xgrid(:,i)-l0)/sqrt(kappa(i))))/Nscal(n);
            elseif i == m
                u(:,i) = u(:,i) + expc ...
                    * (coeffs(2*i-1,n)*sin(lambda*(xgrid(:,i)-l(m-1))/sqrt(kappa(i))) + ...
                    coeffs(2*i,n)*cos(lambda*(xgrid(:,i)-l(m-1))/sqrt(kappa(i))))/Nscal(n);
            else
                u(:,i) = u(:,i) + expc ...
                    * (coeffs(2*i-1,n)*sin(lambda*(xgrid(:,i)-l(i-1))/sqrt(kappa(i))) + ...
                    coeffs(2*i,n)*cos(lambda*(xgrid(:,i)-l(i-1))/sqrt(kappa(i))))/Nscal(n);
            end
        end
    end
    
    usoln(:,cnt) = reshape(u,(NX+1)*m,1);
    
end

u = usoln;
x = reshape(xgrid,(NX+1)*m,1);

end

%% Matrix system arising from eigenvalue problem (want non-trivial solutions)
function A = Amat(lambda,m,kappa,l0,lm,l,aL,bL,aR,bR,H)

A(1,1) = sin(lambda*(l(1)-l0)/sqrt(kappa(1)))+...
    sqrt(kappa(1))*cos(lambda*(l(1)-l0)/sqrt(kappa(1)))*lambda/H(1);
A(1,2) = cos(lambda*(l(1)-l0)/sqrt(kappa(1)))-...
    sqrt(kappa(1))*sin(lambda*(l(1)-l0)/sqrt(kappa(1)))*lambda/H(1);
A(1,4) = -1;

A(2,1) = sqrt(kappa(1))*cos(lambda*(l(1)-l0)/sqrt(kappa(1)))*lambda;
A(2,2) = -sqrt(kappa(1))*sin(lambda*(l(1)-l0)/sqrt(kappa(1)))*lambda;
A(2,3) = -sqrt(kappa(2))*lambda;

for i = 2:m-1
    A(2*i-1,2*i-1) = sin(lambda*(l(i)-l(i-1))/sqrt(kappa(i)))+...
        sqrt(kappa(i))*cos(lambda*(l(i)-l(i-1))/sqrt(kappa(i)))*lambda/H(i);
    A(2*i-1,2*i)   = cos(lambda*(l(i)-l(i-1))/sqrt(kappa(i)))-...
        sqrt(kappa(i))*sin(lambda*(l(i)-l(i-1))/sqrt(kappa(i)))*lambda/H(i);
    A(2*i-1,2*i+1) = 0;
    A(2*i-1,2*i+2) = -1;
    
    A(2*i,2*i-1) = sqrt(kappa(i))*cos(lambda*(l(i)-l(i-1))/sqrt(kappa(i)))*lambda;
    A(2*i,2*i)   = -sqrt(kappa(i))*sin(lambda*(l(i)-l(i-1))/sqrt(kappa(i)))*lambda;
    A(2*i,2*i+1) = -sqrt(kappa(i+1))*lambda;
end

A(2*m-1,2*m-1) = bR*cos(lambda*(lm-l(m-1))/sqrt(kappa(m)))*...
    lambda/sqrt(kappa(m))+aR*sin(lambda*(lm-l(m-1))/sqrt(kappa(m)));
A(2*m-1,2*m)   = aR*cos(lambda*(lm-l(m-1))/sqrt(kappa(m)))...
    -bR*sin(lambda*(lm-l(m-1))/sqrt(kappa(m)))*lambda/sqrt(kappa(m));

A(2*m,1) = bL*lambda/sqrt(kappa(1));
A(2*m,2) = aL;

A = A / (lambda+1); % Scaling

end



