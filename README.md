## MultDiff: A Matlab code for solving the one-dimensional multilayer diffusion problem

``MultDiff`` solves the transient diffusion equation in a one-dimensional composite slab of finite length consisting of `m` layers. The code is applicable to both perfect and imperfect contact at the interfaces between adjacent layers and either Dirichlet, Neumann or Robin boundary conditions at the ends of the slab. The code works well for a large number of layers (large `m`).

## References

If you use ``MultDiff``, we would appreciate that you mention it in your work by citing the following paper:

TBA

<!--- Most approaches for this problem require the solution of a complex transcendental equation arising from the determinant of a `2m by 2m` matrix for the eigenvalues, which is difficult to solve for large `m`. Our approach is based on a semi-analytic method based on the Laplace transform and an orthogonal eigenfunction expansion involving eigenvalues that are obtained either explicitly or by solving simple transcendental equations. -->

## Examples

``MultDiff`` is best demonstrated with some examples.

### Case A

```
% Parameters
m     = 3;                       % Number of layers
l0    = 0.0;                     % Left end of slab
lm    = 1.0;                     % Right end of slab
l     = [0.3,0.7];               % Location of interfaces
kappa = [1,0.1,1];               % Diffusivities 
tvec  = [0.05,0.1,0.2,0.5,1.0];  % Compute solution at these values of t

% Initial condition
u0 = @(x) zeros(size(x));     

% Boundary conditions
bcs.Ltype = 'Dirichlet'; bcs.aL = 1.0; bcs.bL = 0.0; bcs.cL = 1.0;
bcs.Rtype = 'Dirichlet'; bcs.aR = 1.0; bcs.bR = 0.0; bcs.cR = 0.5;

% Solve
[u,x] = multilayer_diffusion(m, kappa, l0, lm, l, u0, bcs, tvec);
```

We can plot the solution to the above problem using the following code:

```
for i = 1:m-1, 
    plot([l(i),l(i)],[-0.1,1.1],'Color',[0.9,0.9,0.9])
    hold on
end
plot(x,u,'b','LineWidth',2.0)
axis([0,1,-0.1,1.1])
xlabel('$x$','Interpreter','LaTeX','FontSize',20)
ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
set(gca,'FontSize',14,'Layer','top')
```

<figure><img src="https://github.com/elliotcarr/MultDiff/raw/master/figures/CaseA.png"></figure>

### Case B

```
% Parameters
m     = 20;                            % Number of layers
l0    = 0.0;                           % Left end of slab
lm    = 1.0;                           % Right end of slab
dx    = (lm-l0)/m; l = dx:dx:lm-dx;    % Location of interfaces
kappa = repmat([1.0,0.1],1,m/2);       % Diffusivities 
tvec  = [0.001,0.1,0.5,1.0,2.0,15.0];  % Compute solution at these values of t

% Initial condition
u0 = @(x) ones(size(x));            

% Boundary conditions
bcs.Ltype = 'Neumann'; bcs.aL = 0.0; bcs.bL = 1.0; bcs.cL = 0.0;
bcs.Rtype = 'Dirichlet'; bcs.aR = 1.0; bcs.bR = 0.0; bcs.cR = 0.1;

% Solve
[u,x] = multilayer_diffusion(m, kappa, l0, lm, l, u0, bcs, tvec);
```

<figure><img src="https://github.com/elliotcarr/MultDiff/raw/master/figures/CaseB.png"></figure>


### Case C

`MultDiff` also allows for imperfect contact at the interfaces between adjacent layers.

```
% Parameters
m     = 20;                          % Number of layers
l0    = 0.0;                         % Left end of slab
lm    = 1.0;                         % Right end of slab
dx    = (lm-l0)/m; l = dx:dx:lm-dx;  % Location of interfaces
kappa = ones(m,1);                   % Diffusivities 
tvec  = [0.01,0.1,0.2,0.5,1.0];      % Compute solution at these values of t

% Initial condition
u0 = @(x) zeros(size(x));     

% Boundary conditions
bcs.Ltype = 'Dirichlet'; bcs.aL = 1.0; bcs.bL = 0.0; bcs.cL = 1.0;
bcs.Rtype = 'Dirichlet'; bcs.aR = 1.0; bcs.bR = 0.0; bcs.cR = 1.0;

% Contact transfer coefficients at interfaces
options.H = 30*ones(m-1,1);

% Solve
[u,x] = multilayer_diffusion(m, kappa, l0, lm, l, u0, bcs, tvec, options);
```

<figure><img src="https://github.com/elliotcarr/MultDiff/raw/master/figures/CaseC.png"></figure>

### Case D

`MultDiff` can can also be used to solve the single layer problem:

```
% Parameters
m     = 3;                   % Number of layers
l0    = 0.0;                 % Left end of slab
lm    = 1.0;                 % Right end of slab
l     = [0.3,0.7];           % Location of interfaces
kappa = [1,1,1];             % Diffusivities 
tvec  = [0.05,0.1,0.2,1.0];  % Compute solution at these values of t

% Initial condition
u0 = @(x) zeros(size(x));     

% Boundary conditions
bcs.Ltype = 'Dirichlet'; bcs.aL = 1.0; bcs.bL = 0.0; bcs.cL = 1.0;
bcs.Rtype = 'Dirichlet'; bcs.aR = 1.0; bcs.bR = 0.0; bcs.cR = 0.5;

% Solve
[u,x] = multilayer_diffusion(m, kappa, l0, lm, l, u0, bcs, tvec);
```

<figure><img src="https://github.com/elliotcarr/MultDiff/raw/master/figures/CaseD.png"></figure>

## Installation

``MultDiff`` can be downloaded from

https://github.com/elliotcarr/MultDiff/archive/master.zip

After unzipping, you will need to add the directory to the MATLAB path. You can do
this via the command:
```
addpath(multdiffroot)
```
where `multdiffroot` is the path to the unzipped directory.

## License

See `LICENSE.txt` for licensing information.
