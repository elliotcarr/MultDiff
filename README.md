## MultDiff: A Matlab code for solving the one-dimensional multilayer diffusion problem

``MultDiff`` solves the transient diffusion equation in a one-dimensional composite slab of finite length consisting of `m` layers. The code is applicable to both perfect and imperfect contact at the interfaces between adjacent layers and either Dirichlet, Neumann or Robin boundary conditions at the ends of the slab. The code works well for a large number of layers (large `m`).

The code features two solution approaches:
 - ``multdiff`` uses a semi-analytic method based on the Laplace transform and an orthogonal eigenfunction expansion.
 - ``multdiff_fvm`` uses a finite volume method and a backward Euler time discretisation.

## References

If you use ``MultDiff``, we would appreciate that you mention it in your work by citing the following paper:

*To be announced*

<!--- Most approaches for this problem require the solution of a complex transcendental equation arising from the determinant of a `2m by 2m` matrix for the eigenvalues, which is difficult to solve for large `m`. Our approach is based on a semi-analytic method based on the Laplace transform and an orthogonal eigenfunction expansion involving eigenvalues that are obtained either explicitly or by solving simple transcendental equations. -->

## Examples

``MultDiff`` is best demonstrated with some examples.

### Example A

``MultDiff`` can be used to solve multilayer diffusion problems with three or more layers.

```
m     = 3;                            % Number of layers
kappa = [1,0.1,1];                    % Diffusivities 
l0    = 0.0;                          % Left end of slab
lm    = 1.0;                          % Right end of slab
l     = [0.3,0.7];                    % Location of interfaces
u0    = @(x) zeros(size(x));          % Initial condition
Lbnd  = {'Dirichlet',1.0,0.0,1.0};    % Boundary condition (x = l0)
Rbnd  = {'Dirichlet',1.0,0.0,0.5};    % Boundary condition (x = lm)
tspan = [0.02,0.05,0.1,0.2,0.5,1.0];  % Times at which to compute solution
[u,x] = multdiff(m,kappa,l0,lm,l,u0,Lbnd,Rbnd,tspan,'Perfect');
```

We can plot the solution to the above problem using the following code:

```
% Plot
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

<figure><img src="https://github.com/elliotcarr/MultDiff/raw/master/figures/ExampleA.png"></figure>

The solution can be verified against the Finite Volume solution:

```
dt = 0.001;       % Time step
options.NX = 10;  % Number of divisions within each slab 
[uf,xf] = multdiff_fvm(m,kappa,l0,lm,l,u0,Lbnd,Rbnd,tspan,'Perfect',dt,options);

figure;
for i = 1:m-1, 
    plot([l(i),l(i)],[-0.1,1.1],'Color',[0.9,0.9,0.9])
    hold on
end
p1 = plot(x,u,'b-','LineWidth',2.0);
p2 = plot(xf,uf,'b.','MarkerSize',24);
axis([0,1,-0.1,1.1])
xlabel('$x$','Interpreter','LaTeX','FontSize',20)
ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
set(gca,'FontSize',14,'Layer','top')
h = legend([p1(1),p2(1)],'Semi-Analytic','Finite Volume','Location','NorthOutside','Orientation','Horizontal');
set(h,'Interpreter','LaTex','FontSize',14)
```

<figure><img src="https://github.com/elliotcarr/MultDiff/raw/master/figures/ExampleA_fvm.png"></figure>

### Example B

``MultDiff`` can be used to solve multilayer diffusion problems with many layers.

```
m     = 100;                           % Number of layers
kappa = 1.1 + cos(1:m);                % Diffusivities 
l0    = 0.0;                           % Left end of slab
lm    = 1.0;                           % Right end of slab
dx    = (lm-l0)/m; l = dx:dx:lm-dx;    % Location of interfaces
u0    = @(x) ones(size(x));            % Initial condition
Lbnd  = {'Neumann',0.0,1.0,0.0};       % Boundary condition (x = l0)
Rbnd  = {'Dirichlet',1.0,0.0,0.1};     % Boundary condition (x = lm)
tspan = [0.001,0.1,0.5,1.0,2.0,15.0];  % Times at which to compute solution
[u,x] = multdiff(m,kappa,l0,lm,l,u0,Lbnd,Rbnd,tspan,'Perfect');
```

<figure><img src="https://github.com/elliotcarr/MultDiff/raw/master/figures/ExampleB.png"></figure>


### Example C

`MultDiff` also allows for imperfect contact at the interfaces between adjacent layers.

```
m     = 20;                          % Number of layers
kappa = ones(m,1);                   % Diffusivities 
l0    = 0.0;                         % Left end of slab
lm    = 1.0;                         % Right end of slab
dx    = (lm-l0)/m; l = dx:dx:lm-dx;  % Location of interfaces
u0    = @(x) zeros(size(x));         % Initial condition
Lbnd  = {'Dirichlet',1.0,0.0,1.0};   % Boundary condition (x = l0)
Rbnd  = {'Dirichlet',1.0,0.0,1.0};   % Boundary condition (x = lm)
tspan = [0.01,0.1,0.2,0.5,1.0];      % Times at which to compute solution
H     = 30*ones(m-1,1);              % Contact transfer coefficients at interfaces
[u,x] = multdiff(m,kappa,l0,lm,l,u0,Lbnd,Rbnd,tspan,'Imperfect',H);
```

<figure><img src="https://github.com/elliotcarr/MultDiff/raw/master/figures/ExampleC.png"></figure>

### Example D

`MultDiff` can can also be used to solve the single layer problem.

```
m     = 3;                          % Number of layers
kappa = [1,1,1];                    % Diffusivities 
l0    = 0.0;                        % Left end of slab
lm    = 1.0;                        % Right end of slab
l     = [0.3,0.7];                  % Location of interfaces
u0    = @(x) zeros(size(x));        % Initial condition
Lbnd  = {'Dirichlet',1.0,0.0,1.0};  % Boundary condition (x = l0)
Rbnd  = {'Robin',1.0,0.5,0.5};      % Boundary condition (x = lm)
tspan = [0.02,0.05,0.1,0.2,1.0];    % Times at which to compute solution
[u,x] = multdiff(m,kappa,l0,lm,l,u0,Lbnd,Rbnd,tspan,'Perfect');
```

<figure><img src="https://github.com/elliotcarr/MultDiff/raw/master/figures/ExampleD.png"></figure>

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

See `LICENSE.md` for licensing information.
