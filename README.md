# MultDiff

MultDiff is a Matlab code for solving the one-dimensional multilayer diffusion problem.

# Examples

## Case A: 3 Layers

```
m     = 3;                       % Number of layers
l0    = 0.0;                     % Left end of slab
lm    = 1.0;                     % Right end of slab
l     = [0.3,0.7];               % Location of interfaces
kappa = [1,0.1,1];               % Diffusivities 
tvec  = [0.05,0.1,0.2,0.5,1.0];  % Compute solution at these values of t
u0    = @(x) zeros(size(x));     % Initial condition

% Boundary condition at x = l0
bconds.Ltype = 'Dirichlet'; 
bconds.aL    = 1.0; 
bconds.bL    = 0.0; 
bconds.cL    = 1.0;

% Boundary condition at x = lm
bconds.Rtype = 'Dirichlet'; 
bconds.aR    = 1.0; 
bconds.bR    = 0.0; 
bconds.cR    = 0.5;

[u,x] = multilayer_diffusion(m, kappa, l0, lm, l, bconds, u0, tvec);

plot(x,u,'b','LineWidth',2.0)
xlabel('$x$','Interpreter','LaTeX','FontSize',20)
ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
set(gca,'Xtick',[l0,l,lm],'Xgrid','on','FontSize',14)
set([gca,h1,h2],'FontSize',fs);
```

![Case A](https://github.com/elliotcarr/MultDiff/raw/master/figures/exampleA.png)
<figure><img src="https://github.com/elliotcarr/MultDiff/raw/master/figures/exampleA.png"></figure>

## Case B: 20 Layers

```
m     = 20;                            % Number of layers
l0    = 0.0;                           % Left end of slab
lm    = 1.0;                           % Right end of slab
dx    = (lm-l0)/m;
l     = dx:dx:lm-dx;                   % Location of interfaces
kappa = repmat([1.8,0.2],1,m/2);       % Diffusivities 
tvec  = [0.001,0.1,0.5,1.0,2.0,10.0];  % Compute solution at these values of t
u0    = @(x) ones(size(x));            % Initial condition

% Boundary condition at x = l0
bconds.Ltype = 'Neumann'; 
bconds.aL    = 0.0; 
bconds.bL    = 1.0; 
bconds.cL    = 0.0;

% Boundary condition at x = lm
bconds.Rtype = 'Dirichlet'; 
bconds.aR    = 1.0; 
bconds.bR    = 0.0; 
bconds.cR    = 0.0;

[u,x] = multilayer_diffusion(m, kappa, l0, lm, l, bconds, u0, tvec);

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

![Case B](https://github.com/elliotcarr/MultDiff/raw/master/figures/exampleB.png)

## Case C: Imperfect contact at the interfaces

## Case D: 1000 Layers

# Installation

MultDiff can be downloaded from

https://github.com/elliotcarr/MultDiff/archive/master.zip

After unzipping, you will need to add MultDiff to the MATLAB path. You can do
this by adding the commands:
```
addpath(multdiffroot)
```
where `multdiffroot` is the path to the unzipped directory.

# License

See `LICENSE.txt` for MultDiff's licensing information.
