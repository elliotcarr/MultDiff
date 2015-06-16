## MultDiff

MultDiff is a Matlab code for solving the one-dimensional multilayer diffusion problem.

## Examples

##### Example A: Three layer slab (exampleA.m)

```matlab
m     = 3;                   % Number of layers
l0    = 0.0;                 % Left end of slab
lm    = 1.0;                 % Right end of slab
l     = [1/3,2/3];           % Location of interfaces
kappa = [1,0.1,0.2];         % Diffusivities 
H     = [0.5,0.5];           % Contact transfer coefficients at interfaces
tvec  = [0.007,2.0,10.0];    % Compute solution at these values of t
phi0  = @(x) zeros(size(x)); % Initial condition

% Boundary condition at x = l0
bconds.Ltype = 'Dirichlet'; bconds.aL = 1.0; bconds.bL = 0.0; bconds.cL = 1.0;

% Boundary condition at x = lm
bconds.Rtype = 'Dirichlet'; bconds.aR = 1.0; bconds.bR = 0.0; bconds.cR = 0.0;

u = multilayer_diffusion(m, kappa, l0, lm, l, bconds, H, tvec, phi0, options)

plot(x,u)
```

## Installation

MulDiff can be downloaded from

https://github.com/elliotcarr/MultDiff/archive/master.zip

After unzipping, you will need to add MultDiff to the MATLAB path. You can do
this by adding the commands:
```
addpath(multdiffroot)
```
where `multdiffroot` is the path to the unzipped directory.

## License

See `LICENSE.txt` for MultDiff's licensing information.
