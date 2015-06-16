close all, clear all, clc

% Parameters
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

% Solve
[u,x] = multilayer_diffusion(m, kappa, l0, lm, l, bconds, u0, tvec);

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