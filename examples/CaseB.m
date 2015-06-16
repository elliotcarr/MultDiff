close all, clear all, clc

addpath('..')

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

% Plot
for i = 1:m-1, 
    plot([l(i),l(i)],[-0.1,1.1],'Color',[0.9,0.9,0.9])
    hold on
end
plot(x,u,'b','LineWidth',2.0)
axis([0,1,-0,1.1])
xlabel('$x$','Interpreter','LaTeX','FontSize',20)
ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
set(gca,'FontSize',14,'Layer','top')