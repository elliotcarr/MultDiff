close all, clear all, clc

addpath('..')

% Parameters
m     = 3;                           % Number of layers
l0    = 0.0;                         % Left end of slab
lm    = 1.0;                         % Right end of slab
l     = [0.3,0.7];                   % Location of interfaces
kappa = [1,0.1,1];                   % Diffusivities 
tspan = [0.02,0.05,0.1,0.2,0.5,1.0]; % Times at which to compute solution

% Initial condition
u0 = @(x) zeros(size(x));     

% Boundary conditions
bcs.Ltype = 'Dirichlet'; bcs.aL = 1.0; bcs.bL = 0.0; bcs.cL = 1.0;
bcs.Rtype = 'Dirichlet'; bcs.aR = 1.0; bcs.bR = 0.0; bcs.cR = 0.5;

% Solve
[u,x] = multdiff(m,kappa,l0,lm,l,u0,bcs,tspan,'Perfect');

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