close all, clear all, clc

% Add directory to current path
addpath('..')

% Parameters
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

% Plot
figure;
for i = 1:m-1, 
    plot([l(i),l(i)],[-0.1,1.1],'Color',[0.9,0.9,0.9])
    hold on
end
plot(x,u,'b','LineWidth',2.0)
axis([0,1,-0.1,1.1])
xlabel('$x$','Interpreter','LaTeX','FontSize',20)
ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
set(gca,'FontSize',14,'Layer','top')