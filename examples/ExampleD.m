close all, clear all, clc

% Add directory to current path
addpath('..')

% Parameters
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

% Plot
figure;
for i = 1:m-1, 
    plot([l(i),l(i)],[-0.1,1.1],'Color',[0.9,0.9,0.9])
    hold on
end
plot(x,u,'b-','LineWidth',2.0)
axis([0,1,-0.1,1.1])
xlabel('$x$','Interpreter','LaTeX','FontSize',20)
ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
set(gca,'FontSize',14,'Layer','top')