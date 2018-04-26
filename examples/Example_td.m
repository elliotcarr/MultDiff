close all, clear all, clc

% Add directory to current path
addpath('..')

% Parameters
m     = 3;                            % Number of layers
gamma = [1,0.1,1];                    
D     = [1,0.1,1];                    % Diffusivities 
l0    = 0.0;                          % Left end of slab
lm    = 1.0;                          % Right end of slab
l     = [0.3,0.7];                    % Location of interfaces
u0    = @(x,i) zeros(size(x));        % Initial condition
% Boundary conditions require Laplace transform of boundary functions
a = 1; b = 0.5; c = 2;
Lbnd  = {'Dirichlet',1.0,0.0,@(t)a + b*cos(c*t),@(s)a/s + b*s / (s^2 + c^2)}; % Boundary condition (x = l0)
Rbnd  = {'Neumann',0.0,1.0,@(t)0.0,@(s)0.0/s}; % Boundary condition (x = lm)
tspan = [0.02,0.05,0.1,0.2,0.5,1.0];  % Times at which to compute solution
H     = Inf*ones(m-1,1);
theta = ones(m-1,1);
[u,x] = multdiff_td(m,D,gamma,theta,l0,lm,l,u0,Lbnd,Rbnd,tspan,H);

% Plot
figure;
for i = 1:m-1, 
    plot([l(i),l(i)],[0.0,2.1],'Color',[0.9,0.9,0.9])
    hold on
end
for k = 1:length(tspan)
    plot(x,u(:,:,k),'b','LineWidth',2.0)
end
axis([0,1,0.0,2.1])
xlabel('$x$','Interpreter','LaTeX','FontSize',20)
ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
set(gca,'FontSize',14,'Layer','top')