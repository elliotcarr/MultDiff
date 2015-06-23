close all, clear all, clc

% Add directory to current path
addpath('..')

% Parameters
m     = 3;                            % Number of layers
kappa = [1,0.1,1];                    % Diffusivities 
l0    = 0.0;                          % Left end of slab
lm    = 1.0;                          % Right end of slab
l     = [0.3,0.7];                    % Location of interfaces
u0    = @(x) zeros(size(x));          % Initial condition
Lbnd  = {'Dirichlet',1.0,0.0,1.0};    % Boundary condition (x = l0)
Rbnd  = {'Dirichlet',1.0,0.0,0.5};    % Boundary condition (x = lm)
tspan = [0.02,0.05,0.1,0.2,0.5,1.0];  % Compute solution at these values of t

% Compute Semi-Analytic solution
[u,x] = multdiff(m,kappa,l0,lm,l,u0,Lbnd,Rbnd,tspan,'Perfect');

% Compute FVM solution
dt = 0.001;       % Time step
options.NX = 10;  % Number of divisions within each slab 
[uf,xf] = multdiff_fvm(m,kappa,l0,lm,l,u0,Lbnd,Rbnd,tspan,'Perfect',dt,options);

% Plot
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