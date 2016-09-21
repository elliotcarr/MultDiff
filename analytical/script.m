clear all, clc, close all, commandwindow

% Add directory to chebfun package
% addpath(</path/to/chebfun>)

% Problem = '1';
% Problem = '2';
% Problem = '3';
Problem = '4';

%% Parameters

if strcmp(Problem,'1')
    
    m         = 4; % Number of layers
    kappa     = [0.2,0.01,0.1,1]; % Diffusivities
    l0        = 0.0; % Left end of slab
    lm        = 1.0; % Right end of slab
    l         = [0.25,0.5,0.75]; % Location of interfaces
    u0        = @(x) ones(size(x)); % Initial condition
    Lbnd      = {'Dirichlet',1.0,0.0,0.2}; % Boundary condition (x = l0)
    Rbnd      = {'Robin',0.4,kappa(4),0.04}; % Boundary condition (x = lm)
    tspan     = [0.02,0.1,0.5,1.0,2.0,10.0]; % Compute solution at these values of t
    interface = 'Perfect';
    
elseif strcmp(Problem,'2')
    
    m         = 10; % Number of layers
    kappa     = repmat([1,0.1],1,m/2); % Diffusivities
    l0        = 0.0; % Left end of slab
    lm        = 1.0; % Right end of slab
    l         = linspace(l0,lm,m+1); l([1,m+1]) = []; % Location of interfaces
    u0        = @(x) zeros(size(x)); % Initial condition
    Lbnd      = {'Dirichlet',1.0,0.0,1.0}; % Boundary condition (x = l0)
    Rbnd      = {'Dirichlet',1.0,0.0,0.0}; % Boundary condition (x = lm)
    interface = 'Perfect';
    tspan     = [0.01,0.2,1];  % Compute solution at these values of t
    
elseif strcmp(Problem,'3')
    
    m         = 10; % Number of layers
    kappa     = repmat([1,0.1],1,m/2); % Diffusivities
    l0        = 0.0; % Left end of slab
    lm        = 1.0; % Right end of slab
    l         = linspace(l0,lm,m+1); l([1,m+1]) = []; % Location of interfaces
    u0        = @(x) zeros(size(x)); % Initial condition
    Lbnd      = {'Dirichlet',1.0,0.0,1.0}; % Boundary condition (x = l0)
    Rbnd      = {'Neumann',0.0,1.0,0.0}; % Boundary condition (x = lm)
    tspan     = [0.007,2,10.0]; % Compute solution at these values of t
    interface = 'Imperfect';
    H         = 0.5*ones(m-1,1); % Contact transfer coefficients at interfaces
    
elseif strcmp(Problem,'4')
    
    m         = 4; % Number of layers
    kappa     = [1.0,0.1,1.0,0.1]; % Diffusivities
    l0        = 0.0; % Left end of slab
    lm        = 1.0; % Right end of slab
    l         = [0.25,0.5,0.75]; % Location of interfaces
    u0        = @(x) 1.0*(x<0.5) + 0.0; % Initial condition
    Lbnd      = {'Neumann',0.0,1.0,0.0}; % Boundary condition (x = l0)
    Rbnd      = {'Neumann',0.0,1.0,0.0}; % Boundary condition (x = lm)
    tspan     = [0.002,0.1,0.5,1.0]; % Compute solution at these values of t
    interface = 'Perfect';
    
end
 
%% Compute analytical solution
options.N = 100; % Number of eigenvalues
options.NX = 50; % Number of divisions in each slab
options.prob = 'C';
if strcmp(interface,'Perfect')
    [u,x] = multdiff_global(m,kappa,l0,lm,l,u0,Lbnd,Rbnd,tspan,interface,options);
elseif strcmp(interface,'Imperfect')
    [u,x] = multdiff_global(m,kappa,l0,lm,l,u0,Lbnd,Rbnd,tspan,interface,H,options);
end

%% Plot
figure;
for i = 1:m-1, 
    plot([l(i),l(i)],[-0.1,1.1],'Color',[0.9,0.9,0.9])
    hold on
end
p1 = plot(x,u,'b-','LineWidth',2.0);
axis([0,1,-0.1,1.1])
xlabel('$x$','Interpreter','LaTeX','FontSize',20)
ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
set(gca,'FontSize',14,'Layer','top')
h = legend([p1(1)],'Analytical','Location','NorthOutside','Orientation','Horizontal');
set(h,'Interpreter','LaTex','FontSize',14)
drawnow;