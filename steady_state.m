function w = steady_state(m,kappa,l0,lm,l,Lbnd,Rbnd)

aL = Lbnd{2};
bL = Lbnd{3};
cL = Lbnd{4};
aR = Rbnd{2};
bR = Rbnd{3};
cR = Rbnd{4};

H = 1e16*ones(m-1,1);

A = zeros(2*m,2*m);
b = zeros(2*m,1);

% Interface conditions
for i = 1:m-1
    A(2*i-1,2*i-1) = 1.0;
    A(2*i-1,2*i)   = (kappa(i)/H(i)) + l(i);
    A(2*i-1,2*i+1) = -1.0;
    A(2*i-1,2*i+2) = -l(i);
    b(2*i-1)       = 0.0;
    
    A(2*i,2*i)   = kappa(i);
    A(2*i,2*i+2) = -kappa(i+1);
    b(2*i)       = 0.0;
end

% Right boundary condition
A(2*m-1,2*m)   = bR+aR*lm;
A(2*m-1,2*m-1) = aR;
b(2*m-1)       = cR;

% Left boundary condition
A(2*m,1) = aL;
A(2*m,2) = -bL+aL*l0;
b(2*m)   = cL;

w = A\b;