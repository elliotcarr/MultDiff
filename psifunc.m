function psi_func = psifunc(i,j,Lbnd,Rbnd,x,m,l0,lm,l,gamma)
% Defines the psi functions that homogeneous the boundary conditions

Ltype = Lbnd{1};
Rtype = Rbnd{1};
aL    = Lbnd{2};
bL    = Lbnd{3};
aR    = Rbnd{2};
bR    = Rbnd{3};

if i == 1
    switch Ltype
        case 'Neumann'
            if j == 1
                psi_func = x.*(x-2*l(1)) / (2*bL*(l(1)-l0));
            elseif j == 2
                psi_func = x.*(x-2*l0) / (2*gamma(1)*(l(1)-l0));
            end 
        otherwise
            if j == 1
                psi_func = 1/aL;
            elseif j == 2
                psi_func = (x-l0+bL/aL) / gamma(1);
            end            
    end
elseif i == m
    switch Rtype
        case 'Neumann'
            if j == 1
                psi_func = x.*(x-2*lm) / (2*gamma(m)*(l(m-1)-lm));
            elseif j == 2
                psi_func = x.*(x-2*l(m-1)) / (2*bR*(lm-l(m-1)));
            end 
        otherwise
            if j == 1
                psi_func = (x-lm-bR/aR) / gamma(m);
            elseif j == 2
                psi_func = 1/aR;
            end 
    end
else
    if j == 1
        psi_func = x.*(x-2*l(i)) / (2*gamma(i)*(l(i-1)-l(i)));
    elseif j == 2
        psi_func = x.*(x-2*l(i-1)) / (2*gamma(i)*(l(i)-l(i-1)));
    end
end

