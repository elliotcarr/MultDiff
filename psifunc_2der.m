function psi_func = psifunc_2der(i,j,Lbnd,Rbnd,x,m,l0,lm,l,gamma)
% Defines the second spatial derivative of the psi functions that homogeneous 
% the boundary conditions

Ltype = Lbnd{1};
Rtype = Rbnd{1};
bL    = Lbnd{3};
bR    = Rbnd{3};

if i == 1
    switch Ltype
        case 'Neumann'
            if j == 1
                psi_func = 1.0 / (bL*(l(1)-l0));
            elseif j == 2
                psi_func = 1.0 / (gamma(1)*(l(1)-l0));
            end 
        otherwise
            if j == 1
                psi_func = 0.0;
            elseif j == 2
                psi_func = 0.0;
            end            
    end
elseif i == m
    switch Rtype
        case 'Neumann'
            if j == 1
                psi_func = 1.0 / (gamma(m)*(l(m-1)-lm));
            elseif j == 2
                psi_func = 1.0 / (bR*(lm-l(m-1)));
            end 
        otherwise
            if j == 1
                psi_func = 0.0;
            elseif j == 2
                psi_func = 0.0;
            end 
    end
else
    if j == 1
        psi_func = 1.0 / (gamma(i)*(l(i-1)-l(i)));
    elseif j == 2
        psi_func = 1.0 / (gamma(i)*(l(i)-l(i-1)));
    end
end

