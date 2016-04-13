function eig_func = eigfunc(lambda,i,Ltype,Rtype,x,m,l0,lm,l)
% Defines the eigenfunctions for the different layers and boundary
% conditions. 

if i == 1
    switch Ltype
        case 'Dirichlet'
            eig_func = sin(lambda*(x-l0));
        case 'Neumann'
            eig_func = cos(lambda*(x-l0));
        case 'Robin'
            eig_func = cos(lambda*(l(1)-x));
    end
elseif i == m
    switch Rtype
        case 'Dirichlet'
            eig_func = sin(lambda*(lm - x));
        case 'Neumann'
            eig_func = cos(lambda*(lm - x));
        case 'Robin'
            eig_func = cos(lambda*(x-l(m-1)));
    end
else
    eig_func = cos(lambda*(l(i)-x));
end

