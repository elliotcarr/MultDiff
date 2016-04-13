function eigs = eigvals(a,b,c,N)
% Solves the nonlinear equation:
%
% f(lambda) = lambda * tan(lambda*(a-b)) = c
%
% numerically for the first N non-negative values of lambda, using a
% combination of the bisection method and Newton's method.
%
% Solution to this equation is required to compute the eigenvalues when
% Robin conditions are applied at either the left or right boundaries.

f = @(lambda) lambda * sin(lambda*(a-b)) - c*cos(lambda*(a-b));
for i = 1:10
    g = chebfun(f,[0,10^i]);
    r = roots(g);
    numroots = length(r);
    if numroots >= N
        r = r(1:N);
        break;
    end
end
eigs = r;

% pause;
% 
% MaxIters = 20;
% tol      = 1e-8;
% eigs     = zeros(N,1);
% 
% if c > 0
%     n = 0;
%     converged = false;
%     
%     lambda_left  = 0.0;
%     lambda_right = (2*n+1) * pi / (2 * (a-b));
%     
%     % Initial guess (midpoint)
%     lambda = (lambda_left+lambda_right)/2;
%     
%     f = compute_func(lambda,a,b,c);
%     absf0 = abs(f);
%     
%     for k = 1:MaxIters
%         fdash = compute_fdash(lambda,a,b,c);
%         lambda = lambda - f / fdash;
%         f = compute_func(lambda,a,b,c);
%         if abs(f)/absf0 < tol
%             converged = true;
%             break;
%         end
%     end
%     
%     if converged == true
%         eigs(1) = lambda;
%     else
%         warning('Failed to converge when finding eigenvalue (n=%i)\n',n+1);
%     end
%     
% end
% 
% nu = 0;
% n = 0;
% while nu < N
%     
%     n = n + 1;
%     
%     converged = false;
%     
%     if c < 0
%         lambda_left  = max((2*n+1) * pi / (2 * (a-b)),0.0);
%         lambda_right = (2*(n+1)+1) * pi / (2 * (a-b));
%     else
%         lambda_left  = max((2*(n-1)+1) * pi / (2 * (a-b)),0.0);
%         lambda_right = (2*n+1) * pi / (2 * (a-b));
%     end
%     
%     % Initial guess (midpoint)
%     lambda = (lambda_left + lambda_right)/2;
%     
%     f = compute_func(lambda,a,b,c);
%     for k = 1:MaxIters
%         fdash = compute_fdash(lambda,a,b,c);
%         lambda = lambda - f / fdash;
%         f = compute_func(lambda,a,b,c);
%         if abs(f) < tol
%             converged = true;
%             break;
%         end
%     end
%     
%     if converged == true
%         unique_eig = true;
%         for i = 0:nu-1
%             if abs(lambda - eigs(i+1)) < 1e-4*eigs(i+1)
%                 unique_eig = false;
%                 break;
%             end
%         end
%         if unique_eig
%             eigs(nu+1) = lambda;
%             nu = nu + 1; % Unique eigenvalue counter
%             if nu == N
%                 break;
%             end
%         end
%     else
%         warning('Failed to converge when finding eigenvalue (n=%i)\n',n+1);
%     end
%     
% end
% 
% eigs = sort(eigs);
% 
% end
% 
% % Sub-functions
% % Compute function f(lambda)
% function f = compute_func(lambda,a,b,c)
% f = lambda * sin(lambda*(a-b)) - c*cos(lambda*(a-b));
% end
% % Compute derivative f'(lambda)
% function fdash = compute_fdash(lambda,a,b,c)
% fdash = sin(lambda*(a-b)) + ...
%     lambda*cos(lambda*(a-b))*(a-b)+c*sin(lambda*(a-b))*(a-b);
% end
