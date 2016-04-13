function wfunc = wfunc(i,Ltype,Rtype,x,m,w)
% Form of the non-transient component of the solution 
% (steady state if it exists)

% if strcmp(Ltype,'Neumann') && strcmp(Rtype,'Neumann')
%     if i == 1
%         wfunc = w(1)*x + w(2)*x.^2;
%     elseif i == m
%         wfunc = w(2*m-1)*x + w(2*m)*x.^2;
%     else
%         wfunc = w(2*i-1)*x + w(2*i)*x.^2;
%     end
% else
    wfunc = w(2*i-1) + w(2*i)*x;
% end