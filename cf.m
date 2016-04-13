function [zk,ck] = cf(n)
% Computes the poles zk and residues ck by the Caratheodory-Fejer method
% for the type (n,n) best approximation to exp(z) on the negative real 
% line.
%
% For full details see:
% J. A. C. Weideman L. N. Trefethen and T. Schmelzer (2006) Talbot
% quadratures and rational approximations. BIT Numer. Math., 46:653-670.
%
% The following code is given in Figure 4.1 of the above paper.

K = 75;                     % no of Cheb coeffs
nf = 1024;                  % no of pts for FFT
w = exp(2i*pi*(0:nf-1)/nf); % roots of unity
t = real(w);                % Cheb pts (twice over)
scl = 9;                    % scale factor for stability
F = exp(scl*(t-1)./(t+1+1e-16)); % exp(x) transpl. to [-1,1]
c = real(fft(F))/nf;        % Cheb coeffs of F
f = polyval(c(K+1:-1:1),w); % analytic part f of F
[U,S,V] = svd(hankel(c(2:K+1))); % SVD of Hankel matrix
s = S(n+1,n+1);             % singular value
u = U(K:-1:1,n+1)'; v = V(:,n+1)'; % singular vector
zz = zeros(1,nf-K);         % zeros for padding
b = fft([u zz])./fft([v zz]); % finite Blaschke product
rt = f-s*w.^K.*b;           % extended function r-tilde
rtc = real(fft(rt))/nf;     % its Laurent coeffs
zr = roots(v); qk = zr(abs(zr)>1); % poles
qc = poly(qk);              % coeffs of denominator
pt = rt.*polyval(qc,w);     % numerator
ptc = real(fft(pt)/nf);     % coeffs of numerator
ptc = ptc(n+1:-1:1); ck = 0*qk;

for k = 1:n                 % calculate residues
    q = qk(k); q2 = poly(qk(qk~=q));
    ck(k) = polyval(ptc,q)/polyval(q2,q);
end

zk = scl*(qk-1).^2./(qk+1).^2; % poles in z-plane
ck = 4*ck.*zk./(qk.^2-1);   % residues in z-plane
