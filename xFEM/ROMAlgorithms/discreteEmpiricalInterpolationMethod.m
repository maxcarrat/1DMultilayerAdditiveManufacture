function [U, P, p] = discreteEmpiricalInterpolationMethod(V, m) 
%DISCRETEEMPIRICALINTERPOLATIONMETHOD evaluate the DEIM indeces and
%reduction opertaor. See S. Chaturantabut and D. C. Sorensen. 
%Nonlinear model reduction via discrete empirical interpolation, 2011.
%   V = DEIM snapshot matrix
%   m = number of DEIM indeces we want to evaluate

[~, ip] = max(abs(V(:,1)));
U = V(:,1);
ep = zeros(size(V,1), 1);
ep(ip) = 1.0;

P = ep;
p = ip;

for l=2:m
    c = (P'*U)\(P'*V(:,l));
    r = V(:,l) - U*c;
    [~, ip] = max(abs(r));
    U = [U, V(:,l)];
    ep = zeros(size(V,1), 1);
    ep(ip) = 1.0;
    P = [P, ep];
    p = [p; ip];
end

end

