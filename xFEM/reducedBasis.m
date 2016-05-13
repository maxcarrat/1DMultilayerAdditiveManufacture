function [ reduceBasis, uNormalized ] = reducedBasis( T, dofs )
%% Reduce Order Basis
% generate a reduce order basis functions, from J.G.Michopulous
% "Performance of Reduced Order Models of Moving Heat Source Deposition
% Problems for efficient inverse analysis", 2014
%T = matrix of solution snapshots
%dofs = degrees of freedom of the local problem.

Y = T(:,2:end)'*T(:,2:end);

% single value decomposition, return D and V are the matrices of the
% eigenvectors and S is the diagonal matrix of the eigenvalues

[D, S, ~] = svd(Y);

k = extractK(S, 10e-2);                         % number of highly energetic eigenvalues
u = zeros(dofs, k);                             % basis vectors matrix

for i=1:k
    for j=1:k
         sum = D(j,i)*T(:,j);
    end
    u(:,i) = 1/sqrt(S(i,i)) * sum;
end

% normalization of the basis vector
uNorm = norm(u);
uNormalized = u/uNorm;

% new vector of reduced order basis functions
reduceBasis =@(x, i) problem.basis_fun(x, i, 0)*uNormalized;

end


function [ k ] = extractK(S, toll)
    k = [];
    for i=1:size(S)
        if(abs(S(i,i)/max(max(S))) >= toll)
            k = [k i];
        end
    end
end

