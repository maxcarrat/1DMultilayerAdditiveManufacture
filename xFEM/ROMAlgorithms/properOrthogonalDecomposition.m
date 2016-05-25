function [ uNormalized, numberOfModes ] = properOrthogonalDecomposition( T )
%% Proper Orthogonal decomposition
% generate a reduced operator, from J.G.Michopulous
% "Performance of Reduced Order Models of Moving Heat Source Deposition
% Problems for efficient inverse analysis", 2014
%T = matrix of solution snapshots

Y = T(:,1:end)'*T(:,1:end);

% single value decomposition, return D and V are the matrices of the
% eigenvectors and S is the diagonal matrix of the eigenvalues

[D, S, ~] = svd(Y);

k = extractK(S, 10e-10);                          % number of highly energetic eigenvalues
u = zeros(size(T,1), size(k,2));                  % basis vectors matrix

for i=1:size(k,2)
    for j=1:size(k,2)
         sum = D(j,i)*T(:,j);
    end
    u(:,i) = 1/sqrt(S(i,i)) * sum;
end

numberOfModes = size(k,2); 

% normalization of the basis vector
uNorm = norm(u);
uNormalized = u/uNorm;

end


function [ k ] = extractK(S, toll)
    k = [];
    for i=1:size(S)
        if(abs(S(i,i)/max(max(S))) >= toll)
            k = [k i];
        end
    end
end

