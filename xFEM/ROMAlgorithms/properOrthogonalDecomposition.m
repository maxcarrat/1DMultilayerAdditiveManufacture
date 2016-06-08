function [ uNormalized, numberOfModes ] = properOrthogonalDecomposition( T, numberOfModes )
%% Proper Orthogonal decomposition
% generate a reduced operator, from J.G.Michopulous
% "Performance of Reduced Order Models of Moving Heat Source Deposition
% Problems for efficient inverse analysis", 2014
%T = matrix of solution snapshots

Y = T(:,:)'*T(:,:);

% eigendecomposition, return D as matrix of the
% eigenvectors and S as diagonal matrix of the eigenvalues.

% [D, S, ~] = svd(Y);
% 
% if numberOfModes==0
%     k = extractK(S, 10e-8);                          % number of highly energetic eigenvalues
%     u = zeros(size(T,1), size(k,2));                  % basis vectors matrix
%     numberOfModes = size(k,2);
% end
% 
% 
% for i=1:numberOfModes
%     for j=1:numberOfModes
%          sum = D(j,i)*T(:,j);
%     end
%     u(:,i) = 1/sqrt(S(i,i)) * sum;
% end
% 
% % normalization of the basis vector
% uNorm = norm(u);
% uNormalized = u/uNorm;
% 
% %check if the results are orthogonal basis
% checkOrtogonality(uNormalized);
% 
% end

[V,D] = eig(Y);
N = size(Y, 1);
clear Y
[L,I] = sort(diag(D)/N);

nL = length(L);
L = L(nL:-1:1);
I = I(nL:-1:1);

%Diagonal matrix containing the square roots of the eigenvalues:
S = sqrt(diag(D));
S = S(I);
V = V(:,I);

if numberOfModes==0
    numberOfModes = numel(extractK(D, 10e-6));
end

a = diag(S(1:numberOfModes))*(V(:,1:numberOfModes)');
uNormalized = T*V(:,1:numberOfModes)*diag(1./S(1:numberOfModes));

checkOrtogonality(uNormalized);

end

function [ k ] = extractK(S, toll)
    k = [];
    for i=1:size(S)
        if(abs(S(i,i)/max(max(S))) >= toll)
            k = [k i];
        end
    end
end

function checkOrtogonality(A)

OrthogonalIdentity = A*A';
tol=1.0e-08;

identityNorm = norm(OrthogonalIdentity);

if(identityNorm<1-tol) || (identityNorm>1+tol)
     'Warning: the POD is not orthogonal!!!' ;
end

end


