function [ uNormalized, numberOfModes ] = properOrthogonalDecomposition( T )
%% Proper Orthogonal decomposition
% generate a reduced operator, from J.G.Michopulous
% "Performance of Reduced Order Models of Moving Heat Source Deposition
% Problems for efficient inverse analysis", 2014
%T = matrix of solution snapshots

Y = T(:,:)'*T(:,:);

% % eigendecomposition, return D as matrix of the
% % eigenvectors and S as diagonal matrix of the eigenvalues.
% 
% [D, S, ~] = svd(Y);
% 
% k = extractK(S, 10e-2);                          % number of highly energetic eigenvalues
% u = zeros(size(T,1), size(k,2));                  % basis vectors matrix
% 
% for i=1:size(k,2)
%     for j=1:size(k,2)
%          sum = D(j,i)*T(:,j);
%     end
%     u(:,i) = 1/sqrt(S(i,i)) * sum;
% end
% 
% numberOfModes = size(k,2); 
% 
% % normalization of the basis vector
% uNorm = norm(u);
% uNormalized = u/uNorm;
% 
% end
% 
% 
% function [ k ] = extractK(S, toll)
%     k = [];
%     for i=1:size(S)
%         if(abs(S(i,i)/max(max(S))) >= toll)
%             k = [k i];
%         end
%     end
% end

[V,D] = eig(Y);
N = size(Y, 1);
clear Y
[L,I] = sort(diag(D)/N);

nL = length(L);
L = L(nL:-1:1);
I = I(nL:-1:1);

% Diagonal matrix containing the square roots of the eigenvalues:
S = sqrt(diag(D));
S = S(I);
V = V(:,I);
m = max(extractK(S, 10e-2));
a = diag(S(1:m))*(V(:,1:m)');
uNormalized = T*V(:,1:m)*diag(1./S(1:m));

numberOfModes = m;
end

function [ k ] = extractK(S, toll)
    k = [];
    for i=1:size(S)
        if(abs(S(i)/max(max(S))) >= toll)
            k = [k i];
        end
    end
end


