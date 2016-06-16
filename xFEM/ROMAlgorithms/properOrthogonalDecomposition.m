function [ uNormalized, numberOfModes ] = properOrthogonalDecomposition( T, numberOfModes )
%% Proper Orthogonal decomposition
% generate a reduced operator, from J.G.Michopulous
% "Performance of Reduced Order Models of Moving Heat Source Deposition
% Problems for efficient inverse analysis", 2014
%T = matrix of solution snapshots

N = size(T,2);
Y = T(:,:)'*T(:,:)*1/N;

% eigendecomposition, return D as matrix of the
% eigenvectors and S as diagonal matrix of the eigenvalues.

% [D, S, ~] = svd(Y);
% [D, S] = eig(Y);
% 
% [ S, D, numberOfModesFiltered ] = sortEigenvalues(S, D);
% 
% % if numberOfModes==0
% %     k = extractK(S, 10e-6);                          % number of highly energetic eigenvalues
% %     u = zeros(size(T,1), size(k,2));                  % basis vectors matrix
% %     numberOfModes = size(k,2);
% % end
% 
% if numberOfModes==0
%     numberOfModes = numberOfModesFiltered;
% end
% 
% u = zeros(size(T,1), numberOfModes);
% uNormalized = zeros(size(T,1), numberOfModes);
% for i=1:numberOfModes
%     for k=1:N
%         u(:,i)=  D(k,i)*T(:,k) + u(:,i);
%     end
%     uNormalized(:,i) = 1/(N*S(i))*u(:,i);
% end
% 
% 
% 
% % sum = zeros(size(T,1), 1);
% % for i=1:numberOfModes
% %     for j=1:numberOfModes
% % %          sum = D(j,i)*T(:,j);
% %         sum = sum + D(j,i)*T(:,j);
% %     end
% % %     u(:,i) = 1/sqrt(S(i,i)) * sum;
% % u(:,i) = sum;
% % end
% % 
% % % normalization of the basis vector
% % % uNorm = norm(u);
% % uNorm = 0.0;
% % for j=1:numberOfModes
% %     for i=1:size(u,1)
% %         uNorm = uNorm + u(i,j)^2;
% %     end
% % end
% % uNorm = sqrt(uNorm);
% % 
% % for i=1:numberOfModes
% %     uNormalized(:,i) = u(:,i)/uNorm;
% % end
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
    numberOfModes = numel(extractK(D, 10e-8));
end

a = diag(S(1:numberOfModes))*(V(:,1:numberOfModes)');
uNormalized = T*V(:,1:numberOfModes)*diag(1./S(1:numberOfModes));

checkOrtogonality(uNormalized);

end


function [ S, D, numberOfModes ] = sortEigenvalues(S, D)
    k = extractK(S, 1.0e-08);
    numberOfModes = numel(k);
    
    S = diag(S);
    
    % Bubble sort backward
    n = numel(S);
    
    while (n > 0)
        % Iterate through S
        nnew = 0;
        for i = 2:n
            % Swap eigenvalues and eigenvectors in wrong order
            if (S(i) > S(i - 1))
                S = swapEigenvalues(S,i,i - 1);
                D = swapEigenvectors(D,i,i-1);
                nnew = i;
            end
        end
        n = nnew;
    end
    
end


function S = swapEigenvalues(S,i,j)
% Swap Eigenvalues

val = S(i);
S(i) = S(j);
S(j) = val;

end

function D = swapEigenvectors(D,i,j)
% Swap Eigenvalues

val = D(:,i);
D(:,i) = D(:,j);
D(:,j) = val;

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


