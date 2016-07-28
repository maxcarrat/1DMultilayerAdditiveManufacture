function [ N, B ] = BsplinesShapeFunctionsAndDerivatives( x, p, knotVector )
%BSPLINESSHAPEFUNCTIONSANDDERIVATIVES Generate the Bsplines shape function
%vector N and its derivative B

% B-SPline basis functions
N = BSplineBasis(x, p, knotVector);

% B-SPline basis functions and their derivative
B = BSplineBasisDerivative(x, p, knotVector);

end


function N = BSplineBasis(x, pmax, knotVector)
% Compute the nonvanishing basis function
% input:
% x = locl coordinate
% pmax = max polynomial order
% knotVector = knot vector
% output:
% N = BSpline basis function

m=length(knotVector);

N(m-1,pmax+1)=0;

% Order p=0
for i=1:m-1
    if x>=knotVector(i) && x<=knotVector(i+1)
        N(i,1)=1;
    end
    
end

%Order p>0
for p=2:pmax+1
    
    for i=1:m-p
       
        if abs(knotVector(i+p-1)-knotVector(i))>2*eps
            a=(x-knotVector(i))/(knotVector(i+p-1)-knotVector(i));
        else
            a=0;
        end
        
        if abs(knotVector(i+p)-knotVector(i+1))>2*eps
        b=(knotVector(i+p)-x)/(knotVector(i+p)-knotVector(i+1));      
        else
        b=0;
        end
        
        N(i,p)=a*N(i,p-1)+b*N(i+1,p-1);
        N(i,p)=min(N(i,p),1);
    
    end
end

N=N(1:m-pmax-1,pmax+1)'; %select only shape functions of the order pmax
end


function dN=BSplineBasisDerivative(x, p, knotVector)
% Compute the non-zero basis function derivative
% input:
% x = locl coordinate
% p = polynomial order
% knotVector = knot vector
% output:
% dN = BSpline basis function derivatives

m=length(knotVector);

dN(m-p-1)=0;

if p > 0
    NOneOrderLower=BSplineBasis(x, p-1, knotVector);
    
    for i=1:m-p-1
        
        if abs(knotVector(i+p)-knotVector(i))>2*eps
            a=p/(knotVector(i+p)-knotVector(i));
        else
            a=0;
        end
        
        if abs(knotVector(i+p+1)-knotVector(i+1))>2*eps
            b=p/(knotVector(i+p+1)-knotVector(i+1));
        else
            b=0;
        end
        
        dN(i)=a*NOneOrderLower(i)-b*NOneOrderLower(i+1);
        
    end
    
end

end

