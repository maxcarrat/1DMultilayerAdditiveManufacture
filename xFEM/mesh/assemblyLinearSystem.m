function [M, K, f] = assemblyLinearSystem(problem, time, integrationOrder)
%   [M, K, f] = ASSEMBLY(problem) assembles the mass and the conductivity matrix and load vector 
%   problem = definition of the boundary value problem

    %global conductivity matrix
    K = zeros(problem.gdof,problem.gdof);
    %global capacity matrix
    M = zeros(problem.gdof,problem.gdof);
    %load vector
    f = zeros(problem.gdof, 1);
    %gauss points
    [rGP, wGP] = gaussPoints( integrationOrder );
    
    numberOfIntegrationPoints = length(rGP);
    
    
    for e=1:problem.N
        
        ldof = 2;
        X1 = problem.coords(e);
        X2 = problem.coords(e+1);
        
        for iGP = 1:numberOfIntegrationPoints
            
            [N, B] = shapeFunctionsAndDerivatives(rGP(iGP));
            
            %extrnal heat source
            f(problem.LM(e,1:ldof)) = f(problem.LM(e,1:ldof)) + problem.F_map(X1,X2) * N' * problem.rhs(mapLocalToGlobal(rGP(iGP), X1, X2), time) * wGP(iGP);
            
            %Capacity matrix
            M(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = M(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
                                                    problem.F_map(X1,X2) * problem.heatCapacity * N' * N * wGP(iGP);
            
            %Diffusion matrix
            K(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = K(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
                                                    problem.B_map(X1,X2) * problem.k * B' * B * wGP(iGP);
            
        end
    end
end