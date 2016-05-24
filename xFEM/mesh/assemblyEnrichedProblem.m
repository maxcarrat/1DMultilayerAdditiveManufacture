function [M, K, f] = assemblyEnrichedProblem(problem)
%   [M, K, f] = ASSEMBLYENRICHEDPROBLEM(problem) assembles the mass and the conductivity matrix and load vector 
%   problem = definition of the boundary value problem

    %global conductivity matrix of the coarse problem
    Kc = zeros(problem.cdof,problem.cdof);
    %global capacity matrix of the coarse problem
    Mc = zeros(problem.cdof,problem.cdof);
    %load vector of the coarse problem
    fc = zeros(problem.cdof, 1);

    KE = localConductivityMatrix(problem);
    ME = localCapacityMatrix(problem);
    
    %First assembly the standard DOFs of the coarse mesh...
    for e=1:problem.N
        
        ldof = 2;
        X1 = problem.coords(e);
        X2 = problem.coords(e+1);
        
        fel = localLoadVector(e, problem);

        fc(problem.LM(e,1:ldof)) = fc(problem.LM(e,1:ldof)) + problem.F_map(X1,X2) * fel;
        Mc(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = Mc(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof))...
            + problem.F_map(X1,X2) * ME(1:ldof, 1:ldof);
        Kc(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = Kc(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof))...
            + problem.B_map(X1,X2) * KE(1:ldof, 1:ldof);
   
    end

    %...secondly, assemble the enrichment DOFs and the coupling entries of
    %the matrices. The corresponding enrichment fucntions are supported by 
    %the last two elements of the mesh, thus the assembly has to be 
    %splitted into two parts:
    %(1) assemble the enrichment overlapping meshes
    
    %conductivity matrix of the enriched problem
    Ke = zeros(problem.edof*problem.modes,problem.edof*problem.modes);
    %capacity matrix of the enriched problem
    Me = zeros(problem.edof*problem.modes,problem.edof*problem.modes);
    %load vector of the enriched problem
    fe = zeros(problem.edof*problem.modes, 1);
    
    enrichedElementCoords = linspace(problem.N-2, problem.N, 2^problem.refinementDepth);
    
    for iMode=1:problem.modes
        KE = rbLocalConductivityMatrix(problem, iMode, enrichedElementCoords);
    	ME = rbLocalCapacityMatrix(problem, iMode, enrichedElementCoords);
    
        for modesSupports=1:2
            ldof = 2;
            X1 = problem.coords(problem.N-2);
            X2 = problem.coords(problem.N-1);
            
            fel = rbLocalLoadVector(problem, iMode, enrichedElementCoords, problem.N);
            
            rbLM = problem.rbLM;
            
            fe(problem.rbLM(modesSupports,1:ldof)) = fe(problem.rbLM(modesSupports,1:ldof)) + problem.F_map(X1,X2) * fel;
            Me(problem.rbLM(modesSupports, 1:ldof), problem.rbLM(modesSupports, 1:ldof)) = Me(problem.rbLM(modesSupports, 1:ldof), problem.rbLM(modesSupports, 1:ldof))...
                + problem.F_map(X1,X2) * ME(1:ldof, 1:ldof);
            Ke(problem.rbLM(modesSupports, 1:ldof), problem.rbLM(modesSupports, 1:ldof)) = Ke(problem.rbLM(modesSupports, 1:ldof), problem.rbLM(modesSupports, 1:ldof))...
                + problem.B_map(X1,X2) * KE(1:ldof, 1:ldof);
        end
    end
    
  
end