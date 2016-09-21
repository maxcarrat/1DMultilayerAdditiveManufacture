function [ K, f ] = applyWeakDirichletBCs(problem, constrainedNodeGlobal, dirichletBCs, K, f)
%APPLYWEAKDIRICHLETBCS apply Weak Dirichlet BCs at a given node
%   problem = transient Poisson problem struct
%   constrainedNode = node of the mesh to be constrained
%   element =
%   K = conductivity matrix
%   f = temperature source vector


%% Constrain conductivity matrix K
ldof = 2;
numberOfModes = problem.modes;

enrichedElementCoords = linspace(problem.coords(end-1), problem.coords(end), 2^problem.refinementDepth+1);
constrainedNode = mapGlobalToLocal(constrainedNodeGlobal, problem.coords(end-1), problem.coords(end));

%sub-matrices
M_FEM = zeros(ldof, ldof);
M_enr = zeros(numberOfModes*ldof, numberOfModes*ldof);
M_coupling = zeros(numberOfModes*ldof, ldof);

for i=1:ldof
    for j=1:ldof
        M_FEM(i,j) = problem.basis_fun(constrainedNode, i, 0.0)*...
           problem.basis_fun(constrainedNode, j, 0.0);
    end
end

for iMode = 1:numberOfModes
    for i=1:ldof
        for jMode = 1:numberOfModes
            for j=1:ldof
                M_enr((iMode-1)*ldof + i,(jMode-1)*ldof + j) = (problem.xFEMBasis_fun(constrainedNode, i, iMode, 0.0, 0.0, problem, enrichedElementCoords )*...
                    problem.xFEMBasis_fun(constrainedNode, j, jMode, 0.0, 0.0, problem, enrichedElementCoords));
            end
        end
    end
end

for j=1:ldof
    for i=1:ldof
        for iMode = 1:numberOfModes
            M_coupling((iMode-1)*ldof + i,j) = problem.xFEMBasis_fun(constrainedNode, i, iMode, 0.0, 0.0, problem, enrichedElementCoords )*...
                problem.basis_fun(constrainedNode, j, 0.0);
        end
    end
end

Mpenalty = [M_FEM, M_coupling'; M_coupling, M_enr];

K = K + problem.penalty*Mpenalty;

%% Constrain source vector f
ldof = 2;
numberOfModes = problem.modes;

%sub-matrices
F_FEM = zeros(ldof, 1);
F_enr = zeros(numberOfModes*ldof, 1);

for i=1:ldof
    F_FEM(i,1) = problem.basis_fun(constrainedNode, i, 0.0)*dirichletBCs;
end

for iMode = 1:numberOfModes
    for i=1:ldof
        F_enr((iMode-1)*ldof + i, 1) =...
            problem.xFEMBasis_fun(constrainedNode, i, iMode, 0.0, 0.0, problem, enrichedElementCoords)*dirichletBCs;
    end
end

Fpenalty = [F_FEM; F_enr];

f = f + problem.penalty*Fpenalty;
end

