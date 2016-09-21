%% XIGAIntegrationTest
% Test the XIga integration using 2 linear enrichment modes and quadratic
% C1 continuous IGA basis functions.

 % Generate a coarse IGA mesh
 p = 1;
 numberOfKnotSpan = 4;
 knotVector = getOpenKnotVector( numberOfKnotSpan, p );
 knotSpanDImension = 2.5;
 CPs = getControlPoints( numberOfKnotSpan, knotSpanDImension, p );
 
 
 % Generate the XIGA Poisson problem
 numberOfEnrichedCPs = 1;
 length = 10;
 duration = 1.0;
 numberOfTimeStep = 10;
 refinementDepth = 2;
 
 rhs = @(x, t) externalSource(x, t, length, numberOfKnotSpan, duration,...
     numberOfTimeStep, refinementDepth, 0.0);
 
 dirichletLeftBC = @(t) 0.0;
 dirichletRightBC = @(t) 0.0;
 
 k = @(x, t, T) 1.0;
 
 % heat capacity function
 heatCapacity= @(x, T, T_last)  1.0;
 
 currentTime = 1.0;
 
 solutionReductionOperator = [1.0, 1.0; 1.25, 0.5; 1.5, 0.0; 1.75, -0.5; 2.0, -1.0];
 
 poissonProblem = poissonProblemTransientXIGA(CPs, numberOfEnrichedCPs, rhs,...
     dirichletLeftBC, dirichletRightBC,...
     0.0, k, heatCapacity, currentTime,...
     knotVector, p, refinementDepth, numberOfEnrichedCPs,...
     solutionReductionOperator);
 
 timeStepSize = 1.0;
 integrationOrder = 3.0;
 integrationModesOrder = 3.0;
 tolerance = 1.0e-04;
 maxIterations = 100;
 temperatureSolutions = linspace(0.0, 0.0, 7);
 
 [temperatureSolutions, convergenceFlag] = solveXIGAMultiPhase( poissonProblem, currentTime,...
            timeStepSize, integrationOrder, integrationModesOrder, tolerance,...
            maxIterations, temperatureSolutions' );
        
 temperatureSolutions
        