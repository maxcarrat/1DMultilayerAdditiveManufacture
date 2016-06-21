%% Generate Convergence Study
% Read data from .txt file and plot the convergnce study of h-FEM and
% XPODFEM w.r.t. Reference overkilled solution.

figure(10)
% Create axes
axes1 = axes;

referenceFile = fopen('myReferenceResultsFile.txt','r');
formatSpec = '%f;\n';
referenceTemperatureSolution = fscanf(referenceFile, formatSpec);

% %% hFEM plot
% % plot the rlative error on the temperature L2 norm using diffrent level of
% % refinement for classical FEM discretization
% 
% dofsVector=[41, 81, 161, 321, 641, 1281, 2561, 5121]; % refDepth = 1,..,8
% relError = [];
% for depth=1:4
%     formatSpec = 'myFEMResultsFile_%d.txt';
%     FEMFileID = sprintf(formatSpec,depth);
%     FEMFile = fopen(FEMFileID, 'r');
%     FEMTemperatureSolution = fscanf(FEMFile, '%f;\n');
%     
%     err = 0.0;
%     for i=1:500%numel(FEMTemperatureSolution)
%         err = err + sqrt(( referenceTemperatureSolution(end - 500 + i) - FEMTemperatureSolution(end - 500 + i))^2 ...
%             / referenceTemperatureSolution(end - 500 + i)^2 );
%     end
%     
%     relError = [relError err/numel(FEMTemperatureSolution)];
%     
% end
% % Create multiple lines using matrix input to plot
% loglog(dofsVector(1:4),relError,'-kd');
% hold on

%% XPODFEM plot
% plot the rlative error on the temperature L2 norm using diffrent level of
% refinement for XPODFEM discretization

dofsVectorPOD=[22, 23, 24, 25, 26, 27, 28, 29]; % PODmodes = 1,..,8
relErrorPOD = [];

referenceFile = fopen('myFEMResultsFile_6.txt','r');
formatSpec = '%f;\n';
FEMTemperatureSolution = fscanf(referenceFile, formatSpec);

for modes=1:4
    formatSpec = 'myPODXFEMResultsFile_%d.txt';
    XFEMFileID = sprintf(formatSpec,modes);
    XFEMFile = fopen(XFEMFileID, 'r');
    XFEMTemperatureSolution = fscanf(XFEMFile, '%f;\n');
    
    err = 0.0;
    for i=1:numel(XFEMTemperatureSolution)
        err = err + sqrt(( FEMTemperatureSolution(i) - XFEMTemperatureSolution(i))^2 ...
            / FEMTemperatureSolution(i)^2 );
    end
    
    relErrorPOD = [relErrorPOD err/numel(XFEMTemperatureSolution)];
    
end

% Create multiple lines using matrix input to plot 
loglog(dofsVectorPOD(1:4),relErrorPOD,'-.r*');
hold on

% Title
title({'\centering{\quad \quad Convergence plot XPODFEM}'; '\centering{using refined hFEM (depth 6) as reference}'}...
    , 'Interpreter','latex');


% Axis labels
ylabel('\fontname{Latin Modern Math} Relative error on temperature');
xlabel('\fontname{Latin Modern Math} DOFs');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');


hold off
