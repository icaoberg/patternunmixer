function results = callSolver(probStruct,optStruct)
%callSolver Call appropriate solver and return the results structure.
%   The input arguments 'probStruct' and 'optStruct' are problem and options structure.
%   The output 'results' is a structure containing all the outputs from solver.
%
%   Private to OPTIMTOOL

%   Copyright 2005-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2007/11/09 20:56:41 $

solver = probStruct.solver;
% Update the problem structure 'probStruct' for the solver
probStruct = createProblemStruct(solver,[],probStruct);
% Add options to the problem structure (which is a required field)
probStruct.options = createOptionsStruct(solver,optStruct); 
% Create result structure for the solver
results = createResultsStruct(solver);
resultsFields = fieldnames(results);
numfields = length(resultsFields);
% Initialize cell array to capture the output from the solver
solverOutput = cell(1,numfields);
% Set 'optimization running' message in the GUI
startMessage = sprintf('%s\n','Optimization running.');
optimGUI = com.mathworks.toolbox.optim.OptimGUI.getOptimGUI; % Get a handle to the GUI
if ~isempty(optimGUI)
    optimGUI.appendResults(startMessage);
end
% Call solver and put results in the structure
[solverOutput{:}] = feval(str2func(solver),probStruct);
for i = 1:numfields
    results.(resultsFields{i}) = solverOutput{i};
end
% Set 'optimization terminated' message in the GUI
endMessage = sprintf('%s\n','Optimization terminated.');
optimGUI = com.mathworks.toolbox.optim.OptimGUI.getOptimGUI; % Get a handle to the GUI
if ~isempty(optimGUI)
    optimGUI.appendResults(endMessage);
end

