function [problemStruct, err] = getoptimrandstates(problemStruct,hashProb)
%GETOPTIMRANDSTATES get states of random number generator for solvers
%
%   Private to OPTIMTOOL

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2007/12/10 21:50:30 $

err = '';
% Update or create the appdata structure before getting it
setoptimrandstates(problemStruct,hashProb,false);
if isappdata(0,'optim_rand_states_struct')
    optimSolversRandStates =  getappdata(0,'optim_rand_states_struct');
else % random states are not available; return
    return;
end

problemStruct.randstate = [];
problemStruct.randnstate = [];

switch problemStruct.solver
    case 'ga'
         if optimSolversRandStates.ga.garandchoice
            problemStruct.randstate = optimSolversRandStates.ga.randstate;
            problemStruct.randnstate = optimSolversRandStates.ga.randnstate;
         end
    case 'gamultiobj'
         if optimSolversRandStates.gamultiobj.gamultiobjrandchoice
            problemStruct.randstate = optimSolversRandStates.gamultiobj.randstate;
            problemStruct.randnstate = optimSolversRandStates.gamultiobj.randnstate;
        end
    case 'patternsearch'
         if optimSolversRandStates.patternsearch.patternsearchrandchoice
            problemStruct.randstate = optimSolversRandStates.patternsearch.randstate;
            problemStruct.randnstate = optimSolversRandStates.patternsearch.randnstate;
        end
    case 'threshacceptbnd'
         if optimSolversRandStates.threshacceptbnd.threshacceptbndrandchoice
            problemStruct.randstate = optimSolversRandStates.threshacceptbnd.randstate;
            problemStruct.randnstate = optimSolversRandStates.threshacceptbnd.randnstate;
        end
    case 'simulannealbnd'
         if optimSolversRandStates.simulannealbnd.simulannealbndrandchoice
            problemStruct.randstate = optimSolversRandStates.simulannealbnd.randstate;
            problemStruct.randnstate = optimSolversRandStates.simulannealbnd.randnstate;
        end
end
