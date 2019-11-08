function hashModel = setoptimrandstates(problemStruct,hashModel,randchoice)
%SETOPTIMRANDSTATES set states of random number generator for solvers
%
%   Private to OPTIMTOOL

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2007/12/10 21:50:46 $

% Get structure that stores the states
if isappdata(0,'optim_rand_states_struct')
    optimSolversRandStates = getappdata(0,'optim_rand_states_struct');
else % Create a structure that can store all the random states
    optimSolversRandStates.ga.randstate = [];
    optimSolversRandStates.ga.randnstate = [];
    optimSolversRandStates.ga.garandchoice = false;

    optimSolversRandStates.gamultiobj.randstate = [];
    optimSolversRandStates.gamultiobj.randnstate = [];
    optimSolversRandStates.gamultiobj.gamultiobjrandchoice = false;

    optimSolversRandStates.patternsearch.randstate = [];
    optimSolversRandStates.patternsearch.randnstate = [];
    optimSolversRandStates.patternsearch.patternsearchrandchoice = false;

    optimSolversRandStates.threshacceptbnd.randstate = [];
    optimSolversRandStates.threshacceptbnd.randnstate = [];
    optimSolversRandStates.threshacceptbnd.threshacceptbndrandchoice = false;

    optimSolversRandStates.simulannealbnd.randstate = [];
    optimSolversRandStates.simulannealbnd.randnstate = [];
    optimSolversRandStates.simulannealbnd.simulannealbndrandchoice = false;

    % Save in the appdata
    setappdata(0,'optim_rand_states_struct',optimSolversRandStates);
end

% Return if the hashModel does not have 'randstate' and 'randnstate' keys
% Any modification in the hashModel is okay when nargout == 0; In this
% case, simply set the states in the appdata for the correct solver
if nargout > 0 && ~hashModel.containsKey('randstate') && ~hashModel.containsKey('randnstate')
    return;
else
    % Remove randstate and randnstate from the hashModel. The GUI only
    % needs to know if randstates are present or not (a boolean flag
    % 'randchoice') 
    hashModel.remove('randstate');
    hashModel.remove('randnstate');
    if randchoice
        useCurrentRandState = true;
    else
        useCurrentRandState = false;
    end
end

% Get the random states from the problem structure and save in the appdata
randstate = problemStruct.randstate;
randnstate = problemStruct.randnstate;

if ~isempty(randstate) && ~isempty(randnstate)
    switch problemStruct.solver
        case 'ga'
            optimSolversRandStates.ga.randstate = randstate;
            optimSolversRandStates.ga.randnstate = randnstate;
        case 'gamultiobj'
            optimSolversRandStates.gamultiobj.randstate = randstate;
            optimSolversRandStates.gamultiobj.randnstate = randnstate;
        case 'patternsearch'
            optimSolversRandStates.patternsearch.randstate = randstate;
            optimSolversRandStates.patternsearch.randnstate = randnstate;
        case 'threshacceptbnd'
            optimSolversRandStates.threshacceptbnd.randstate = randstate;
            optimSolversRandStates.threshacceptbnd.randnstate = randnstate;
        case 'simulannealbnd'
            optimSolversRandStates.simulannealbnd.randstate = randstate;
            optimSolversRandStates.simulannealbnd.randnstate = randnstate;
    end
end

switch problemStruct.solver
    case 'ga'
        if hashModel.containsKey('garandchoice') 
            if strcmpi('true',hashModel.get('garandchoice'))
                optimSolversRandStates.ga.garandchoice = true;
            elseif strcmpi('false',hashModel.get('garandchoice'))
                optimSolversRandStates.ga.garandchoice = false;
            end
        elseif useCurrentRandState
            optimSolversRandStates.ga.garandchoice = true;
            hashModel.put('garandchoice','true');
        end
    case 'gamultiobj'
        if hashModel.containsKey('gamultiobjrandchoice')
            if strcmpi('true',hashModel.get('gamultiobjrandchoice'))
                optimSolversRandStates.gamultiobj.gamultiobjrandchoice = true;
            elseif strcmpi('false',hashModel.get('gamultiobjrandchoice'))
                optimSolversRandStates.gamultiobj.gamultiobjrandchoice = false;
            end
        elseif useCurrentRandState
            optimSolversRandStates.gamultiobj.gamultiobjrandchoice = true;
            hashModel.put('gamultiobjrandchoice','true');
        end
    case 'patternsearch'
        if hashModel.containsKey('patternsearchrandchoice')
            if strcmpi('true',hashModel.get('patternsearchrandchoice'))
                optimSolversRandStates.patternsearch.patternsearchrandchoice = true;
            elseif strcmpi('false',hashModel.get('patternsearchrandchoice'))
                optimSolversRandStates.patternsearch.patternsearchrandchoice = false;
            end
        elseif useCurrentRandState
            optimSolversRandStates.patternsearch.patternsearchrandchoice = true;
            hashModel.put('patternsearchrandchoice','true');
        end
    case 'threshacceptbnd'
        if hashModel.containsKey('threshacceptbndrandchoice') 
            if strcmpi('true',hashModel.get('threshacceptbndrandchoice'))
                optimSolversRandStates.threshacceptbnd.threshacceptbndrandchoice = true;
            elseif strcmpi('false',hashModel.get('threshacceptbndrandchoice'))
                optimSolversRandStates.threshacceptbnd.threshacceptbndrandchoice = false;
            end
        elseif useCurrentRandState
            optimSolversRandStates.threshacceptbnd.threshacceptbndrandchoice = true;
            hashModel.put('threshacceptbndrandchoice','true');
        end
    case 'simulannealbnd'
        if hashModel.containsKey('simulannealbndrandchoice') 
            if strcmpi('true',hashModel.get('simulannealbndrandchoice'))
                optimSolversRandStates.simulannealbnd.simulannealbndrandchoice = true;
            elseif strcmpi('false',hashModel.get('simulannealbndrandchoice'))
                optimSolversRandStates.simulannealbnd.simulannealbndrandchoice = false;
            end
        elseif useCurrentRandState
            optimSolversRandStates.simulannealbnd.simulannealbndrandchoice = true;
            hashModel.put('simulannealbndrandchoice','true');
        end
end
% Save in the appdata
setappdata(0,'optim_rand_states_struct',optimSolversRandStates);
