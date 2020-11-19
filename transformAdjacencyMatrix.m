%--------------------------------------------------------------%
% Author: Miel Sharf
% Date 06/08/20
% 
% Description: Transform the adjacency matrix created by CreateTransitionNetwork
%  into an representation description that is easier to compute with.
% 
% Returns:
% SV_states - maps SV values (columns) to a linear index (in rows).
% time_from - time step at of each 'from' transition
% dT - number of time steps of the transition
% from_state_map - maps number assigned by SV_states to the 'from' node
% from_state_map - maps number assigned by SV_states to the 'to' node
% time_map - maps time step (not used)
%--------------------------------------------------------------%

function [SV_states, time_from, n_tsteps, from_state_map, to_state_map] = ...
          transformAdjacencyMatrix(state_from, state_to, svToStateNumber)
% Move the s=0,v=0 to the end of the state list.
nStates = size(svToStateNumber,1);
SV_states = [svToStateNumber(2:end, :); svToStateNumber(1, :)];
%node = 1, node = last are special nodes - shift all others by 1 and use
%svToStateNumber to find out the order of the others.

SourceNode = 1;
TerminalNode = max(max(state_from), max(state_to));

from_state_map = zeros(size(state_from));
to_state_map = zeros(size(state_to));
FromSourceNodeIndices = (state_from == SourceNode);
ToTerminalNodeIndices = (state_to == TerminalNode);
FromRegularNodeIndices = ~FromSourceNodeIndices;
ToRegularNodeIndices = ~ToTerminalNodeIndices;

from_state_map(FromSourceNodeIndices) = nStates; %(0,0)
to_state_map(ToTerminalNodeIndices) = nStates; %(0,0)

%Shift state_from and state_to:
state_from(FromRegularNodeIndices) = state_from(FromRegularNodeIndices) - 1;
state_to(ToRegularNodeIndices) = state_to(ToRegularNodeIndices) - 1;

from_state_map(FromRegularNodeIndices) = mod(state_from(FromRegularNodeIndices)-1, nStates);
to_state_map(ToRegularNodeIndices) = mod(state_to(ToRegularNodeIndices)-1, nStates);
from_state_map(from_state_map == 0) = nStates;
to_state_map(to_state_map == 0) = nStates;

% Determine the starting and ending times of a transition.
time_from = zeros(size(state_from));
time_to = zeros(size(state_to));

time_from(FromSourceNodeIndices) = 0;
time_to(ToTerminalNodeIndices) = intmax('uint16');

time_from(FromRegularNodeIndices) = floor((state_from(FromRegularNodeIndices) - 0.5)/nStates) + 1;
time_to(ToRegularNodeIndices) = floor((state_to(ToRegularNodeIndices) - 0.5)/nStates) + 1;

n_tsteps = time_to - time_from; %dT
end

% SV_states - all possible (s,v) pairs. off is equal to (0,0)
% from_state_map - the index (in SV_states) of the from state of the edge.
% to_state_map - the index (in SV_states) of the to state of the edge.
% time_from - time of from state of the edge.
% dT - Transition length