function [Adj] = buildStateAdjacencyMatrix(Dictionary_uint16, T_Startup, T_Shutdown)
%%%%%% buildStateAdjacencyMatrix.
%{
Input:
Dictionary: a States x 3 array. The first column defines the s-variable
of the state (engine speed), the second column defines the v-variable
of the state (bypass valve), and the third number describes the key,
which is a uint32. The (0,0) state corresponds to 'off'.
If Dictionary only has 2 columns, add a third one as 1:States.

T_Startup: Amount of time [???] it takes to complete a startup.

T_Shutdown: Amount of time [???] it takes to complete a shutdown.

Output:
Adj: A States x States x T array. A(:,:,t) is the adjacency matrix of
all transitions which take exact t units of time. Namely, A(i,j,t) == 1
if and only if the transition from the state with key i to the state
with key j is valid and takes t units of time.

Main idea - first build an adjacency matrix with lexicographic order,
then use the dictionary to perform a "permutation" giving the desired
adjacency matrix.

Lexicographic order: (s_1,v_1) is before (s_2,v_2) if s_1<s_2 or if
s_1=s_2 and v_1 < v_2. In that case, if there are S total speeds
(1,2,...,S) and V total bypass valve positions (0,1,...,V-1), then the
state (s,v) is in the k-th position where k = V*(s-1)+v+1.
%}

Dictionary = double(Dictionary_uint16);
%If Dictionary only has 2 columns - Add a key column 1:States.
if (size(Dictionary, 2) == 2)
  Dictionary(:, 3) = (1:size(Dictionary, 1))';
end

S = max(Dictionary(:, 1));
V = max(Dictionary(:, 2)) + 1;
%Check that all numbers between 1 and S*V+1 are used as keys. There are a
%total of S*V on states and one off state.
assert(sum(sort(Dictionary(:, 3)) ~= (1:S * V + 1).') == 0);

%% State all possible transition laws. This part can be changed if a new turbine model is used.
% The on state i corresponds to (s,v) where s = floor((OnIndexArray-1)./V)+1;
% and v = mod(OnIndexArray-1,V)+1. Conversely, i = V*(s-1)+v+1.

%Transition laws between two on states.
StayTheSame = @(i, j) (j == i); %(s,v) -> (s,v)
IncreaseSByOne = @(i, j) (j == i + V); %(s,v) -> (s+1,v)
DecreaseSByOne = @(i, j) (j == i - V); %(s,v) -> (s-1,v)
%Tricks with V must be done carefully to assure that (1,V-1) is not
%conneted to (2,0).
IncreaseVByOne = @(i, j) ((j == i + 1) & (mod(i-1, V) <= V - 2)); %(s,v) -> (s,v+1)
DecreaseVByOne = @(i, j) ((j == i - 1) & (mod(i-1, V) >= 1)); %(s,v) -> (s,v-1)

IncreaseSByOneIncreaseVByOne = @(i, j) ((j == i + V + 1) & (mod(i-1, V) <= V - 2)); %(s,v) -> (s+1,v+1)
IncreaseSByOneIncreaseVByTwo = @(i, j) ((j == i + V + 2) & (mod(i-1, V) <= V - 3)); %(s,v) -> (s+1,v+2)
IncreaseSByOneDecreaseVByOne = @(i, j) ((j == i + V - 1) & (mod(i-1, V) >= 1)); %(s,v) -> (s+1,v-1)
IncreaseSByOneDecreaseVByTwo = @(i, j) ((j == i + V - 2) & (mod(i-1, V) >= 2)); %(s,v) -> (s+1,v-2)
DecreaseSByOneIncreaseVByOne = @(i, j) ((j == i - V + 1) & (mod(i-1, V) <= V - 2)); %(s,v) -> (s-1,v+1)
DecreaseSByOneDecreaseVByOne = @(i, j) ((j == i - V - 1) & (mod(i-1, V) >= 1)); %(s,v) -> (s-1,v-1)

%Transition laws between an on state and the off state.
StartUp = @(j) (j >= V*(S - 1)+1); %Startup to max engine speed.
ShutDown = @(j) (j == 1); %Shutdown only from (s=1,v=0) state.

num_states = S * V + 1;
Adj = zeros(num_states, num_states, max([T_Startup, T_Shutdown, 2]), 'uint16');

% Transitions of Length 1
Adj(:, :, 1) = Adj(:, :, 1) + LexOrderLawAdjMatrix_NoOffState(S, V, StayTheSame);
Adj(:, :, 1) = Adj(:, :, 1) + LexOrderLawAdjMatrix_NoOffState(S, V, DecreaseSByOne);
Adj(:, :, 1) = Adj(:, :, 1) + LexOrderLawAdjMatrix_NoOffState(S, V, IncreaseVByOne);
Adj(:, :, 1) = Adj(:, :, 1) + LexOrderLawAdjMatrix_NoOffState(S, V, DecreaseVByOne);
Adj(:, :, 1) = Adj(:, :, 1) + LexOrderLawAdjMatrix_NoOffState(S, V, DecreaseSByOneIncreaseVByOne);
Adj(:, :, 1) = Adj(:, :, 1) + LexOrderLawAdjMatrix_NoOffState(S, V, DecreaseSByOneDecreaseVByOne);
% Add a valid transition from the off state to itself of length 1.
Adj(end, end, 1) = 1;


%Transitions of Length 2
Adj(:, :, 2) = Adj(:, :, 2) + LexOrderLawAdjMatrix_NoOffState(S, V, IncreaseSByOne);
Adj(:, :, 2) = Adj(:, :, 2) + LexOrderLawAdjMatrix_NoOffState(S, V, IncreaseSByOneIncreaseVByOne);
Adj(:, :, 2) = Adj(:, :, 2) + LexOrderLawAdjMatrix_NoOffState(S, V, IncreaseSByOneIncreaseVByTwo);
Adj(:, :, 2) = Adj(:, :, 2) + LexOrderLawAdjMatrix_NoOffState(S, V, IncreaseSByOneDecreaseVByOne);
Adj(:, :, 2) = Adj(:, :, 2) + LexOrderLawAdjMatrix_NoOffState(S, V, IncreaseSByOneDecreaseVByTwo);

%Startup Transitions
Adj(:, :, T_Startup) = Adj(:, :, T_Startup) + LexOrderLawAdjMatrix_ToFromOffState(S, V, StartUp, 0);

%Shutdown Transitions
Adj(:, :, T_Shutdown) = Adj(:, :, T_Shutdown) + LexOrderLawAdjMatrix_ToFromOffState(S, V, ShutDown, 1);

%% Permute the adjacency matrix to correspond to the dictionary.
% The permuted matrix satifies the following equation:
% Adj_Permut(key_i,key_j,:) = Adj(lex_i,lex_j,:)
% where the keys key_i,key_j correspond to elements x_i,x_j, which have
% lexicographic positions lex_i,lex_j.

% Use the formula V*(s-1)+v+1 for the lexicographic position of (s,v)
LexicographicOrder = V * (Dictionary(:, 1) - 1) + Dictionary(:, 2) + 1;
% This formula does not work for the off state, for which it gives a
% negative number. We recall that it was put in the last position.
LexicographicOrder(LexicographicOrder < 0) = S * V + 1;

%Permute.
Adj(Dictionary(:, 3), Dictionary(:, 3), :) = Adj(LexicographicOrder, LexicographicOrder, :);
end

%%% Create a partial adjacency matrix which corresponds to some law.
%%% Namely, A(i,j) = 1 if and only if Law(i,j) == 1, where
%%% i = V*(s1-1)+v1+1 and j = V*(s2-1)+v2+1. This law is applied
%%% only to `on' states and not to the `off' states.
function AdjLaw = LexOrderLawAdjMatrix_NoOffState(S, V, Law)
%Enumerates all possible states.
num_states = double(S*V);
%IndexPairArray stores inside the index of all pairs (i,j) where
%i,j=1,2,...,S*V are on states.
IndexPairArray = (1:num_states^2)';
AdjPartial = zeros(num_states, 'uint16');
% Matlab matrices can be called as arrays, which gives them columnwise
% (i.e., if A is a m x m matrix then A(2) is the 2,1 entry). Thus
% AdjPartial(k) corresponds to the pair of states (i,j) such that k =
% num_states*(j-1)+i (as k,i,j start from 1). Thus, we get A(k) = Law(i,j).
AdjPartial(:) = uint16(Law(mod(IndexPairArray-1, num_states)+1, floor((IndexPairArray - 1)./num_states)+1));

% Add the zero-state in the adjacency matrix... - it is added in the end.
AdjLaw = [AdjPartial, zeros(num_states, 1); zeros(1, num_states+1)];
end

%%% Create a Partial Adjacency Matrix corresponding to some law. Only
%%% consider transitions regarding the off state. If IsTo == 1, we consider
%%% transitions into the off stae. If IsTo == 0, we consider transitiosn
%%% from the off state.
function AdjLaw = LexOrderLawAdjMatrix_ToFromOffState(S, V, Law, IsTo)
%Enumerates all possible states.
num_states = S * V;
%IndexPairArray stores inside the index of all pairs (i,j) where
%i,j=1,2,...,S*V are on states.
IndexPairArray = (1:num_states)';
% For each on state, check if the transitions between it and the off state
% is valid or not.
AdjPartial = uint16(Law(IndexPairArray));
%Add zeros to makes this a proper adjacency matrix.
%In both cases, we put a zero in the last position at the final row even
%though we can %move from the off state to itself. We do this as this
%function will be called to specify the transitios which take T_shutdown
%or T_startup time units.
if (IsTo) %To
  AdjLaw = uint16([zeros(num_states), AdjPartial; zeros(1, num_states+1)]);
  return;
end
%Else - From
AdjLaw = uint16([zeros(num_states, num_states+1); AdjPartial.', 0]);
end
