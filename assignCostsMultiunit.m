%--------------------------------------------------------------%
% File: assignCostsMultiunit.m (script)
% Author: Miel Sharf
% Date 08/8/2020
% v2.0
% Description: Assign costs to edges for the shortest path solver/
% Included transition penalty to deal with oscillations.
% Costs only model fuel costs for the turbines - the demand is dealt with
% later.
%--------------------------------------------------------------%
function [decided_costs] = assignCostsMultiunit(dt, fuel_map, mdot_fuel_SU, total_nodes, sol_select, n_tsteps, from_state_map, to_state_map, price_kg_f, is_transition)

%% Definitions:
Startup_cost = 3.75; % SU and SD consts are the same
Transition_Penalty = 0; %0.01; % Higher costs for transitions.
Weight1 = 1;
Weight2 = 1;

warning off 'MATLAB:strrep:InvalidInputType'

% average mass flow rate
%%% If you get an error here - check whether the length of mdot_fuel_SU is
%%% equal to 10 - otherwise, update the fourth row in the following
%%% computation accordingly by changing the 10 to the length of
%%% mdot_fuel_SU.

av_mdot = [zeros(length(n_tsteps), 1), ... %1st column (what is it?)
  ...
  [zeros(total_nodes, 1); ...
  strrep(sol_select(total_nodes+1:end-total_nodes).', 2*ones(1, 10), mdot_fuel_SU.').'; ...
  zeros(total_nodes, 1)], ... 2nd column ->Startup
  ...
  zeros(length(n_tsteps), 1), ... 3rd column
  ...
  (Weight1 * [zeros(total_nodes, 1); ...
  fuel_map(from_state_map(total_nodes+1:end-total_nodes)); ...
  zeros(total_nodes, 1)] ...
  +Weight2 * [zeros(total_nodes, 1); ...
  fuel_map(to_state_map(total_nodes+1:end-total_nodes)); ...
  zeros(total_nodes, 1)] ... 4th column
  ) / (Weight1 + Weight2)];
%% Compute cost of fuel.
% Multiply the average mass flow in each transition by the transition
% length. Then multiply it by dt and the price of the mass.

fuel_cost = bsxfun(@times, av_mdot, double(n_tsteps)) .* dt .* price_kg_f; %fuel cost is totally correct, need to choose right column!
SU_SD_penalties = [zeros(total_nodes, 1); ...
  Startup_cost * ones(length(n_tsteps)-2*total_nodes, 1); ...
  zeros(total_nodes, 1)];


final_cost = fuel_cost;
final_cost(:, 2:3) = bsxfun(@plus, final_cost(:, 2:3), SU_SD_penalties);

decided_costs = final_cost(sub2ind(size(final_cost), ...
  1:size(final_cost, 1), ...
  sol_select.') ...
  ).' +Transition_Penalty * is_transition;
% decided_costs_1 = 0.01*transition_penalty;