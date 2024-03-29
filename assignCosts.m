%--------------------------------------------------------------%
% File: assignCosts.m (script)
% Author: Miguel Dias
% Date 10/08/16
% v1.0
% Description: Assign costs to edges for the shortest path solver/
% Included transition penalty to deal with oscillations
%--------------------------------------------------------------%
function [decided_costs] = assignCosts(...
  dt, power_map, heat_map, fuel_map, mdot_fuel_SU, total_nodes, sol_select, ...
  time_from, tStepsRequired, from_state_map, to_state_map, power_demand, heat_demand, ...
  elec_tariff, heat_tariff, price_kg_f, transition_penalty_indicator, transition_penalty,...
  timestepsPerDay)

%% Definitions:
STARTUP_COST = 3.75; % SU and SD consts are the same
ZERO_COST = zeros(total_nodes, 1); % The first and last total_nodes transitions have zero cost -> transitions "from start" and "to end"
Joule_TO_kWh = 1 / 3.6e6; % 1kWh = 3.6e6J
HOURS_PER_DAY = 24;
W1 = 1; % Weight of the source node's price in the cost of a transition
W2 = 1; % Weight of the destination node's price in the cost of a transition
%% Compute the four possible prices (Off-Off, Start up, Shut down, Remaining transitions)
%Compute average demands, tariff and mdot - these are general for the
%four cases

tAdvance = tStepsRequired( total_nodes+1:end-total_nodes);
timestepOfDay = time_from( total_nodes+1:end-total_nodes);
sh = floor(  timestepOfDay             / timestepsPerDay) + 1; % transition start's hour-of-day index (1...25)
eh = floor( (timestepOfDay + tAdvance) / timestepsPerDay) + 1; % transition end's hour-of-day index (1...25)

av_P_demand = (... 
     W1*[ZERO_COST; power_demand( sh ); ZERO_COST] ...
   + W2*[ZERO_COST; power_demand( eh ); ZERO_COST] ...
     )/(W1+W2);

av_H_demand = (... 
     W1*[ZERO_COST; heat_demand( sh ); ZERO_COST] ...
   + W2*[ZERO_COST; heat_demand( eh ); ZERO_COST] ...
     )/(W1+W2); 

% NOTE: here we force the electric tariff to be the same at midnight on both sides, 
% which is slightly incorrect when changing seasons (i.e. the next day has another tariff).
% This helps us avoid having to test each day for seasonal tariff transitions.
av_tariff = (... 
     W1*[ZERO_COST; elec_tariff.tariff( sh ); ZERO_COST] ...
   + W2*[ZERO_COST; elec_tariff.tariff( mod(eh-1,HOURS_PER_DAY)+1 ); ZERO_COST] ...
     )/(W1+W2);

warning off 'MATLAB:strrep:InvalidInputType'

% average mass flow rate
%%% If you get an error here - check whether the length of mdot_fuel_SU is
%%% equal to 10 - otherwise, update the fourth row in the following
%%% computation accordingly by changing the 10 to the length of
%%% mdot_fuel_SU.
nMDF = numel(mdot_fuel_SU);
av_mdot = [zeros(numel(tStepsRequired), 1), ... %1st column (what is it?)
  ...
  [ZERO_COST; ...
  strrep(sol_select(total_nodes+1:end-total_nodes).', 2*ones(1, nMDF), mdot_fuel_SU.').'; ...
  ZERO_COST], ... 2nd column -> Startup
  ...
  zeros(numel(tStepsRequired), 1), ... 3rd column
  ...
  (W1 * [ZERO_COST; ...
         fuel_map(from_state_map(total_nodes+1:end-total_nodes)); ...
         ZERO_COST] + ...
   W2 * [ZERO_COST; ...
         fuel_map(to_state_map(total_nodes+1:end-total_nodes)); ...
         ZERO_COST] ... 4th column
  ) / (W1 + W2)];

% Compute average power and heat produced
av_H_produced = [zeros(numel(tStepsRequired), 3), ...
  (W1 * [ZERO_COST; ...
         heat_map(from_state_map(total_nodes+1:end-total_nodes)); ...
         ZERO_COST] + ...
   W2 * [ZERO_COST; ...
         heat_map(to_state_map(total_nodes+1:end-total_nodes)); ...
         ZERO_COST] ...
  ) / (W1 + W2)];
%%

av_P_produced = [zeros(numel(tStepsRequired), 3), ...
  (W1 * [ZERO_COST; ...
         power_map(from_state_map(total_nodes+1:end-total_nodes)); ...
         ZERO_COST] + ...
   W2 * [ZERO_COST; ...
         power_map(to_state_map(total_nodes+1:end-total_nodes)); ...
         ZERO_COST] ...
  ) / (W1 + W2)];

%% Compute cost parcels
%First multiply tariff, 2nd to multiply timesteps and
%3rd to compute power balance

elec_cost = (av_P_demand - av_P_produced) .* (tStepsRequired .* av_tariff * (dt * Joule_TO_kWh));

delta_heat = av_H_demand - av_H_produced;
heat_cost = delta_heat .* (delta_heat >= 0) .* (tStepsRequired * (dt * Joule_TO_kWh * heat_tariff));

fuel_cost = av_mdot .* (tStepsRequired .* (dt * price_kg_f)); %fuel cost is totally correct, need to choose right column!
SU_SD_penalties = [... SU_SD => startup / shutdown
  ZERO_COST; ...
  STARTUP_COST * ones(numel(tStepsRequired)-2*total_nodes, 1); ...
  ZERO_COST];

final_cost = fuel_cost + elec_cost + heat_cost;
final_cost(:, 2:3) = final_cost(:, 2:3) + SU_SD_penalties;

decided_costs = final_cost(sub2ind(size(final_cost), ...
  1:size(final_cost, 1), ...
  sol_select.') ...
  ).' +transition_penalty * transition_penalty_indicator;