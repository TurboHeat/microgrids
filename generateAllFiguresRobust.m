%--------------------------------------------------------------%
% File: generateAllFiguresRobust.m (script)
% Author: Miel Sharf
% Date 05/09/2020
% v1.0
% Description: Analysis of the solution found by the shortest_path solver
% and comparison to robust shortest path solvers. 
% Plots the figures of tariff, demand and solution for all 36 cases. 
% Also compares the solutions on a test case. 
%
% To turn this function into MC, one also needs to iterate on all MC
% cases...
%--------------------------------------------------------------%

%% Startup:
% The startup used to load the variable "decided_costs". This is changed as
% the function assignAllCosts now returns the value...
close all;
clc;

% d_index; 1-3-> Large; 4-6->FSR; 7-9->Small restaurant; % 10-12->Residential
%d_index=mod(index,12)+[mod(index,12)==0]*12;
end_time = 24; %end time in h
dt = 15; %in s
n_lines = 3600 / dt;
T = end_time * n_lines; %number of time steps
t_plot = (1:n_lines * end_time)' ./ n_lines;

%% Parameters for Robust Algorithms
alpha_Linfty = 1; %Parameter for first uncertainty set.
% In equation (6) in the paper, we take Delta_P(t) = alpha_Linfty*std_P(t) and
% Delta_H(t) = alpha_Linfty*std_H(t).

alpha_Mixed = 0.1; %Parameter for second uncertainty set (Algorithm 1 in paper).
alpha_Spikes = 1*sqrt(T); %Parameter for second uncertainty set (Algorithm 1 in paper).
% In equation (7) in the paper, we take Delta_P(t) = alpha_Mixed*std_P(t) and
% Delta_H(t) = alpha_Mixed*std_H(t). Moreover, we take delta_{P,t} = std_P(t),
% delta_{H,t} = std_H(t) and mu_1 = alpha_spikes.

% See paper for more explanation about the meaning of these parameters.

Mixed_approximation_type = 0; %0 = No approximation - might be slow.
                              %1 = Additive approximation - solution costs
                              %no more than epsilon + cost of optimal solution.
                              %2 = Multilicative approximation - solution
                              %costs no more than (1+epsilon) times the
                              %cost of the optimal solution.
                              %3 = Maximum Iterations - Maximum number of
                              %shortest-path runs done by the algorithm.
                              %Runs the additive approximation scheme with
                              %corresponding epsilon.
                              
Epsilon = 1e-2;
MaxIterations = 100;

%%
joule2kWh = 1 / 3.6e6; %1kWh=3.6e6J
fuel_index = [1 * ones(12, 1); 2 * ones(12, 1); 3 * ones(12, 1)];

%% Initialize all vectors for the loop
power_nom_MGT = zeros(T, numel(fuel_index));
heat_nom_MGT = zeros(T, numel(fuel_index));
mdot_nom_MGT = zeros(T, numel(fuel_index));
path_nom_MGT = cell(1,numel(fuel_index));
path_cost_nom = zeros(1, numel(fuel_index));
path_edge_nom = cell(1,numel(fuel_index));

power_robust_Linfty_MGT = zeros(T, numel(fuel_index));
heat_robust_Linfty_MGT = zeros(T, numel(fuel_index));
mdot_robust_Linfty_MGT = zeros(T, numel(fuel_index));
path_robust_Linfty = cell(1,numel(fuel_index));
path_cost_robust_Linfty = zeros(1, numel(fuel_index));
path_edge_robust_Linfty = cell(1,numel(fuel_index));

power_robust_Mixed_MGT = zeros(T, numel(fuel_index));
heat_robust_Mixed_MGT = zeros(T, numel(fuel_index));
mdot_robust_Mixed_MGT = zeros(T, numel(fuel_index));
path_robust_Mixed = cell(1,numel(fuel_index));
path_cost_robust_Mixed = zeros(1, numel(fuel_index));
path_edge_robust_Mixed = cell(1,numel(fuel_index));

%Cost of schedule when the true demand is revealed.
nom_cost = zeros(1,numel(fuel_index));
robust_Linfty_cost = zeros(1,numel(fuel_index));
robust_mixed_cost = zeros(1,numel(fuel_index));
%% Compute decided_costs for all algorithms
%%%TODO: Wrap everything up in a for-loop over all days.

%Demand profile for comparing the different solutiosn.
%%TODO - connect to CHP2004Provider or a similar class.
power_demand_for_comparison = zeros(T,numel(fuel_index));
heat_demand_for_comparison = zeros(T,numel(fuel_index));

%%%TODO: Update calls to assignAllCostsMoving with correct interface.
decided_costs_nominal = assignAllCostsMoving(0); %alpha = 0.

decided_costs_robust_Linfty = assignAllCostsMoving(alpha_Linfty);

decided_costs_robust_mixed_nospike = assignAllCostsMoving(alpha_Mixed);
decided_costs_robust_mixed_withspike = assignAllCostsMoving(alpha_Mixed+alpha_Spikes);


%Retrieve the index list of the edges - saves time later when changing
%weights.
g = digraph(state_from, state_to, 1:size(decided_costs_nominal,1));
EdgePermutationMap = g.Edges.Weight;

for i = 1:numel(fuel_index) 
  %% Run nominal algorithm
  g_nominal = digraph(state_from, state_to, decided_costs_nominal(:, i));
  [path_nom_MGT{i}, path_cost_nom(i), path_edge_nom{i}] = shortestpath(g_nominal, 1, max(state_to), 'Method', 'acyclic');
  [power_nom_MGT(:, i), heat_nom_MGT(:, i), mdot_nom_MGT(:, i)] = extractPath(path_nom_MGT{i}, power_map, heat_map, fuel_map, SV_states);
  %% Run robust algorithm with Linfty uncertainty set
  g_robust_Linfty = digraph(state_from, state_to, decided_costs_robust_Linfty(:, i));
  [path_robust_Linfty{i}, path_cost_robust_Linfty(i), path_edge_robust_Linfty{i}] = shortestpath(g_robust_Linfty, 1, max(state_to), 'Method', 'acyclic');
  [power_robust_Linfty_MGT(:, i), heat_robust_Linfty_MGT(:, i), mdot_robust_Linfty_MGT(:, i)] = extractPath(path_robust_Linfty{i}, power_map, heat_map, fuel_map, SV_states);
  %% Run robust algorithm with mixed uncertainty set
  W_spike = decided_costs_robust_mixed_withspike(:,i) - decided_costs_robust_mixed_nospike(:,i);
  max_W_spike = max(W_spike);
  min_W_spike = min(W_spike);
  unique_W_spike = unique(W_spike);
  
  if(Mixed_approximation_type == 3)
      Epsilon = (max_W_spike - min_W_spike) / (max(MaxIterations,2) - 1); %MaxIterations should never be smaller than 2, but if it is, change its value to 2.
  end
  
  switch Mixed_approximation_type
      case 0
          Thresholds = unique_W_spike;
      case {1,3}
          Thresholds = [min_W_spike:Epsilon:max_W_spike,max_W_spike]'; %Add max(W_spike) at the end to assure that the maximum is also taken into account
      case 2
          Thresholds = exp([log(min_W_spike):log(1+Epsilon):log(max_W_spike),log(max_W_spike)])'; %Add max(W_spike) as before. This is a geometric series.
      otherwise
          error('Unsupported Approximation Type');
  end
  if(numel(Thresholds) > numel(unique_W_spike))
      Thresholds = unique_W_spike;
  end 
  V_costs = zeros(numel(Thresholds),1);
  V_paths = cell(numel(Thresholds),1);
  V_edge_paths = cell(numel(Thresholds),1);
  g_Mixed = digraph(state_from, state_to, decided_costs_robust_mixed_nospike(:, i));
  for j=1:numel(Thresholds)
      Weights = (W_spike <= Thresholds(j)).*decided_costs_robust_mixed_nospike(:, i) + (W_spike > Thresholds(j)).*inf;
      g_Mixed.Edges.Weight = Weights(EdgePermutationMap);
      [V_paths{j}, path_cost, edge_path] = shortestpath(g_Mixed, 1, max(state_to), 'Method', 'acyclic');
      V_costs(j) = path_cost + max(W_spike(edge_path));
      V_edge_paths{j} = edge_path;
  end
  
  [path_cost_robust_Mixed,index_path_Mixed] = min(V_costs);
  path_robust_Mixed = V_paths(index_path_Mixed);
  path_edge_robust_Mixed = V_edge_paths(index_path_Mixed);
  [power_robust_Mixed_MGT(:, i), heat_robust_Mixed_MGT(:, i), mdot_robust_Mixed_MGT(:, i)] = extractPath(path_robust_Mixed, power_map, heat_map, fuel_map, SV_states);
  %% Plot solutions
  %Electricity rate
  d_index = mod(i, 12) + (mod(i, 12) == 0) * 12;
  figure(1);
  plot(t_plot, elec_tariff(:, d_index));
  xlabel('Time (h)');
  ylabel('Rate ($/kWh)');
  title('Applicable electric rate');
  
  %Power generation and demand
  figure(2);
  plot(t_plot,nominal_power_demand(:,i)/1e3,'b','LineWidth', 2);
  hold all; grid on;
  plot(t_plot,power_nom_MGT(:,i)/1e3, t_plot, power_robust_Linfty_MGT(:,i)/1e3, t_plot,power_robust_Mixed_MGT(:,i) ,'LineWidth', 2); %in kW
%   plot(t_plot,wc_power_demand(:,i)/1e3,'b-.','LineWidth', 2);
%   plot(t_plot,2*nominal_power_demand(:,i)/1e3-wc_power_demand(:,i)/1e3,'b-.','LineWidth', 2);
  legend('Power demand','Nominal Commitment', 'Robust (L_\infty) Commitement','Robust (Mixed) Commitement');
  xlim([0 25]); ylim([0 Inf]); grid on;
  xlabel('Time [h]'); ylabel('Demand and commitment [kW]'); hold off;
  
  %Heat generation and demand
  figure(3);
  plot(t_plot,nominal_heat_demand(:,i)/1e3,'b','LineWidth', 2);
  hold all; grid on;
  plot(t_plot,heat_nom_MGT(:,i)/1e3, t_plot, heat_robust_Linfty_MGT(:,i)/1e3, t_plot,heat_robust_Mixed_MGT(:,i) ,'LineWidth', 2); %in kW
%   plot(t_plot,wc_heat_demand(:,i)/1e3,'b-.','LineWidth', 2);
%   plot(t_plot,2*nominal_heat_demand(:,i)/1e3-wc_heat_demand(:,i)/1e3,'b-.','LineWidth', 2);
  legend('Heat demand','Nominal Commitment', 'Robust Commitement');
  xlim([0 25]); ylim([0 Inf]); grid on;
  xlabel('Time [h]'); ylabel('Demand and commitment [kW]'); hold off;
  
  %% Compare schedules for a single scenario.
  
  % Replace old cost calculation method, which uses a ZOH for the demands
  % and generations, with a new edge-based method that linearly
  % interpolates the demand and generation.
  
  %%%TODO: If needed (ask Beni) - separate costs to bought electricity,
  %%%sold electricity, bought fuel, and bought heat. Do this by looking at
  %%%the individual components which generate the edge cost in
  %%%assignCosts.m
  
  nom_cost(i) = sum(assignCosts(double(dt), power_map, heat_map, fuel_map, mdot_fuel_SU, 1,... %1 - only one node in this path is connected to the artificial start/end node.
                    sol_select(path_edge_nom), time_from(path_edge_nom), n_tsteps(path_edge_nom),...
                    from_state_map(path_edge_nom), to_state_map(path_edge_nom),power_demand_for_comparison(:, d_index), heat_demand_for_comparison(:, d_index),...
                    elec_tariff(:, d_index), heat_tariff(fuel_index(i)), price_kg_f(fuel_index(i)), transition_penalty(path_edge_nom),0));  
                
  robust_Linfty_cost(i) = sum(assignCosts(double(dt), power_map, heat_map, fuel_map, mdot_fuel_SU, 1,... %1 - only one node in this path is connected to the artificial start/end node.
                    sol_select(path_edge_robust_Linfty), time_from(path_edge_robust_Linfty), n_tsteps(path_edge_robust_Linfty),...
                    from_state_map(path_edge_robust_Linfty), to_state_map(path_edge_robust_Linfty),power_demand_for_comparison(:, d_index), heat_demand_for_comparison(:, d_index),...
                    elec_tariff(:, d_index), heat_tariff(fuel_index(i)), price_kg_f(fuel_index(i)), transition_penalty(path_edge_robust_Linfty),0));  
                
  robust_mixed_cost(i) = sum(assignCosts(double(dt), power_map, heat_map, fuel_map, mdot_fuel_SU, 1,... %1 - only one node in this path is connected to the artificial start/end node.
                    sol_select(path_edge_robust_Mixed), time_from(path_edge_robust_Mixed), n_tsteps(path_edge_robust_Mixed),...
                    from_state_map(path_edge_robust_Mixed), to_state_map(path_edge_robust_Mixed),power_demand_for_comparison(:, d_index), heat_demand_for_comparison(:, d_index),...
                    elec_tariff(:, d_index), heat_tariff(fuel_index(i)), price_kg_f(fuel_index(i)), transition_penalty(path_edge_robust_Mixed),0));               

end
% disp(['All data saved to folder ', savepath]);

%% Save economic data
% savepath = '..\Data';
% save([savepath, 'econ_data.mat'], 'elec_tariff', 'heat_tariff', ...
%   'power_demand', 'heat_demand', 'power_MGT', 'heat_MGT', 'MGT_cost', ...
%   'electric_energy_kWh', 'heat_energy_kWh', 'base_electricity_charge', ...
%   'base_heat_charge', 'total_charge', 'bought_elec', 'sold_energy', ...
%   'bought_heat', 'bought_fuel', 'MGT_cost', 'savings', 'new_demand', ...
%   'FC', 'ut_PDC', 'ut_IDC', 'MGT_IDC', 'MGT_PDC', 'mdot_MGT');