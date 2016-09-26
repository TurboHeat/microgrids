%--------------------------------------------------------------%
% File: assign_costs_fcn.m (script)
% Author: Miguel Dias
% Date 10/08/16
% v1.0
% Description: Assign costs to edges for the shortest path solver/
% Included transition penalty to deal with oscillations
%--------------------------------------------------------------%
function [decided_costs]=assign_costs_fcn(dt,power_map, heat_map, fuel_map,mdot_fuel_SU, total_nodes, sol_select,time_from,n_tsteps,from_state_map,to_state_map, power_demand, heat_demand, elec_tariff, heat_tariff, price_kg_f, transition_penalty)
%% Definitions:
Startup_cost  =3.75; % SU and SD consts are the same
joule2kWh=1/3.6e6; %1kWh=3.6e6J
%% Compute the four possible prices (Off-Off, Start up, Shut down, Remaining transitions)
    %Compute average demands, tariff and mdot - these are general for the
    %four cases
    Weight1 = 1; Weight2 = 1;   
    %av_P_demand=(P_demand(t)+P_demand(t+1))/2
av_P_demand = (... 
     Weight1*[zeros(total_nodes,1)
              power_demand(time_from(total_nodes+1:end-total_nodes))   %first 46 and last 46 transitions have zero cost-> transitions from start and to end
              zeros(total_nodes,1)]...
   + Weight2*[zeros(total_nodes,1)
              power_demand((time_from(total_nodes+1:end-total_nodes)...
            + n_tsteps(total_nodes+1:end-total_nodes)))
              zeros(total_nodes,1)]...
     )/(Weight1+Weight2);
%av_H_demand=(H_demand(t)+H_demand(t+1))/2
 av_H_demand = (... 
     Weight1*[zeros(total_nodes,1)
              heat_demand(time_from(total_nodes+1:end-total_nodes))   
              zeros(total_nodes,1)] ...
   + Weight2*[zeros(total_nodes,1)
              heat_demand((time_from(total_nodes+1:end-total_nodes)...
            + n_tsteps(total_nodes+1:end-total_nodes)))
              zeros(total_nodes,1)]...
     )/(Weight1+Weight2);                   
%av_tariff=(tariff(t)+tariff(t+1))/2
av_tariff = (... 
     Weight1*[zeros(total_nodes,1)
              elec_tariff(time_from(total_nodes+1:end-total_nodes))   
              zeros(total_nodes,1)] ...
   + Weight2*[zeros(total_nodes,1)
              elec_tariff((time_from(total_nodes+1:end-total_nodes)...
            + n_tsteps(total_nodes+1:end-total_nodes)))
              zeros(total_nodes,1)]...
     )/(Weight1+Weight2);

warning off 'MATLAB:strrep:InvalidInputType'
av_mdot=[zeros(numel(n_tsteps),1),... %1st column
         ...
        [zeros(total_nodes,1)
           strrep(sol_select(total_nodes+1:end-total_nodes).',2*ones(1,5),mdot_fuel_SU.').';
           zeros(total_nodes,1)],... 2nd column ->Startup
         ...
         zeros(numel(n_tsteps),1),... 3rd column
          ...
         (Weight1*[zeros(total_nodes,1)
          fuel_map(from_state_map(total_nodes+1:end-total_nodes))   
          zeros(total_nodes,1)] ...
        + Weight2*[zeros(total_nodes,1)
          fuel_map(to_state_map(total_nodes+1:end-total_nodes))
          zeros(total_nodes,1)]... 4th column
    )/(Weight1+Weight2)];  

    % Compute average power and heat produced
av_H_produced = [zeros(numel(n_tsteps),3),...
          (Weight1*[zeros(total_nodes,1)
      heat_map(from_state_map(total_nodes+1:end-total_nodes))   
      zeros(total_nodes,1)] ...
+ Weight2*[zeros(total_nodes,1)
      heat_map(to_state_map(total_nodes+1:end-total_nodes))
      zeros(total_nodes,1)]...
)/(Weight1+Weight2)];   
 
Weight1 = 1; Weight2 = 1;
av_P_produced = [zeros(numel(n_tsteps),3),...
              (Weight1*[zeros(total_nodes,1)
          power_map(from_state_map(total_nodes+1:end-total_nodes))   
          zeros(total_nodes,1)] ...
+ Weight2*[zeros(total_nodes,1)
          power_map(to_state_map(total_nodes+1:end-total_nodes))
          zeros(total_nodes,1)]...
 )/(Weight1+Weight2)];

  
%% Compute cost parcels
%First call to bsxfun to multiply tariff, 2nd to multiply timesteps and
%3rd to compute power balance
elec_cost = bsxfun(@times,...
            bsxfun(@times,...
            bsxfun(@minus,av_P_demand,av_P_produced), double(n_tsteps)),av_tariff)...
            *dt*joule2kWh;
        
delta_heat=bsxfun(@minus,av_H_demand,av_H_produced);
heat_cost=bsxfun(@times,delta_heat.*bsxfun(@ge,delta_heat,0),...
          double(n_tsteps))*dt*joule2kWh*heat_tariff;

fuel_cost=bsxfun(@times,av_mdot,double(n_tsteps)).*dt.*price_kg_f; %fuel cost is totally correct, need to choose right column!
SU_SD_penalties=[zeros(total_nodes,1)
                 Startup_cost*ones(numel(n_tsteps)-2*total_nodes,1)
                 zeros(total_nodes,1)];
             
final_cost=fuel_cost+elec_cost+heat_cost;
final_cost(:,2:3)=bsxfun(@plus,final_cost(:,2:3),SU_SD_penalties);

decided_costs = final_cost(sub2ind(size(final_cost),...
                                   1:size(final_cost,1),...
                                   sol_select.')...
                           ).'+0.01*transition_penalty;
