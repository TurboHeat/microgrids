function [Output] = runRobustLinftyAlgorithm(timeStepSize, endTime, kwargs)
% Analysis of the solution found by the shortestpath solver for the robust
% Linfty uncertainty set.
%
% Shortest path solvers are used according to algorithms in the paper [].
% We look for the path which minimizes the worst-case cost over all
% possible demand profiles modeled. As the price is higher when the
% demand is higher, the "worst-case" considered by the algorithm is
% achieved when the demand is the maximal possible. For that reason,
% we only consider the case in which the demand is equal to 
% "mean + alpha*std", and not the case in which it is equal to
% "mean - alpha*std".
%
% For alpha = 0, we get the nominal algorithm from (Rist et al., 2017)
%% Handling inputs:
arguments
  timeStepSize (1,1) double {mustBePositive} = 15 % length of a timestep in [s]
  endTime (1,1) double {mustBePositive} = 24; % final time in [h]
                            
  kwargs.alpha (1,1) double = 0.45; 
  % Parameter for first uncertainty set. 
  % In equation (6) in the paper, we take:
  %   Delta_P(t) = alpha_Linfty*std_P(t) 
  %   Delta_H(t) = alpha_Linfty*std_H(t).
  
  kwargs.PriceIndex (1,1) double {mustBePositive} = 1;
  kwargs.BuildingType (1,1) double {mustBePositive} = 1;
    
  kwargs.dataPath (1,1) string = "../Data"
  kwargs.transitionPenalty (1,1) double = 0.01;
end
% Unpack kwargs:
transitionPenalty = kwargs.transitionPenalty;
iP = kwargs.PriceIndex;
iB = kwargs.BuildingType;

%% Load Parameters 
loadParametersForRobustAlgorithms; % this is a script - VERY BAD PRACTICE!!!
% Calling this in a loop is even worse, because many files are read from the hard-drive inside

%% Initialize Variables
% Output Data
Output.Power_Generation = zeros(NUM_WINDOWS,endTime); %Does not save the power generation at all times, but just at ``XX:00" times.
Output.Heat_Generation = zeros(NUM_WINDOWS,endTime);
Output.Fuel_Consumption = zeros(NUM_WINDOWS,endTime);
Output.EstimatedCost = zeros(NUM_WINDOWS,1);
Output.TrueCost = zeros(NUM_WINDOWS,1);
Output.AlgorithmType = 1; % 0 - Nominal (not robust).
                            % 1 - L_infty robust uncertainty set.
                            % 2 - Mixed robust uncertainty set.

AlgorithmParameters.alpha = kwargs.alpha;
Output.AlgorithmParameters{1} = AlgorithmParameters;

%% Run Algorithm
% NUM_WINDOWS = 3; %Debug

fuelPrice = PRICE_kg_f(iP);
heatTariff = HEAT_TARIFF(iP);
for iW = 1:NUM_WINDOWS
    d = demands_estimate(iW, iB);
    mElec = 1e3*d.valMean(:,1); %1e3* - conversion from kWh to W.
    mHeat = 1e3*d.valMean(:,2);
    sElec = 1e3*d.valStd(:,1);
    sHeat = 1e3*d.valStd(:,2);
    
    decided_costs = assignCostsInternal(...
        sol_select, stateFromMap, stateToMap, nTotalNodes, ...
        mElec, sElec, mHeat, sHeat, kwargs.alpha, ...
        elecTariffs(iW, iB), heatTariff, fuelPrice,...
        POWER_MAP, HEAT_MAP, FUEL_MAP, MDOT_FUEL_SU, ...
        transitionPenaltyFlag, transitionPenalty, ...
        timeFrom, nTimesteps, stepsPerHour, timeStepSize);
    
    g = digraph(g.Edges.EndNodes(:,1), g.Edges.EndNodes(:,2), decided_costs);
    [path_MGT, path_cost, path_edge] = shortestpath(g, 1, max(g.Edges.EndNodes(:,2)), 'Method', 'acyclic');
    [power_MGT, heat_MGT, mdot_MGT] = extractPath(path_MGT, POWER_MAP, HEAT_MAP, FUEL_MAP, SV_states);
    
    
    %% Check Performance of the Schedule on True Demand
    % Replace old cost calculation method, which uses a ZOH for the demands
    % and generations, with a new edge-based method that linearly
    % interpolates the demand and generation.
    
    d = demands_true(iW, iB);
    
    true_cost = sum(assignCostsInternal(...
        sol_select(path_edge), stateFromMap(path_edge), stateToMap(path_edge), NODES_CONNECTED_TO_ARTIFICIAL_START, ...
        1e3*d.valMean(:,1), 0*d.valMean(:,1), 1e3*d.valMean(:,2), 0*d.valMean(:,2), 0, ...
        elecTariffs(iW, iB), heatTariff, fuelPrice,...
        POWER_MAP, HEAT_MAP, FUEL_MAP, MDOT_FUEL_SU, ...
        transitionPenaltyFlag(path_edge), NO_STATE_TRANSITION_PENALTY, ...
        timeFrom(path_edge), nTimesteps(path_edge), stepsPerHour, timeStepSize));
    
    %% "Parse" Output Data
    Output.Power_Generation(iW,:) = power_MGT(RoundHourIndices);
    Output.Heat_Generation(iW,:) = heat_MGT(RoundHourIndices);
    Output.Fuel_Consumption(iW,:) = mdot_MGT(RoundHourIndices);
    Output.EstimatedCost(iW) = path_cost;
    Output.TrueCost(iW) = true_cost;
end


end

function [decided_costs] = assignCostsInternal(...
  sol_select, stateFromMap, stateToMap, nTotalNodes, ...
  elecDemandMean, elecDemandStd, heatDemandMean, heatDemandStd, demandStandardEnvelope, ...
  elecTariff, heatTariff, fuelPrice,...
  powerMap, heatMap, fuelMap, mdot_fuel_SU, ...
  transitionPenaltyFlag, transitionPenalty, ...
  time_from, nTimesteps, stepsPerHour, timeStepSize)

% Apply alpha:
% 1kWh = 3.6e6J
elecDemand = elecDemandMean + demandStandardEnvelope * elecDemandStd;
heatDemand = heatDemandMean + demandStandardEnvelope * heatDemandStd;

% Assign edge costs
decided_costs = assignCosts(...
  timeStepSize, powerMap, heatMap, fuelMap, ...
  mdot_fuel_SU, nTotalNodes, sol_select, time_from, nTimesteps,...
  stateFromMap, stateToMap,...
  elecDemand, heatDemand, elecTariff, heatTariff, fuelPrice, ...
  transitionPenaltyFlag, transitionPenalty, stepsPerHour);
end

