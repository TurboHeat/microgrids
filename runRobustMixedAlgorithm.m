function [Output] = runRobustMixedAlgorithm(timeStepSize, endTime, kwargs)
% Analysis of the solution found by the Algorithm 1 in (), for a L1/Linfty
% mixed uncertainty set. The solution is based on running the shortest path
% algorithm multiple times, for slightly different graphs.
%% Handling inputs:
arguments
  timeStepSize (1,1) double {mustBePositive} = 15 % length of a timestep in [s]
  endTime (1,1) double {mustBePositive} = 24; % final time in [h]
                            
  kwargs.alphaMixed (1,1) double = 0.1; 
  % Parameter for second uncertainty set (Algorithm 1 in the paper).
  
  kwargs.alphaSpikes (1,1) double = 4*sqrt( getNumTimesteps(timeStepSize, endTime) ); 
  % Parameter for second uncertainty set (Algorithm 1 in paper).
  % In equation (7) in the paper, we take Delta_P(t) = alpha_Mixed*std_P(t) and
  % Delta_H(t) = alpha_Mixed*std_H(t). Moreover, we take delta_{P,t} = std_P(t),
  % delta_{H,t} = std_H(t) and mu_1 = alpha_spikes.
  
  kwargs.PriceIndex (1,1) double {mustBePositive} = 1;
  kwargs.BuildingType (1,1) double {mustBePositive} = 1;

  kwargs.approx (1,1) ApproximationType = ApproximationType.MaxIter;
  kwargs.epsilon (1,1) double = 1E-1;
  kwargs.maxIter (1,1) double {mustBeGreaterThanOrEqual(kwargs.maxIter, 2)} = 30;
  
  kwargs.dataPath (1,1) string = "../Data"
  kwargs.transitionPenalty (1,1) double = 0.01;
end
% Unpack kwargs:
transitionPenalty = kwargs.transitionPenalty;
epsilon = kwargs.epsilon;
iP = kwargs.PriceIndex;
iB = kwargs.BuildingType;
%% Load Parameters 
loadParametersForRobustAlgorithms; % this is a script - VERY BAD PRACTICE!!!
% Calling this in a loop is even worse, because many files are read from the hard-drive inside

% NWI = 3; %Debug
NWI = NUM_WINDOWS-1; % Number of Windows of Interest
% Explanation of the "-1" above:
%   The {averaging} windows are used to predict "the next day".
%   The last "next day" happens to be in the next year (01/01/2005).
%   We are not interested in predicting that day even though we have enough data for it, 
%   because we have no ground-truth to compare it with later on.
%=> As a result we use NUM_WINDOWS-1 to reflect this

%% Initialize Variables
% Output Data
Output.Power_Generation = zeros(NWI,endTime); %Does not save the power generation at all times, but just at ``XX:00" times.
Output.Heat_Generation = zeros(NWI,endTime);
Output.Fuel_Consumption = zeros(NWI,endTime);
Output.EstimatedCost = zeros(NWI,1);
Output.TrueCost = zeros(NWI,1);
Output.AlgorithmType = 2; % 0 - Nominal (not robust).
                          % 1 - L_infty robust uncertainty set.
                          % 2 - Mixed robust uncertainty set.
                          
switch Output.AlgorithmType
    case 0
        Output.AlgorithmParameters{1} = [];
    case 1
        AlgorithmParameters.alpha = kwargs.alphaLin;
        Output.AlgorithmParameters{1} = AlgorithmParameters;
    case 2
        AlgorithmParameters.alphaMixed = kwargs.alphaMixed;
        AlgorithmParameters.alphaSpikes = kwargs.alphaSpikes;
        AlgorithmParameters.ApproximationType = kwargs.approx;
        switch kwargs.approx
        case {ApproximationType.Additive, ApproximationType.Multiplicative}
            AlgorithmParameters.epsilon = kwargs.epsilon;
        case ApproximationType.MaxIter
            AlgorithmParameters.maxIterations = kwargs.maxIter;
        case ApproximationType.Exact
            %Do nothing.
            otherwise
          error('Unsupported Approximation Type');
        end
        Output.AlgorithmParameters{1} = AlgorithmParameters;
    otherwise
        error('Unsupported Algorithm Type');
end
Output.AlgorithmParameters{1} = AlgorithmParameters;

%% Run Algorithm
DATES_OF_DATA = (datetime(2004,1,15):datetime(2004,12,31)).'; 
isWeekend = weekday(DATES_OF_DATA) == 7 | weekday(DATES_OF_DATA) == 1;

fuelPrice = PRICE_kg_f(iP);
heatTariff = HEAT_TARIFF(iP);

for iW = 1:NWI
    d = demands_estimate(iW, iB);    
    mElec = 1e3*d.valMean(:,1, 1+isWeekend(iW)); %1e3* - conversion from kWh to W.
    mHeat = 1e3*d.valMean(:,2, 1+isWeekend(iW));
    sElec = 1e3*d.valStd(:,1, 1+isWeekend(iW));
    sHeat = 1e3*d.valStd(:,2, 1+isWeekend(iW));
    
    decided_costs_nospike = assignCostsInternal(...
        sol_select, stateFromMap, stateToMap, nTotalNodes, ...
        mElec, sElec, mHeat, sHeat, kwargs.alphaMixed, ...
        elecTariffs(iW, iB), heatTariff, fuelPrice,...
        POWER_MAP, HEAT_MAP, FUEL_MAP, MDOT_FUEL_SU, ...
        transitionPenaltyFlag, transitionPenalty, ...
        timeFrom, nTimesteps, stepsPerHour, timeStepSize);
    
    decided_costs_withspike = assignCostsInternal(...
        sol_select, stateFromMap, stateToMap, nTotalNodes, ...
        mElec, sElec, mHeat, sHeat, kwargs.alphaMixed + kwargs.alphaSpikes, ...
        elecTariffs(iW, iB), heatTariff, fuelPrice,...
        POWER_MAP, HEAT_MAP, FUEL_MAP, MDOT_FUEL_SU, ...
        transitionPenaltyFlag, transitionPenalty, ...
        timeFrom, nTimesteps, stepsPerHour, timeStepSize);
    
    W_spike = decided_costs_withspike - decided_costs_nospike;
    max_W_spike = max(W_spike);
    min_W_spike = min(W_spike);
    unique_W_spike = unique(W_spike);
    nUWS = numel(unique_W_spike);
    
    if (kwargs.approx == ApproximationType.MaxIter)
        epsilon = (max_W_spike - min_W_spike) / (kwargs.maxIter - 1);
    end
    
    switch kwargs.approx
        case ApproximationType.Exact
            Thresholds = unique_W_spike;
        case {ApproximationType.Additive, ApproximationType.MaxIter}
            Thresholds = [min_W_spike:epsilon:max_W_spike,max_W_spike].'; %Add max(W_spike) at the end to assure that the maximum is also taken into account
        case ApproximationType.Multiplicative
            Thresholds = exp([log(min_W_spike):log(1+epsilon):log(max_W_spike),log(max_W_spike)]).'; %Add max(W_spike) as before. This is a geometric series.
    end
    
    nTrs = numel(Thresholds);
    if nTrs > nUWS
        Thresholds = unique_W_spike;
        nTrs = numel(Thresholds);
    end
    V_costs = zeros(nTrs,1);
    V_paths = cell(nTrs,1);
    V_edge_paths = cell(nTrs,1);
    g = digraph(g.Edges.EndNodes(:,1), g.Edges.EndNodes(:,2), 1:numel(g.Edges.EndNodes(:,1)));
    % Retrieve the index list of the edges - saves time later when changing weights.
    edgePermutationMap = g.Edges.Weight;
    for iT = 1:nTrs
        Weights =  decided_costs_nospike;
        Weights(W_spike > Thresholds(iT)) = inf;
        g.Edges.Weight = Weights(edgePermutationMap);
        [V_paths{iT}, path_cost, edge_path] = shortestpath(g, 1,max(g.Edges.EndNodes(:,2)), 'Method', 'acyclic');
        V_costs(iT) = path_cost + max(W_spike(edge_path));
        V_edge_paths{iT} = edge_path;
    end
    
    [path_cost, index_path] = min(V_costs);
    path_MGT = V_paths{index_path};
    path_edge = V_edge_paths{index_path};
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
        transitionPenaltyFlag(path_edge), transitionPenalty, ...
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

function T = getNumTimesteps(dt, endTime)
  SECONDS_PER_MINUTE = 60;
  MINUTES_PER_HOUR = 60;
  SECONDS_PER_HOUR = SECONDS_PER_MINUTE * MINUTES_PER_HOUR;
  T = endTime * SECONDS_PER_HOUR / dt;    
end