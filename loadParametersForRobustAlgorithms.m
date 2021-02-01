function [FUEL_MAP, HEAT_MAP, HEAT_TARIFF, MDOT_FUEL_SU, NODES_CONNECTED_TO_ARTIFICIAL_START, ...
 NUM_WINDOWS, POWER_MAP, PRICE_kg_f, RoundHourIndices, SV_states, demands_estimate, ...
 demands_true, elecTariffs, g, nTimesteps, nTotalNodes, sol_select, stateFromMap, ...
 stateToMap, stepsPerHour, timeFrom, transitionPenaltyFlag] = ...
 loadParametersForRobustAlgorithms(endTime, iB, dataPath, psf, timeStepSize)
%% Constants:
% Fuel
[PRICE_kg_f, HEAT_TARIFF] = NATURAL_GAS_PARAMS();

SECONDS_PER_MINUTE = 60;
MINUTES_PER_HOUR = 60;
SECONDS_PER_HOUR = SECONDS_PER_MINUTE * MINUTES_PER_HOUR;
NODES_CONNECTED_TO_ARTIFICIAL_START = 1;
% NO_DEMAND = zeros(T,1); - defined later, as T is only defined later.

%% Initialization
stepsPerHour = SECONDS_PER_HOUR / timeStepSize; % number of time steps in 1h

%% Loading
% Possible transitions:
tmp = load(fullfile(dataPath, 'graph_24h.mat'));
g = tmp.g; SV2State = tmp.svToStateNumber; clear tmp;
nTotalNodes = size(SV2State, 1); %(smax-smin)*(vmax-vmin+1)+1
% Engine maps:
[POWER_MAP, HEAT_MAP, FUEL_MAP, MDOT_FUEL_SU] = loadEngineMaps(nTotalNodes, dataPath);
% Tariffs and demands
maxPower = max(POWER_MAP)/1000; %[W] -> [kW]
[elecTariffs, demands_estimate, demands_true, NUM_WINDOWS] = getTariffsAndScaledDemands(dataPath, iB, maxPower.*psf);

%% Preliminary computations
T = getNumTimesteps(timeStepSize, endTime);
IndicesInHour = T/endTime;
RoundHourIndices = IndicesInHour*(1:endTime);

[SV_states, timeFrom, nTimesteps, stateFromMap, stateToMap] = ...
  transformAdjacencyMatrix(g.Edges.EndNodes(:,1), g.Edges.EndNodes(:,2), SV2State);

% Select the correct column with the final price
sol_select = [~SV_states(stateFromMap, 1) & ~SV_states(stateToMap, 1), ... % Off-off
              ~SV_states(stateFromMap, 1) & SV_states(stateToMap, 1), ... % Start up
              ~SV_states(stateToMap, 1), ... % Shut down
              true(numel(nTimesteps), 1)]; % Remaining transitions
[~, sol_select] = max(sol_select, [], 2);

% Assign a small penalty to every input (s,v) change
transitionPenaltyFlag = [zeros(nTotalNodes, 1); ...
  ~(SV_states(stateFromMap(nTotalNodes+1:end-nTotalNodes), 1) == SV_states(stateToMap(nTotalNodes+1:end-nTotalNodes), 1) & ... %checks equality of S values
    SV_states(stateFromMap(nTotalNodes+1:end-nTotalNodes), 2) == SV_states(stateToMap(nTotalNodes+1:end-nTotalNodes), 2));
   zeros(nTotalNodes,1)]; % checks equality of V values
end % function

function [elecTariffs, demands_averaged, demands_true, NUM_WINDOWS] = getTariffsAndScaledDemands(dataPath, buildings, newMaxPower)
arguments
  dataPath (1,1) string
  buildings (1,:) BuildingType = BuildingType(1:4)
  newMaxPower (1,:) double = NaN % technically, several are supported, but we're only using one at a time
end
% Try to load precomputed results
if isnan(newMaxPower)
  suffix = sprintf('B%1u', buildings);
elseif isscalar(newMaxPower)
  suffix = sprintf('B%1u_x%4.2f', buildings, newMaxPower);
else
  suffix = sprintf('B%1u_xM', buildings);
end
resultsFilename = fullfile(dataPath, sprintf('tarriffs_and_demands_%s.mat', suffix));
if isfile(resultsFilename)
  load(resultsFilename, ...
    'elecTariffs', 'demands_averaged', 'demands_true', 'NUM_WINDOWS');
  return
end
%% 
[CHP_averaged, CHP_nonAveraged, NUM_WINDOWS] = loadDemandDatasets(dataPath, buildings);

%% Get all tariff & demand combinations
nSc = numel(newMaxPower);
nBuildings = numel(buildings);
elecTariffs = cell(nBuildings, NUM_WINDOWS);
demands_averaged = cell(nBuildings, NUM_WINDOWS, nSc);
demands_true = cell(nBuildings, NUM_WINDOWS, nSc);
% Get tariffs ONLY:
for b = 1:nBuildings
  chp = CHP_averaged(b);
  for d = 1:NUM_WINDOWS
    do = chp.next();
    elecTariffs{b,d} = getSeasonalElectricityTariff(buildings(b), do.timeEnd);
  end
  % Rewind the object
  chp.fastForward(chp.timestamps(1), 'first');  
end

for s = 1:nSc
  % Perform scaling:
  arrayfun(@(x)x.rescalePower(newMaxPower(s)), CHP_averaged);
  arrayfun(@(x)x.fastForward(chp.timestamps(1), 'first'), CHP_nonAveraged);
  arrayfun(@(x)x.rescalePower(newMaxPower(s)), CHP_nonAveraged);
  arrayfun(@(x)x.fastForward(datetime(2004, 1, 15), 'last'), CHP_nonAveraged);
  % Get moving-averaged demands:
  for b = 1:nBuildings
    chp = CHP_averaged(b);
    for d = 1:NUM_WINDOWS
      do = chp.next(); 
      demands_averaged{b,d,s} = do; % demand object
    end
  end
  % Get true demands:
  for b = 1:nBuildings
    chp = CHP_nonAveraged(b);
    for d = 1:NUM_WINDOWS
      do = chp.next(); 
      demands_true{b,d,s} = do; % demand object
    end
  end
end
% Unbox:
demands_averaged = reshape([demands_averaged{:}], nBuildings, NUM_WINDOWS).';
demands_true = reshape([demands_true{:}], nBuildings, NUM_WINDOWS).';
elecTariffs = reshape([elecTariffs{:}], nBuildings, NUM_WINDOWS).';

% Make sure we don't need to repeat this computation
save(resultsFilename, 'buildings', 'elecTariffs', 'demands_averaged', 'demands_true', ...
  'NUM_WINDOWS', 'newMaxPower');
end

function [chp_averaged,chp_notAveraged, nWindows] = loadDemandDatasets(dataPath, buildings)
% The code below creates 2-week averaging windows for the requested building types, where
% the first window is [01-Jan-2004 00:00:00, 15-Jan-2004 00:00:00] (because we don't have
% data from the end of 2003 and we don't want to use the end of 2004 as a substitute).
SINGLE_DAY_WINDOW = 1;
PPD = CHP2004Provider.DATAPOINTS_PER_DAY;

% DO NOT REORDER THE CSV LOADING (this will require updating the BuildingType enumeration accordingly)
CSV_NAMES = ["RefBldgLargeHotelNew2004_v1.3_7.1_4A_USA_MD_BALTIMORE.csv";
             "RefBldgFullServiceRestaurantNew2004_v1.3_7.1_4A_USA_MD_BALTIMORE.csv";
             "RefBldgSmallHotelNew2004_v1.3_7.1_4A_USA_MD_BALTIMORE.csv";
             "USA_NY_New.York-Central.Park.725033_TMY3_HIGH.csv";
             "RefBldgHospitalNew2004_v1.3_7.1_4A_USA_MD_BALTIMORE.csv"];
csvPaths = fullfile(dataPath, CSV_NAMES);

for iC = numel(buildings):-1:1
  iB = buildings(iC);
  chp_averaged(iC) = CHP2004Provider(csvPaths(iB), 'windowSize', CHP2004Provider.DEFAULT_WINDOW*PPD);
  chp_notAveraged(iC) = CHP2004Provider(csvPaths(iB), 'windowSize', SINGLE_DAY_WINDOW*PPD);
end

startDay = datetime(2004, 1, 15); % skip the first two weeks
arrayfun(@(x)x.fastForward(startDay, 'last'), chp_averaged);
arrayfun(@(x)x.fastForward(startDay, 'last'), chp_notAveraged);
% Count how many times "next" can be called:
nextCnt = 0; % Should be equal to: 366-14+1 = 353 (1 day extra into the new year)
while chp_averaged(1).hasNext()
  nextCnt = nextCnt + 1;
  [~] = chp_averaged(1).next();
end
% Rewind again:
arrayfun(@(x)x.fastForward(startDay, 'last'), chp_averaged);
nWindows = nextCnt;

% chp = struct2array(load('../Data/CHP2004.mat', 'chp'));
% nWindows = struct2array(load('../Data/CHP2004.mat', 'nextCnt'));
end

function [powerMap, heatMap, fuelMap, mdotFuelSU] = loadEngineMaps(expectedTotalNodes, dataPath)
arguments 
  expectedTotalNodes (1,1) double {mustBeInteger, mustBePositive} = NaN
  dataPath (1,1) string = "../Data"
end
%% Load (and optionally rename) mappings
if isfile(fullfile(dataPath, 'CHP.mat')) % Engine similar to Capstone C65 (incl. extrapolation)
  tmp = load(fullfile(dataPath, 'CHP.mat'));
  powerMap = tmp.PowerSurf;
  heatMap = tmp.PheatSurfExt;
  m_dot_f_map = tmp.FFSurfExt;
else
  tmp = load(fullfile(dataPath, 'sv_mappings.mat')); % "Virtual" engine
  powerMap = tmp.power_map;
  heatMap = tmp.heat_map;
  m_dot_f_map = tmp.m_dot_f_map;
end
% Verify that we have the right mappings:
assert( isequal( numel(powerMap), numel(heatMap), numel(m_dot_f_map), expectedTotalNodes-1 ), ...
  'Mismatch between matrix sizes! Try using a different ''sv_mappings'' file.');

% Reshape the power maps
powerMap = [reshape(powerMap.', [], 1); 0];
heatMap = [reshape(heatMap.', [], 1); 0];
fuelMap = [reshape(m_dot_f_map.', [], 1); 0];
mdotFuelSU = sum(m_dot_f_map(:, 1:end-1)+m_dot_f_map(:, 2:end), 2); % Computes mdot of fuel consumed in startup for all the different max states
end

function T = getNumTimesteps(dt, endTime)
  SECONDS_PER_MINUTE = 60;
  MINUTES_PER_HOUR = 60;
  SECONDS_PER_HOUR = SECONDS_PER_MINUTE * MINUTES_PER_HOUR;
  T = endTime * SECONDS_PER_HOUR / dt;    
end