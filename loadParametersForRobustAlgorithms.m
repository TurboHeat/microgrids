%% Constants:
% Fuel
[PRICE_kg_f, HEAT_TARIFF] = NATURAL_GAS_PARAMS();
NUM_FUEL_PRICES = numel(PRICE_kg_f);
% Buildings
BUILDINGS = BuildingType(1:4);
NUM_BUILDINGS = numel(BUILDINGS);
FUEL_INDEX = repelem(1:NUM_FUEL_PRICES, 1, NUM_FUEL_PRICES*NUM_BUILDINGS).';

SECONDS_PER_MINUTE = 60;
MINUTES_PER_HOUR = 60;
SECONDS_PER_HOUR = SECONDS_PER_MINUTE * MINUTES_PER_HOUR;
NO_ENVELOPE_ALPHA = 0;
NODES_CONNECTED_TO_ARTIFICIAL_START = 1;
NO_STATE_TRANSITION_PENALTY = 0;
% NO_DEMAND = zeros(T,1); - defined later, as T is only defined later.

%% Initialization
stepsPerHour = SECONDS_PER_HOUR / timeStepSize; % number of time steps in 1h

%% Loading
tmp = load(fullfile(kwargs.dataPath, 'graph_24h.mat'));
g = tmp.g; SV2State = tmp.svToStateNumber; clear tmp;

[elecTariffs, demands_estimate, demands_true, NUM_WINDOWS] = getTariffsAndDemands(BUILDINGS);
nTotalNodes = size(SV2State, 1); %(smax-smin)*(vmax-vmin+1)+1
% Engine maps:
[POWER_MAP, HEAT_MAP, FUEL_MAP, MDOT_FUEL_SU] = loadEngineMaps(nTotalNodes, kwargs.dataPath);

%% Preliminary computations
T = getNumTimesteps(timeStepSize, endTime);
IndicesInHour = T/endTime;
RoundHourIndices = IndicesInHour*(1:endTime);
nF = numel(FUEL_INDEX);
NO_DEMAND = zeros(T,1);

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

%% Preallocation
nPrices = numel(PRICE_kg_f);
nScenarios = nPrices*NUM_BUILDINGS;

power_MGT = zeros(T,1);
heat_MGT = zeros(T,1);
mdot_MGT = zeros(T,1);


function [elecTariffs, demands_averaged, demands_true, NUM_WINDOWS] = getTariffsAndDemands(buildings)
% We save the results for the most common way to call this function
EXISTING_DATA_PATH = "../Data/tarriffs_and_demands.mat";
if isequal(buildings, 1:4) && isfile(EXISTING_DATA_PATH)
  load(EXISTING_DATA_PATH, ...
    'elecTariffs', 'demands_averaged', 'demands_true', 'NUM_WINDOWS');
  return
end
%% 
[CHP_averaged, CHP_nonAveraged, NUM_WINDOWS] = LOAD_DEMAND_DATASETS();

%% Get all tariff & demand combinations
nBuildings = numel(buildings);
elecTariffs = cell(nBuildings, NUM_WINDOWS);
demands_averaged = cell(nBuildings, NUM_WINDOWS);
demands_true = cell(nBuildings, NUM_WINDOWS);

for b = 1:nBuildings
  chp = CHP_averaged(b);
  for d = 1:NUM_WINDOWS
    do = chp.next(); demands_averaged{b,d} = do; % demand object
    elecTariffs{b,d} = getSeasonalElectricityTariff(b, do.timeEnd);
  end
end
for b = 1:nBuildings
  chp = CHP_nonAveraged(b);
  for d = 1:NUM_WINDOWS
    do = chp.next(); demands_true{b,d} = do; % demand object
  end
end
% Unbox:
demands_averaged = reshape([demands_averaged{:}], nBuildings, NUM_WINDOWS).';
demands_true = reshape([demands_true{:}], nBuildings, NUM_WINDOWS).';
elecTariffs = reshape([elecTariffs{:}], nBuildings, NUM_WINDOWS).';

% Make sure we don't need to repeat this computation
save(EXISTING_DATA_PATH, 'elecTariffs', 'demands_averaged', 'demands_true', 'NUM_WINDOWS');
end

function [chp_averaged,chp_notAveraged, nWindows] = LOAD_DEMAND_DATASETS()
% The code below creates 2-week averaging windows for the 5 building types, where
% the first window is [01-Jan-2004 00:00:00, 15-Jan-2004 00:00:00] (because we don't have
% data from the end of 2003 and we don't want to use the end of 2004 as a substitute).

chp_averaged = [CHP2004Provider("../Data/RefBldgLargeHotelNew2004_v1.3_7.1_4A_USA_MD_BALTIMORE.csv", 'windowSize', CHP2004Provider.DEFAULT_WINDOW*CHP2004Provider.DATAPOINTS_PER_DAY),...
       CHP2004Provider("../Data/RefBldgFullServiceRestaurantNew2004_v1.3_7.1_4A_USA_MD_BALTIMORE.csv", 'windowSize', CHP2004Provider.DEFAULT_WINDOW*CHP2004Provider.DATAPOINTS_PER_DAY),...
       CHP2004Provider("../Data/RefBldgSmallHotelNew2004_v1.3_7.1_4A_USA_MD_BALTIMORE.csv", 'windowSize', CHP2004Provider.DEFAULT_WINDOW*CHP2004Provider.DATAPOINTS_PER_DAY),...
       CHP2004Provider("../Data/USA_NY_New.York-Central.Park.725033_TMY3_HIGH.csv", 'windowSize', CHP2004Provider.DEFAULT_WINDOW*CHP2004Provider.DATAPOINTS_PER_DAY),...
       CHP2004Provider("../Data/RefBldgHospitalNew2004_v1.3_7.1_4A_USA_MD_BALTIMORE.csv", 'windowSize', CHP2004Provider.DEFAULT_WINDOW*CHP2004Provider.DATAPOINTS_PER_DAY)];     

chp_notAveraged = [CHP2004Provider("../Data/RefBldgLargeHotelNew2004_v1.3_7.1_4A_USA_MD_BALTIMORE.csv", 'windowSize', 1*CHP2004Provider.DATAPOINTS_PER_DAY),...
       CHP2004Provider("../Data/RefBldgFullServiceRestaurantNew2004_v1.3_7.1_4A_USA_MD_BALTIMORE.csv", 'windowSize', 1*CHP2004Provider.DATAPOINTS_PER_DAY),...
       CHP2004Provider("../Data/RefBldgSmallHotelNew2004_v1.3_7.1_4A_USA_MD_BALTIMORE.csv", 'windowSize', 1*CHP2004Provider.DATAPOINTS_PER_DAY),...
       CHP2004Provider("../Data/USA_NY_New.York-Central.Park.725033_TMY3_HIGH.csv", 'windowSize', 1*CHP2004Provider.DATAPOINTS_PER_DAY),...
       CHP2004Provider("../Data/RefBldgHospitalNew2004_v1.3_7.1_4A_USA_MD_BALTIMORE.csv", 'windowSize', 1*CHP2004Provider.DATAPOINTS_PER_DAY)];   
   
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