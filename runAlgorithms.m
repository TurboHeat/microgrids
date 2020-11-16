function [allScenariosData] = runAlgorithms(fuelIdx)
%% Choose Running Algorithms
runBenchmark = 1;

runNominal = 1;

runRobustLinfty = 1;
LinftyAlphas = [0.01:0.01:0.05, 0.1:0.05:0.7];

runRobustMixed = 1;
mixedAlphas = repmat([0.01:0.01:0.05,0.07,0.1,0.13,0.15,0.17,0.2],[1,10]);
mixedSpikeAlphas = [.5*sqrt(5760)*ones(1,11),...
  1*sqrt(5760)*ones(1,11),...
  1.5*sqrt(5760)*ones(1,11),...
  2*sqrt(5760)*ones(1,11),...
  2.5*sqrt(5760)*ones(1,11),...
  3*sqrt(5760)*ones(1,11),...
  3.5*sqrt(5760)*ones(1,11),...
  4*sqrt(5760)*ones(1,11),...
  1*5760*ones(1,11),...
  2*5760*ones(1,11)];

%% Choose Prices and Building Combinations
%i-th scenario has building type BuildingTypes(i) and the index of the gas
%price is PriceIndices(i).
%
%The building type should be in the set {1,2,3,4}, and the price index
%should be in the set {1 ... numFuelPrices}
BUILDINGS_TO_TEST = BuildingType(1:4);
nBuildTypes = numel(BUILDINGS_TO_TEST);
nFuelPrices = numel( NATURAL_GAS_PARAMS() );
[buildingTypes, priceIndices] = getCaseIndices(nBuildTypes, nFuelPrices);

%% Assert and initilize variables

assert(numel(buildingTypes) == numel(priceIndices), ...
  'Incompatible number of building types and gas price indices.');

assert(numel(mixedAlphas) == numel(mixedSpikeAlphas), ...
  'Incompatible number of parameters for robust mixed algorithm.');

numScenarios = numel(buildingTypes);
numRobustLinfty = numel(LinftyAlphas);
numRobustMixed = numel(mixedAlphas);

% Prepare arrays to allow a single itertor to go over all
% algorithm-scenario pairs.
scenarioList = 1:numScenarios;
algorithmList = 1:(1*runBenchmark + 1*runNominal + numRobustLinfty*runRobustLinfty + numRobustMixed *runRobustMixed);
numAlgorithms = numel(algorithmList);

[scenarios, algorithms] = meshgrid(scenarioList, algorithmList);
numPairs = numel(scenarios);
scenarios = reshape(scenarios,[numPairs,1]);
algorithms = reshape(algorithms,[numPairs,1]);

%-1 = Benchmark. 0 = Nominal. 1 = Robust Linfty.  2 = Robust Mixed.
algorithmTypeKey = [-ones(runBenchmark,1);
                     zeros(runNominal,1);
                     ones(numRobustLinfty*runRobustLinfty,1);
                     2*ones(numRobustMixed *runRobustMixed,1)];
algorithmType = algorithmTypeKey(algorithms);

algorithmParameterKey = [1:runBenchmark,1:runNominal,1:numRobustLinfty*runRobustLinfty,1:numRobustMixed *runRobustMixed]';
algorithmParameters = algorithmParameterKey(algorithms);

% Prepare structure for output data
outputData(numPairs).Data.Power_Generation = [];
outputData(numPairs).Data.Heat_Generation = [];
outputData(numPairs).Data.Fuel_Consumption = [];
outputData(numPairs).Data.EstimatedCost = [];
outputData(numPairs).Data.TrueCost = [];
outputData(numPairs).Data.AlgorithmType = [];
outputData(numPairs).Data.AlgorithmParameters{1} = [];

outputData(numPairs).BuildingType = 0;
outputData(numPairs).PriceIndex = 0;
%% Run Algorithms on Scenarios

parfor iter = 1:numel(scenarios)
  jScenario = scenarios(iter);
  iAlgorithm = algorithms(iter);
  AlgType = algorithmType(iter);
  AlgParam = algorithmParameters(iter);
  
  building = buildingTypes(jScenario);
  priceInd = priceIndices(jScenario);
  
  switch AlgType
    case -1
      out = runBenchmarkAlgorithm('PriceIndex',priceInd,...
        'BuildingType', building);
    case 0
      out = runNominalAlgorithm('PriceIndex',priceInd,...
        'BuildingType', building);
    case 1
      out = runRobustLinftyAlgorithm('PriceIndex',priceInd,...
        'BuildingType', building,...
        'alpha', LinftyAlphas(AlgParam));
    case 2
      out = runRobustMixedAlgorithm('PriceIndex',priceInd,...
        'BuildingType',building,...
        'alphaMixed', mixedAlphas(AlgParam),...
        'alphaSpikes', mixedSpikeAlphas(AlgParam));
  end
  
  outputData(iter).Data = out;
  outputData(iter).BuildingType = building;  
  outputData(iter).PriceIndex = priceInd;
end

%% Parse to a new form, which is easier to analyze
% Initialize Structure - the values have no meaning, only the data
% structure.
currentScenarioData(numAlgorithms) = outputData(1).Data;
allScenariosData(numScenarios).Data = currentScenarioData;
allScenariosData(numScenarios).BuildingType = 0;
allScenariosData(numScenarios).PriceIndex = 0;

for jScenario = 1:numScenarios
  building = buildingTypes(jScenario);
  priceInd = priceIndices(jScenario);
  k = 1;
  for iter = 1:numPairs
    if (building == outputData(iter).BuildingType) && ...
       (priceInd == outputData(iter).PriceIndex)
      currentScenarioData(k) = outputData(iter).Data;
      k = k+1;
    end
  end
  allScenariosData(jScenario).Data = currentScenarioData;
  allScenariosData(jScenario).BuildingType = building;
  allScenariosData(jScenario).PriceIndex = priceInd;
end
end

function [buildIdx, priceIdx] = getCaseIndices(nBuildings, nFuelPrices)
buildIdx = repmat(1:nBuildings, 1, nFuelPrices);
priceIdx = repelem(1:nFuelPrices, 1, nBuildings);
end
