%% Choose Running Algorithms
runBenchmark = 1;

runNominal = 1;

runRobustLinfty = 1;
LinftyAlphas = [0.01:0.01:0.05, 0.1:0.05:0.7];

runRobustMixed = 1;
MixedAlphas = repmat([0.01:0.01:0.05,0.07,0.1,0.13,0.15,0.17,0.2],[1,10]);
MixedSpikeAlphas = [.5*sqrt(5760)*ones(1,11),...
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
BuildingTypes = [1,2,3,4,1,2,3,4,1,2,3,4];
PriceIndices  = [1,1,1,1,2,2,2,2,3,3,3,3];

%% Assert and initilize variables

assert(numel(BuildingTypes) == numel(PriceIndices), ...
       'Incompatible number of building types and gas price indices.');
   
assert(numel(MixedAlphas) == numel(MixedSpikeAlphas), ...
       'Incompatible number of parameters for robust mixed algorithm.');
   
numScenarios = numel(BuildingTypes);
numRobustLinfty = numel(LinftyAlphas);
numRobustMixed = numel(MixedAlphas);

% Prepare arrays to allow a single itertor to go over all
% algorithm-scenario pairs.
ScenarioList = 1:numScenarios;
AlgorithmList = 1:(1*runBenchmark + 1*runNominal + numRobustLinfty*runRobustLinfty + numRobustMixed *runRobustMixed);
numAlgorithms = numel(AlgorithmList);

[Scenarios,Algorithms] = meshgrid(ScenarioList,AlgorithmList);
numPairs = numel(Scenarios);
Scenarios = reshape(Scenarios,[numPairs,1]);
Algorithms = reshape(Algorithms,[numPairs,1]);

%-1 = Benchmark. 0 = Nominal. 1 = Robust Linfty.  2 = Robust Mixed.
AlgorithmTypeKey = [-ones(runBenchmark,1);zeros(runNominal,1);ones(numRobustLinfty*runRobustLinfty,1);2*ones(numRobustMixed *runRobustMixed,1)];
AlgorithmType = AlgorithmTypeKey(Algorithms);

AlgorithmParameterKey = [1:runBenchmark,1:runNominal,1:numRobustLinfty*runRobustLinfty,1:numRobustMixed *runRobustMixed]';
AlgorithmParameters = AlgorithmParameterKey(Algorithms);

% Prepare structure for output data
OutputData(numPairs).Data.Power_Generation = [];
OutputData(numPairs).Data.Heat_Generation = [];
OutputData(numPairs).Data.Fuel_Consumption = [];
OutputData(numPairs).Data.EstimatedCost = [];
OutputData(numPairs).Data.TrueCost = [];
OutputData(numPairs).Data.AlgorithmType = [];
OutputData(numPairs).Data.AlgorithmParameters{1} = [];

OutputData(numPairs).BuildingType = 0;
OutputData(numPairs).PriceIndex = 0;
%% Run Algorithms on Scenarios

parfor iter = 1:numel(Scenarios)
    jScenario = Scenarios(iter);
    iAlgorithm = Algorithms(iter);
    AlgType = AlgorithmType(iter);
    AlgParam = AlgorithmParameters(iter);
    
    Building = BuildingTypes(jScenario);
    PriceInd = PriceIndices(jScenario);
    
    switch AlgType
        case -1
            OutputData(iter).Data = runBenchmarkAlgorithm('PriceIndex',PriceInd,...
                'BuildingType',Building);
        case 0
            OutputData(iter).Data = runNominalAlgorithm('PriceIndex',PriceInd,...
                'BuildingType',Building);
        case 1
            OutputData(iter).Data = runRobustLinftyAlgorithm('PriceIndex',PriceInd,...
                'BuildingType',Building,...
                'alpha',LinftyAlphas(AlgParam));
        case 2
            OutputData(iter).Data = runRobustMixedAlgorithm('PriceIndex',PriceInd,...
                'BuildingType',Building,...
                'alphaMixed',MixedAlphas(AlgParam),...
                'alphaSpikes',MixedSpikeAlphas(AlgParam)); 
    end
    
    OutputData(iter).BuildingType = Building;
    OutputData(iter).PriceIndex = PriceInd;
end

%% Parse to a new form, which is easier to analyze
% Initialize Structure - the values have no meaning, only the data
% structure.
CurrentScenarioData(numAlgorithms) = OutputData(1).Data;
AllScenariosData(numScenarios).Data = CurrentScenarioData;
AllScenariosData(numScenarios).BuildingType = 0;
AllScenariosData(numScenarios).PriceIndex = 0;

for jScenario = 1:numScenarios
    Building = BuildingTypes(jScenario);
    PriceInd = PriceIndices(jScenario);
    k = 1;
    for iter = 1:numPairs
        if(Building == OutputData(iter).BuildingType && ...
           PriceInd == OutputData(iter).PriceIndex)
       CurrentScenarioData(k) = OutputData(iter).Data;
       k = k+1;
        end
    end
    AllScenariosData(jScenario).Data = CurrentScenarioData;
    AllScenariosData(jScenario).BuildingType = Building;
    AllScenariosData(jScenario).PriceIndex = PriceInd;
end
