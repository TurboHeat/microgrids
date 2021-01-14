function [] = runAlgorithms(fuelIdx, minId, maxId, W, psf)
arguments
  fuelIdx (1,:) double = []
  minId (1,1) double = -Inf
  maxId (1,1) double =  Inf
  W (1,1) double = 10
  psf (1,1) = NaN % power scaling factor (of the consumer's demand, relative to the turbine's capacity)
end
%% Constants
OUTPUT_FOLDER = "../Data/Results";
if ~isfolder(OUTPUT_FOLDER), mkdir(OUTPUT_FOLDER); end

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
nFuelPrices = numel(NATURAL_GAS_PARAMS());
[buildingTypes, priceIndices] = getCaseIndices(nBuildTypes, nFuelPrices);

%% Reporting on the current run:
if nargin > 0 && ~isempty(fuelIdx)
  tsdisp("Only fuelIdx = " + mat2str(fuelIdx) + " will be considered!")
  keepIdx = priceIndices == fuelIdx;
  priceIndices = priceIndices(keepIdx);
  buildingTypes = buildingTypes(keepIdx);
else
  tsdisp("All " + nFuelPrices + " fuelIdx will be considered!")
end

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
scenarios = scenarios(:);
algorithms = algorithms(:);

%-1 = Benchmark. 0 = Nominal. 1 = Robust Linfty.  2 = Robust Mixed.
algorithmTypeKey = [-ones(runBenchmark,1);
                     zeros(runNominal,1);
                     ones(numRobustLinfty*runRobustLinfty,1);
                     2*ones(numRobustMixed *runRobustMixed,1)];
algorithmType = algorithmTypeKey(algorithms);

algorithmParameterKey = [1:runBenchmark,1:runNominal,1:numRobustLinfty*runRobustLinfty,1:numRobustMixed *runRobustMixed]';
algorithmParameters = algorithmParameterKey(algorithms);

% Prepare structure for output data
%{
outputData(numPairs).Data.Power_Generation = [];
outputData(numPairs).Data.Heat_Generation = [];
outputData(numPairs).Data.Fuel_Consumption = [];
outputData(numPairs).Data.EstimatedCost = [];
outputData(numPairs).Data.TrueCost = [];
outputData(numPairs).Data.AlgorithmType = [];
outputData(numPairs).Data.AlgorithmParameters{1} = [];

outputData(numPairs).BuildingType = 0;
outputData(numPairs).PriceIndex = 0;
%}
%% Run Algorithms on Scenarios
% Initialize parallel pool:
ts = datetime('now');
if isempty(gcp('nocreate'))
  try
    results = string(strsplit(evalc("parpool(" + W + ");"), '\n'));
    tsdisp(results(1), ts);
    tsdisp(results(2));
  catch ex
    fwrite(2, "Unable to initialize computation for fuelIdx=" + fuelIdx + ...
      ", id=[" + minId + ".." + maxId + "]. Additional information:" + string(char([13 10])));
    fwrite(2, ex.getReport());
    return
  end
end
% Suppress some warnings:
spmd % Issue commands to all workers in pool (Single Program, Multiple Data)
  warning('off', 'MATLAB:mir_warning_maybe_uninitialized_temporary');
end

nSc = numel(scenarios);
if batchStartupOptionUsed() || isempty(gcp('nocreate'))
  % In case of batch execution (running on a cluster)
  parfor iter = max(1, minId):min(nSc, maxId)
    % Unpack scenario configurations:
    jScenario = scenarios(iter);
    iAlgorithm = algorithms(iter);
    algType = algorithmType(iter);
    algParam = algorithmParameters(iter);

    building = buildingTypes(jScenario);
    priceInd = priceIndices(jScenario);
    
    % Determine if we need to compute the current file:
    fName = makeFilename(iter, algType, algParam, building, priceInd, psf);
    fPath = fullfile(OUTPUT_FOLDER, fName);
    if isfile(fPath)
      tsdisp(fName + " already exists - skipping.");
      continue
    else      
      tsdisp("Starting work on: " + fName);
      % Create an empty file (so that other workers will skip it):
      emptyFile(fPath); % In most cases, does nothing in parfor.
    end
    
    % Perform the computation
    tv = tic();
    switch algType
      case -1
        out = runBenchmarkAlgorithm('PriceIndex',priceInd,...
          'BuildingType', building,...
          'powerScalingFactor', psf);
      case 0
        out = runNominalAlgorithm('PriceIndex',priceInd,...
          'BuildingType', building,...
          'powerScalingFactor', psf);
      case 1
        out = runRobustLinftyAlgorithm('PriceIndex',priceInd,...
          'BuildingType', building,...
          'alpha', LinftyAlphas(algParam),...
          'powerScalingFactor', psf);
      case 2
        out = runRobustMixedAlgorithm('PriceIndex',priceInd,...
          'BuildingType',building,...
          'alphaMixed', mixedAlphas(algParam),...
          'alphaSpikes', mixedSpikeAlphas(algParam),...
          'powerScalingFactor', psf);
    end
    
    parsave(fPath, out, iter, algType, algParam, building, priceInd, psf);
    %{
    outputData(iter).Data = out;
    outputData(iter).BuildingType = building;  
    outputData(iter).PriceIndex = priceInd;
    %}
    tsdisp("Iteration #" + iter + " took " + round(toc(tv)) + "sec.");
  end
else  
  ppm = ParforProgressbar(nSc);
  parfor iter = 1:nSc
    % Unpack scenario configurations:
    jScenario = scenarios(iter);
    iAlgorithm = algorithms(iter);
    algType = algorithmType(iter);
    algParam = algorithmParameters(iter);

    building = buildingTypes(jScenario);
    priceInd = priceIndices(jScenario);
    
    % Determine if we need to compute the current file:
    fName = makeFilename(iter, algType, algParam, building, priceInd);
    fPath = fullfile(OUTPUT_FOLDER, fName);
    if isfile(fPath)
      tsdisp(fName + " already exists - skipping.");
      continue
    else
      tsdisp("Starting work on: " + fName);
    end
    
    % Perform the computation
    tv = tic();
    switch algType
      case -1
        out = runBenchmarkAlgorithm('PriceIndex',priceInd,...
          'BuildingType', building,...
          'powerScalingFactor', psf);
      case 0
        out = runNominalAlgorithm('PriceIndex',priceInd,...
          'BuildingType', building,...
          'powerScalingFactor', psf);
      case 1
        out = runRobustLinftyAlgorithm('PriceIndex',priceInd,...
          'BuildingType', building,...
          'alpha', LinftyAlphas(algParam),...
          'powerScalingFactor', psf);
      case 2
        out = runRobustMixedAlgorithm('PriceIndex',priceInd,...
          'BuildingType',building,...
          'alphaMixed', mixedAlphas(algParam),...
          'alphaSpikes', mixedSpikeAlphas(algParam),...
          'powerScalingFactor', psf);
    end
    fName = makeFilename(iter, algType, algParam, building, priceInd);
    parsave(fullfile(OUTPUT_FOLDER, fName), out, iter, algType, algParam, building, priceInd);
    %{
    outputData(iter).Data = out;
    outputData(iter).BuildingType = building;  
    outputData(iter).PriceIndex = priceInd;
    %}
    disp("Iteration #" + iter + " took " + toc(tv));
    ppm.increment(); %#ok<PFBNS>
  end
  delete(ppm);  
end
end

function emptyFile(fPath)
% Method 1:
fclose(fopen(fPath, 'w'));
% Method 2:
% parsave(fPath, [], [], [], [], [], []);
end

function parsave(fName, outputData, iter, algType, algParam, buildingId, fuelPriceId, psf)  
  save(fName, 'outputData', 'iter', 'algType', 'algParam', 'buildingId', 'fuelPriceId', 'psf');
end

function [buildIdx, priceIdx] = getCaseIndices(nBuildings, nFuelPrices)
buildIdx = repmat(1:nBuildings, 1, nFuelPrices);
priceIdx = repelem(1:nFuelPrices, 1, nBuildings);
end

function [nameStr] = makeFilename(iter, algType, algParams, bldgTypeId, fuelPriceId, powerScalingFactor)
% Returns the name of the file into which to save intermediate results.
% NOTE: The formatspec may need to be modified if the simulation grid is modified
%       significantly (e.g., more than 1E4-1 iterations).
nameStr = compose("I%04u_AT%02d_AP%03u_B%1u_F%1u_PSF%4.2f.mat", ...
  iter, algType, algParams, bldgTypeId, fuelPriceId, powerScalingFactor);
end

function [] = tsdisp(msg, ts)
arguments
  msg (1,1) string
  ts (1,1) datetime = datetime('now')
end
  disp(datestr(ts,"[dd-mm-yyyy HH:MM:SS.FFF] ") + msg);
end
