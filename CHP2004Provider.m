classdef CHP2004Provider < ConsumptionDataProvider
  % A provider class for working with the "Commercial and Residential Hourly Load Profiles 
  % for all TMY3 Locations in the United States" dataset found in:
  %   https://openei.org/datasets/dataset/commercial-and-residential-hourly-load-profiles-for-all-tmy3-locations-in-the-united-states
  %   https://openei.org/datasets/files/961/pub/COMMERCIAL_LOAD_DATA_E_PLUS_OUTPUT/USA_NY_New.York-Central.Park.725033_TMY3/
  %
  % Calling example:
  %{
  chp = CHP2004Provider("E:\RefBldgHospitalNew2004_v1.3_7.1_4A_USA_MD_BALTIMORE.csv", ...
                        'windowSize', CHP2004Provider.DEFAULT_WINDOW*CHP2004Provider.DATAPOINTS_PER_DAY);
  %}
  properties (Access = public, Constant = true)
    DEFAULT_WINDOW = 14 % [days]
    DATAPOINTS_PER_DAY  = 25; % midnight on both sides is used
  end
  
  properties (Access = private, Constant = true)
    YEAR = 2004;
    COP = 0.7; % Absorption chiller coefficient of performance
  end
  
  properties (GetAccess = public, SetAccess = immutable)
    observationWindow (1,1) double {mustBeInteger, mustBePositive} = CHP2004Provider.DEFAULT_WINDOW;
  end
  
  properties (Access = private)
    % Stateful cursor for obtaining elements:
    currentWindowPosition (1,1) double {mustBeInteger, mustBePositive} = 1;
  end
  
  methods 
    function cdp = CHP2004Provider(csvFullPath, kwargs)
      arguments
        csvFullPath (1,1) string {mustBeNonempty}
        % Window size, in dataset-specific time units, for computation of statistics:
        kwargs.windowSize (1,1) double {mustBeInteger, mustBePositive} = CHP2004Provider.DEFAULT_WINDOW
      end
      cdp.observationWindow = kwargs.windowSize;
      % Set import up:
      io = detectImportOptions(csvFullPath);
      io.SelectedVariableNames = io.SelectedVariableNames([1,12:14]);
      io.PreserveVariableNames = true;
      % Import:
      rawData = readtable(csvFullPath, io);
      % Process data:
      cdp.data = [rawData{:,2}, rawData{:,3} + CHP2004Provider.COP * rawData{:,4}];
      
      % Fix midnights (because 24:00 is considered ambiguous by MATLAB)
      tmp = datetime(strrep(CHP2004Provider.YEAR + "/" + rawData{:,1}, "24:00", "23:59"),...
                     'InputFormat', 'yyyy/MM/dd  HH:mm:ss');
      isMidnight = tmp.Minute == 59;
      tmp(isMidnight) = tmp(isMidnight) + duration([0 1 0]); % add a minute
      % Add Feb 29th (Sunday) as an average of the previous and next and Sundays
      time29 = datetime(CHP2004Provider.YEAR, 2, 29, 1:24, 0, 0).';
      data29 = 0.5 * (cdp.data(1249:1272,:) + cdp.data(1561:1584,:));      
      cdp.timestamps = [tmp(1:1416); time29; tmp(1417:end)];
      % Build final time/data matrices
    end

    % TODO: test OR try/catch OR circshift to make sure we're not trying to access invalid indices
    function tshObjNew = next(cdpObj)
      % Retrieve data:
      idx = cdpObj.currentWindowPosition : cdpObj.currentWindowPosition + cdpObj.observationWindow;
      tshObjNew = TimestampedStatHolder(cdpObj.timestamps(idx), cdpObj.data(idx,:));
      % Advance cursor:
      % NOTE/TODO: while it is more efficient to update the statistic values instead of
      % recomputing them (e.g., updating the mean involves subtracting from the known mean
      % the values no longer needed multiplied by their weights and adding the new values
      % multiplied by the new weights - this allows to save computations, especially with 
      % means spanning many values), it requires storing additional information in the 
      % TimestampedStatHolder class.
      
      cdpObj.currentWindowPosition = cdpObj.currentWindowPosition + CHP2004Provider.DATAPOINTS_PER_DAY;
    end

  end          
  
end % classdef