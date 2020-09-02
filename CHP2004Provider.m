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
  t = chp.next();
  figure(); 
  yyaxis left; errorbar(t.hourOfVal, t.valMean(:,1), t.valStd(:,1),'DisplayName', 'Power'); 
  ylabel('Electric Power Consumption [kW]')
  yyaxis right; errorbar(t.hourOfVal, t.valMean(:,2), t.valStd(:,2),'DisplayName', 'Heat');
  ylabel('Heating & Cooling Consumption [kW]')
  xlim([0 24]); xticks(0:2:24); grid('on');
  title("Average hourly consumption in the range " + char(t.timeStart, "[dd/MM/yyyy") + char(t.timeEnd, "–dd/MM/yyyy)"));   
  
  chp_month = CHP2004Provider("E:\RefBldgHospitalNew2004_v1.3_7.1_4A_USA_MD_BALTIMORE.csv", ...
                        'windowSize', 2*CHP2004Provider.DEFAULT_WINDOW*CHP2004Provider.DATAPOINTS_PER_DAY);  
  %}
  properties (Access = public, Constant = true)
    DEFAULT_WINDOW = 14 % [days]
    DATAPOINTS_PER_DAY  = 24; % This is a property of the raw data.
  end
  
  properties (Access = private, Constant = true)
    YEAR = 2004;
    COP = 0.7; % Absorption chiller coefficient of performance
  end
  
  properties (GetAccess = public, SetAccess = immutable)
    observationWindow (1,1) double {mustBeInteger, mustBePositive} = CHP2004Provider.DEFAULT_WINDOW;
    nTimeSteps (1,1) double {mustBeInteger, mustBePositive} = realmax();
    nQuantities (1,1) double {mustBeInteger, mustBePositive} = 1;
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
      
      % Add the first midnight
      time0 = datetime(CHP2004Provider.YEAR, 1, 1, 0, 0, 0);
      data0 = cdp.data(end,:); % data is "periodic"
      
      % Add Feb 29th (Sunday) as an average of the previous and next and Sundays
      time29 = datetime(CHP2004Provider.YEAR, 2, 29, 1:24, 0, 0).';
      data29 = 0.5 * (cdp.data(1249:1272,:) + cdp.data(1561:1584,:));      
      
      % Build final time/data matrices
      cdp.timestamps = [time0; tmp(1:1416); time29; tmp(1417:end)];
      cdp.data = [data0; cdp.data(1:1416,:); data29; cdp.data(1417:end,:)];
      
      % Store some information about the data
      cdp.nTimeSteps = numel(cdp.timestamps);
      cdp.nQuantities = size(cdp.data, 2);
    end

    function tshObjNew = next(cdpObj)
      % Retrieve data:
      idx = cdpObj.currentWindowPosition : cdpObj.currentWindowPosition + cdpObj.observationWindow;
      tshObjNew = TimestampedStatHolder(cdpObj.timestamps(idx), cdpObj.data(idx,:),...
        'separateWeekends', true, 'hourlyStats', true, 'weights', ones(numel(idx),1), ...
        'biasCorrection', false, 'friSatWeekend', false, 'periodicOutput', true, ...
        'allowUnequalTimeIntervals', false);
      % NOTE/TODO: while it is more efficient to update the statistic values instead of
      % recomputing them (e.g., updating the mean involves subtracting from the known mean
      % the values no longer needed multiplied by their weights and adding the new values
      % multiplied by the new weights - this allows to save computations, especially with 
      % means spanning many values), it requires storing additional information in the 
      % TimestampedStatHolder class.
      
      % Advance cursor by a day:
      cdpObj.currentWindowPosition = cdpObj.currentWindowPosition + CHP2004Provider.DATAPOINTS_PER_DAY;
    end
    
    function tf = hasNext(cdpObj)
      % This function determines whether a call to next() will succeed. 
      % Can be used as a condition in a while loop that calls next() repeatedly.
      % Alternatively, one may perform a test such as:
      %   assert( hasNext(cdpObj), 'next:noMoreData', 'End of dataset reached!')
      tf = cdpObj.currentWindowPosition + cdpObj.observationWindow <= cdpObj.nTimeSteps;
    end
  end          
  
end % classdef