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
  properties (Access = public, Constant = true, Hidden = true)
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
  
  properties (GetAccess = public, SetAccess = private)
    scaling (1,1) double = NaN + 1i; % New maximum value for re-scaling the data, [kW]
    % The real part contains the new maximum value (NaN means default) and 
    % the imag part contains the scaling factor (1 by default)
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
        kwargs.isResidential (1,1) logical = ~contains(csvFullPath, "RefBldg");
        kwargs.newMaxPowerDemand (1,1) double = NaN; % [kW]
      end
      cdp.observationWindow = kwargs.windowSize;
      isRes = kwargs.isResidential;
      
      % Set import up:
      io = detectImportOptions(csvFullPath);
      if isRes
        io.SelectedVariableNames = io.SelectedVariableNames([1:2,5:14]);
      else
        io.SelectedVariableNames = io.SelectedVariableNames([1:4,6:8,10]);
      end
      io.PreserveVariableNames = true;      
      
      % Import:
      rawData = readtable(csvFullPath, io);
      
      % Process data:
      if isRes
        %{
        Power: Electricity:Facility + Electricity:HVAC + Fans:Electricity + 
               InteriorLights:Electricity + ExteriorLights:Electricity + 
               Appl:InteriorEquipment + Misc:InteriorEquipment - Cooling:Electricity            
        Heat:  Heating:Gas + HVACFan:Fans + Water Heater:WaterSystems + Cooling:Electricity / 0.7
        %}
        cdp.data = [sum(rawData{:,[2,6:11]},2) - rawData{:,4},... 
                    sum(rawData{:,[3,5,12]},2) + rawData{:,4} / CHP2004Provider.COP];
      else
        %{
        Power: Electricity:Facility + Fans:Electricity + 
               InteriorLights:Electricity + InteriorEquipment:Electricity
        Heat:  Gas:Facility - InteriorEquipment:Gas + Cooling:Electricity / 0.7
        %}
        cdp.data = [sum(rawData{:,[2,3,5,6]},2),... 
        diff(rawData{:,[8,7]},1,2) + rawData{:,4} / CHP2004Provider.COP];
      end

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
      
      % Rescaling
      cdp.rescalePower(kwargs.newMaxPowerDemand);
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
    
    function tshObj = fastForward(cdpObj, targetTime, inclusionFlag)
      % This function updates the observation window position of the current object and
      % optionally returns a TimestampedStatHolder (if an output was requested).
      arguments
        cdpObj (1,1) CHP2004Provider
        targetTime (1,1) datetime {mustBeNonempty}        
        % where targetTime will be with respect to the returned timestamps
        inclusionFlag (1,1) string {mustBeMember(inclusionFlag, ["beforeFirst", "first", "last", "afterLast"])} = "beforeFirst"
      end
      
      ts = cdpObj.timestamps;
      try
        switch inclusionFlag
          case "beforeFirst"
            cdpObj.currentWindowPosition = find(targetTime < ts, 1, 'first');
          case "first" % "right before", inclusive
            cdpObj.currentWindowPosition = find(targetTime <= ts, 1, 'first');
          case "last" 
            idx = find(targetTime >= ts, 1, 'last');
            cdpObj.currentWindowPosition = idx - cdpObj.observationWindow;
          case "afterLast"
            idx = find(targetTime > ts, 1, 'last');
            cdpObj.currentWindowPosition = idx - cdpObj.observationWindow;
        end
      catch me
        throw(addCause(me, MException('fastForward:InvalidWindowPosition',...
          ['Fast-forwarding failed because the resulting window position is not a positive value. '...
          'This can happen in the cases of "first" or "last" inclusionFlag (currently: "%s") when '...
          'the provided targetTime is outside the extreme timestamps of the available data.'], inclusionFlag)));
      end
      if nargout > 0 
        tshObj = cdpObj.next();
      end
    end
    
    function rescalePower(cdpObj, newMaxPowerDemand)
      % A function for rescaling the demands. 
      arguments
        cdpObj (1,1) CHP2004Provider
        newMaxPowerDemand (1,1) double {mustBeNonempty, mustBeReal} % new maximum value in [kW]
      end
      if cdpObj.currentWindowPosition ~= 1
        % Only allow if the object is rewinded. This is a failsafe to ensure that the user
        % doesn't change the scaling in the middle of a computation. 
        throw(MException('CHP2004Provider:rescaleUnrewindedObject',...
          'Please rewind/fast-forward the object to the first time step before changing the scaling!'));
      end
      if ~isnan(newMaxPowerDemand)        
        oldPowerScale = max(cdpObj.data(:,1));
        scalingFactor = newMaxPowerDemand ./ oldPowerScale;
        % Store the new maximum (real part) and the scaling factor (imag part) of the
        %   scaling property.
        cdpObj.scaling = newMaxPowerDemand + (scalingFactor*imag(cdpObj.scaling))*1i;
        % Rescale both the heat and the power demand by the same factor:
        cdpObj.data(:) = cdpObj.data .* scalingFactor;
      else
        if ~isnan(cdpObj.scaling)
          unscalingFactor = imag(cdpObj.scaling);
          cdpObj.data(:) = cdpObj.data ./ unscalingFactor;
          cdpObj.scaling = NaN + 1i;
        end
      end
    end
  end % methods          

end % classdef