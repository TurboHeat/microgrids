classdef TimestampedStatHolder < handle
  % This is a data class for holding timestamped statistics of some provided data.
  % This class stores only a summary of the data and not individual values.
  
  properties (GetAccess = public, SetAccess = immutable)
    %% General
    nValues
    
    %% Time
    dt
    timeStart
    timeEnd
    hourOfVal = NaN % This property will be populated if hourly data was requested
    
    %% Statistics
    % These will be either 1xN or 2xN depending on whether weekends are treated
    % separately, where N is the number of quantities being processed (typically 1 or 2,
    % for power and/or heat)
    valMean
    valStd
    % etc ...
  end % properties
  
  methods
    function tshObj = TimestampedStatHolder(timeVec, values, kwargs)
      %% Input testing:
      arguments
        timeVec (:,1) datetime {mustBeNonempty}
        values (:,:) double {mustBeNonempty, mustBeNonNan}
        kwargs.separateWeekends (1,1) logical = true
        kwargs.weights (:,1) double {mustBeNonnegative, mustBeFinite} = ones(numel(timeVec),1)
        kwargs.biasCorrection (1,1) logical = false
        kwargs.friSatWeekend (1,1) logical = false
        kwargs.allowUnequalTimeIntervals (1,1) logical = false
        kwargs.hourlyStats (1,1) logical = true
        kwargs.periodicOutput (1,1) logical = false
      end
      assert( numel(timeVec) == size(values,1), ...
        'TimestampedStatHolder:incorrectInputLengths',...
        'The inputs must have the same number of rows!');
      assert( issorted(timeVec, 'strictascend'), ...
        'TimestampedStatHolder:unsortedTimestamps',...
        'Timestamps must be unique and sorted in increasing order!');
      dt = diff(timeVec); % this results in an array of durations
      if kwargs.allowUnequalTimeIntervals
        tshObj.dt = dt;
      else
        assert( numel(unique(dt)) == 1, ...
          'TimestampedStatHolder:unequalTimeSteps',...
          'Unequal time intervals detected!');
        tshObj.dt = dt(1);
      end
      if kwargs.hourlyStats
        assert( all(~timeVec.Minute) && all(~timeVec.Second),...
          'TimestampedStatHolder:unroundedTimestamps', ...
          'When hourly outputs are requested, input timestamps must correspond to round hours!');
      end
      
      %% Object creation:
      tshObj.nValues = numel(timeVec);
      tshObj.timeStart = timeVec(1);
      tshObj.timeEnd = timeVec(end);
      
      if kwargs.separateWeekends
        we = TimestampedStatHolder.isweekend(timeVec, kwargs.friSatWeekend);
        if kwargs.hourlyStats
          [tshObj.valMean, tshObj.valStd, tshObj.hourOfVal] = ...
            TimestampedStatHolder.hourlyWeightedMeanStd(timeVec(we), values(we, :), kwargs.weights(we));
        else          
          % Mean:
          tshObj.valMean = [TimestampedStatHolder.weightedMean(values( we, :), kwargs.weights( we));
                            TimestampedStatHolder.weightedMean(values(~we, :), kwargs.weights(~we))];
          % Standard deviation:
          tshObj.valStd = [TimestampedStatHolder.weightedStdev(values( we, :), kwargs.weights( we), tshObj.valMean(1), kwargs.biasCorrection);
                           TimestampedStatHolder.weightedStdev(values(~we, :), kwargs.weights(~we), tshObj.valMean(2), kwargs.biasCorrection)];
        end
      else
        if kwargs.hourlyStats
          [tshObj.valMean, tshObj.valStd, tshObj.hourOfVal] = ...
            TimestampedStatHolder.hourlyWeightedMeanStd(timeVec, values, kwargs.weights);          
        else
          % Mean:
          tshObj.valMean = TimestampedStatHolder.weightedMean(values, kwargs.weights);
          % Standard deviation:
          tshObj.valStd = TimestampedStatHolder.weightedStdev(values, kwargs.weights, tshObj.valMean, kwargs.biasCorrection);
        end
      end
      
      if kwargs.periodicOutput
        tshObj.valMean = tshObj.valMean([1:end,1],:);
        tshObj.valStd = tshObj.valStd([1:end,1],:);
        tshObj.hourOfVal(end+1) = 24;
      end
    end % constructor
        
  end % public methods
  
  methods (Access = protected, Static = true)
    function [g, uh] = findSameHourGroups(timestamps)
      h = timestamps.Hour;
      [g,uh] = findgroups(h);
    end
    
    function [meanValsVec, stdValsVec, uh] = hourlyWeightedMeanStd(timestamps, values, weights)
      % This function produces hourly weighted averages of the provided values
      % It is assumed that weekends have been included/excluded in advance     
      [g,uh] = TimestampedStatHolder.findSameHourGroups(timestamps);
      [meanValsVec,stdValsVec] = splitapply(@TimestampedStatHolder.wMeanAndStd, values, weights, g);      
    end    
    
    function [meanValsVec, h] = hourlyWeightedMean(timestamps, values, weights)
      % This function produces hourly weighted averages of the provided values
      % It is assumed that weekends have been included/excluded in advance    
      g = TimestampedStatHolder.findSameHourGroups(timestamps);
      meanValsVec = splitapply(@TimestampedStatHolder.weightedMean, values, weights, g);      
    end
    
    function [stdValsVec, h] = hourlyWeightedStd(timestamps, values, weights)
      % This function produces hourly weighted standard deviations of the provided values
      % It is assumed that weekends have been included/excluded in advance   
      g = TimestampedStatHolder.findSameHourGroups(timestamps);
      stdValsVec = splitapply(@TimestampedStatHolder.weightedStdev, values, weights, g);      
    end
    
    function [meanOfVals, stdOfVals] = wMeanAndStd(values, weights, biasCorrection)
      % This function produces both the weighted averages and standard deviations of the provided values
      % It is assumed that weekends have been included/excluded in advance      
      if nargin < 3
        biasCorrection = false;
      end
      meanOfVals = TimestampedStatHolder.weightedMean(values, weights);
      stdOfVals = TimestampedStatHolder.weightedStdev(values, weights, meanOfVals, biasCorrection);      
    end
    
    function meanOfVals = weightedMean(values, weights)
      meanOfVals = sum(weights .* values, 1) ./ sum(weights);
    end
    
    function stdOfVals = weightedStdev(values, weights, wMean, biasCorrection)
      % https://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/weightsd.pdf
      if nargin < 4
        biasCorrection = false;
      end
      if nargin < 3
        wMean = TimestampedStatHolder.weightedMean(values, weights);
      end
      if biasCorrection
        % Including bias correction:
        n = nnz(weights);
        stdOfVals = sqrt( sum(weights .* (values - wMean).^2, 1) / ...
          ( (n-1)/n * sum(weights)) );
      else
        % No bias correction:
        stdOfVals = sqrt( sum(weights .* (values - wMean).^2, 1) / ...
                      sum(weights));
      end      
    end
    
    function tf = isweekend(dtObj, isWeekendFriSat)
      %ISWEEKEND True for datetimes occurring on a weekend.
      % See also datetime.isweekend.
      arguments
        dtObj (:,:) datetime {mustBeNonempty}
        isWeekendFriSat (1,1) logical = false
      end
      dow = weekday(dtObj);
      if isWeekendFriSat
        tf = (dow == 6) | (dow == 7);
      else % Weekend is Sat-Sun
        tf = (dow == 1) | (dow == 7);
      end
    end
  end
  
end % Classdef