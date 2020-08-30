classdef TimestampedStatHolder < handle
  % This is a data class for holding timestamped statistics of some provided data.
  % This class stores only a summary of the data and not individual values.
  
  properties (GetAccess = public, SetAccess = immutable)
    % General
    nValues
    
    % Time
    dt
    timeStart
    timeEnd
    
    % Statistics
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
      end
      assert( numel(timeVec) == size(values,1), 'TimestampedStatHolder:incorrectInputLengths',...
        'The inputs must have the same number of rows!');
      assert( issorted(timeVec, 'strictascend'), 'TimestampedStatHolder:unsortedTimestamps',...
        'Timestamps must be unique and sorted in increasing order!');
      dt = diff(timeVec); % this results in an array of durations
      if kwargs.allowUnequalTimeIntervals
        tshObj.dt = dt;
      else
        assert( numel(unique(dt)) == 1, 'TimestampedStatHolder:unequalTimeSteps',...
          'Unequal time intervals detected!');
        tshObj.dt = dt(1);
      end
      
      %% Object creation:
      tshObj.nValues = numel(timeVec);
      tshObj.timeStart = timeVec(1);
      tshObj.timeEnd = timeVec(end);
      
      if kwargs.separateWeekends
        we = TimestampedStatHolder.isweekend(timeVec, kwargs.friSatWeekend);
        % Mean:
        tshObj.valMean = [TimestampedStatHolder.weightedMean(values( we, :), kwargs.weights( we));
                          TimestampedStatHolder.weightedMean(values(~we, :), kwargs.weights(~we))];
        % Standard deviation:
        tshObj.valStd = [TimestampedStatHolder.weightedStdev(values( we, :), kwargs.weights( we), tshObj.valMean(1), kwargs.biasCorrection);
                         TimestampedStatHolder.weightedStdev(values(~we, :), kwargs.weights(~we), tshObj.valMean(2), kwargs.biasCorrection)];
      else
        % Mean:
        tshObj.valMean = TimestampedStatHolder.weightedMean(values, :, kwargs.weights);
        % Standard deviation:
        tshObj.valStd = TimestampedStatHolder.weightedStdev(values, :, kwargs.weights, tshObj.valMean, kwargs.biasCorrection);
      end
    end % constructor
        
  end % public methods
  
  methods (Access = protected, Static = true)    
    function valMean = weightedMean(values, weights)
      valMean = sum(weights .* values, 1) ./ sum(weights);
    end
    
    function valStd = weightedStdev(values, weights, wMean, biasCorrection)
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
        valStd = sqrt( sum(weights .* (values - wMean).^2, 1) / ...
          ( (n-1)/n * sum(weights)) );
      else
        % No bias correction:
        valStd = sqrt( sum(weights .* (values - wMean).^2, 1) / ...
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