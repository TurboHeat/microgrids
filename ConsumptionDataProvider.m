classdef (Abstract) ConsumptionDataProvider < handle & matlab.mixin.Heterogeneous
  % An interface for various loader classes.
  
  properties
    data (:,:) double % the actual data being provided
    timestamps (:,1) datetime % timestamps for the provided data
  end
  
  methods (Access = public, Abstract = true)
    % The NEXT method provides a StatHolder object that is "next after" tshObjOld in some sense.
    % It is up to the implementation to test whether a "next" output is possible, and act
    % accordingly if not.
    tshObjNew = next(cdpObj)
    
    % The HASNEXT method returns a logical indicating whether a "next" entry exists.
    tf = hasNext(cdpObj)
  end
  
  methods (Static, Sealed, Access = protected)
    function do = getDefaultScalarElement
      do = EmptyDataProvider();
    end
  end
  
end