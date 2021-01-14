classdef EmptyDataProvider < ConsumptionDataProvider
  %EMPTYDATAPROVIDER Summary of this class goes here
  %   Detailed explanation goes here

  methods
    function tshObjNew = next(~)
      tshObjNew = TimestampedStatHolder.empty();
    end
    
    function tf = hasNext(cdpObj)
      tf = false;
    end
  end
end

