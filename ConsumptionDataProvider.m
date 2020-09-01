classdef (Abstract) ConsumptionDataProvider < handle
  % An interface for various loader classes.
  
  properties
    data(:,:) double % the actual data being provided
    timestamps (:,1) datetime % timestemps for the provided data    
  end
  
  methods (Access = public, Abstract = true)
    % The NEXT method provides a StatHolder object that is "next after" tshObjOld in some sense.
    % Whenever the method is called, lastCdp is updated.
    tshObjNew = next(cdpObj)
  end

end