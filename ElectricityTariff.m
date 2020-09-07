classdef ElectricityTariff < handle
  % A class for storing hourly electric tariff infomation.
  properties (GetAccess = public, SetAccess = immutable)
    tariff (:,1) double    % [$/kW]
    timeOfDay (:,1) double % [hour] corresponding to tariff data
    name (1,1) string
  end
  
  methods
    function eto = ElectricityTariff(tariff, name, timeOfDay)
      arguments
        tariff (:,1) double 
        name (1,1) = ""
        timeOfDay (:,1) double = (0:23).' % [hour] corresponding to tariff data
      end      
      eto.tariff = tariff;
      eto.name = name;
      eto.timeOfDay = timeOfDay;
    end
    
    function tariff = tariffAtTime(eto, queryTime)
      %% 0-order interpolation:
      tariff = eto.tariff( floor(queryTime) + 1 );
%       tariff = reshape(nakeinterp1(eto.timeOfDay, eto.tariff, floor(queryTime(:))), size(queryTime));

      %% 1st order interpolation:
%     tariff = reshape(nakeinterp1(eto.timeOfDay, eto.tariff, queryTime(:)), size(queryTime));
    end
  end   
  
end