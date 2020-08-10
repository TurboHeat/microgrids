
%% Function that creates the demand vectors for power and heat in Wh
% Buildingtypes:
%           1. Large Hotel
%           2. Full Service Restaurant
%           3. Small Hotel
%           4. Residential High
%           5. Hospital
% Days:     1. - 10.th January
%           2. - 10.th April
%           3. - 10.th July
%           4. - 10.th October
% Stepping of vector is defined with dt
function [power_return, heat_return] = createDemandProfileVector(buildingtype, day, dt)
n_lines = 3600 / dt;
n_residential = 20;
power = [];
heat = [];
switch buildingtype
  case 1 % Large Hotel
    load('Demand Profiles\DemandLargeHotel.mat')
    %load('DemandLargeHotel.mat')
    %         load('Demand Profiles/DemandLargeHotel.mat') % mac syntax
    power_return(1, :) = power(:, day);
    heat_return(1, :) = heat(:, day);
  case 2 % Full Service Restaurant
    load('Demand Profiles\DemandFSR.mat')
    %         load('Demand Profiles/DemandFSR.mat') % mac syntax
    power_return(1, :) = power(:, day);
    heat_return(1, :) = heat(:, day);
  case 3 % Small Hotel
    load('Demand Profiles\DemandSmallHotel.mat')
    %         load('Demand Profiles/DemandSmallHotel.mat') % mac syntax
    power_return(1, :) = power(:, day);
    heat_return(1, :) = heat(:, day);
  case 4 % Residential One Unit * n_residential
    load('Demand Profiles\DemandResidentialOneUnit.mat')
    %         load('Demand Profiles/DemandResidentialOneUnit.mat') % mac syntax
    power_return(1, :) = power(:, day) * n_residential;
    heat_return(1, :) = heat(:, day) * n_residential;
  case 5 % Hospital
    load('Demand Profiles\DemandHospital.mat')
    %         load('Demand Profiles/DemandHospital.mat') % mac syntax
    power_return(1, :) = power(:, day);
    heat_return(1, :) = heat(:, day);
end
for line = 2:n_lines
  power_return(line, :) = power_return(1, :); % copy values 240 times, i.e. each time step in an hour has the same demand
  heat_return(line, :) = heat_return(1, :);
end
% Bring vectors together in one long row vector. Merge the 24 column
% vector into a single array
power_return = power_return(:)';
heat_return = heat_return(:)';
% Make kW to W
power_return = power_return' * 1000;
heat_return = heat_return' * 1000;
end
