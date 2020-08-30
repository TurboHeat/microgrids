function [power_return, heat_return] = createDemandProfileVector(buildingtype, day, dt)
%% Function that creates the demand vectors for power and heat in Wh
%{
Buildingtypes:
          1. Large Hotel
          2. Full Service Restaurant
          3. Small Hotel
          4. Residential High
          5. Hospital
Days:     1. - 10.th January
          2. - 10.th April
          3. - 10.th July
          4. - 10.th October
Stepping of vector is defined with dt

Example of converting .csv data to what is found in the .mat files for the "Hospital" case:
1) Open RefBldgHospitalNew2004_v1.3_7.1_4A_USA_MD_BALTIMORE.csv
2) power = L218:L241; heat = M218:M241 + 0.7*N218:N241;
3) (repeat stage 2 for the three other example days)
%}
n_lines = 3600 / dt;
n_residential = 20;
power = [];
heat = [];
switch buildingtype
  case 1 % Large Hotel
    load('../Data/Demand Profiles/DemandLargeHotel.mat')
    power_return(1, :) = power(:, day);
    heat_return(1, :) = heat(:, day);
  case 2 % Full Service Restaurant
    load('../Data/Demand Profiles/DemandFSR.mat')
    power_return(1, :) = power(:, day);
    heat_return(1, :) = heat(:, day);
  case 3 % Small Hotel
    load('../Data/Demand Profiles/DemandSmallHotel.mat')
    power_return(1, :) = power(:, day);
    heat_return(1, :) = heat(:, day);
  case 4 % Residential One Unit * n_residential
    load('../Data/Demand Profiles/DemandResidentialOneUnit.mat')
    power_return(1, :) = power(:, day) * n_residential;
    heat_return(1, :) = heat(:, day) * n_residential;
  case 5 % Hospital
    load('../Data/Demand Profiles/DemandHospital.mat', 'power', 'heat')
    power_return(1, :) = power(:, day);
    heat_return(1, :) = heat(:, day);
end
for line = 2:n_lines
  % TODO: use repelem/repmat, or modify calling code to access elements using floor(t/dt) or interp1
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
