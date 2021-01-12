function [price_kg_f, price_kWh] = NATURAL_GAS_PARAMS()
% This function contains some computations using constants. There is no need to perform
% them each time, so instead the end result is returned.
%{
Qr = 49736500;   % [J/kg]
h_env = 3.015e5; % [J/kg]
h_100 = 3.9748e5;% [J/kg]
eta_HRU = 0.89;
eta_b = 0.98;
mair_mf = 17.2 * 1.2;
Ph_mf = eta_HRU * (mair_mf + 1) * ((Qr * eta_b + mair_mf * h_env) / (mair_mf + 1) - h_100) / 3.6e6; %kWh/kg
% Fuel frice in Dollars per Thousand Cubic Feet, from the US Department of Energy website (eia)
price_ft3 = [18.41, ... Residential price (Aug 20)
              9.12, ... Residential price (Feb 20)
              8.50, ... Commercial price (Aug 20)
              6.87, ... Commercial price (Feb 20)
              3.88, ... Industrial price (Dec 19)
              2.55];  % Industrial price (Dec 19)
price_ft3 = linspace(18.41, 2.55, 4); % We take 4 evenly-spaced prices
density_CH4 = 0.68; % kg/m^3
ft3_m3 = power(0.3048, 3);
price_m3 = price_ft3 ./ (1000 * ft3_m3); %price in $/m^3
price_kg_f = price_m3 / density_CH4; % for MGT costs
price_kWh = price_kg_f / Ph_mf; % price in $/kWh, for heat tariff
%}
price_kg_f = [0.956092668150889,0.681538445502454,0.406984222854018,0.132430000205582];
price_kWh = [0.083413912137329,0.059460541749894,0.035507171362459,0.011553800975024];
end