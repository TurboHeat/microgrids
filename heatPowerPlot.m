function heatPowerPlot(heat_map)

valve_pos = 0:0.2:0.8;
np = size(heat_map,2);
[x, speed_level] = deal(zeros(1,np));
z = heat_map(1:end-1);

for k = 1:np - 1
  valve_pos_level = floor((k -1)/9) + 1;
  x(k) = valve_pos(valve_pos_level);
  speed_level(k) = mod(k-1, 9) + 1;
end

plot3(speed_level, x, z, '.'), xlabel('speed level'), ylabel('valve pos'), zlabel('heat. power')