function [] = steadyStateChecks()
clc;
close all;
load('..\Data\sv_mappings.mat', 'm_dot_f_map', 'power_map');

for k = 1:5
  figure(k)
  marginal_cost = zeros(9, 1);
  for j = 2:8
    marginal_cost(j) = 0.5 * (m_dot_f_map(k, j+1) - m_dot_f_map(k, j-1));
  end
  marginal_cost(1) = m_dot_f_map(k, 2) - m_dot_f_map(k, 1);
  marginal_cost(9) = m_dot_f_map(k, 9) - m_dot_f_map(k, 8);
  plot(power_map(k, :), marginal_cost);
  xlabel('Power generation');
  ylabel('Marginal m\_dot (proportional to marginal cost)')
  title(['Marginal cost when the bypass valve is in position ', num2str(k)]);
end