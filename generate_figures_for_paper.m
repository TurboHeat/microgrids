%--------------------------------------------------------------%
% File: generate_figures_for_paper.m (script)
% Author: Miguel Dias
% Date 14/03/17
% v1.0
% Description: Plots the figures 11-14 for the final paper version. Only
% plots the case study images and used matlabfrag to export it to eps
% file to yield readable images.
%--------------------------------------------------------------%

close all; clearvars; clc;
addpath(genpath('C:\Users\migueld\Dropbox\Technion Grand Energy Program\Miguel Dias\Data'));
load data_for_plotting.mat ;
savepath = '..\Data\Case_Study_Images_Final_Paper\eps_v2\';
end_time = 24; %end time in h
dt = 15; %in s
n_lines = 3600 / dt;
t_plot = (1:n_lines * end_time)' ./ n_lines;
save_fig = 1;

%% d_index; 1-3-> Large; 4-6->FSR; 7-9->Small restaurant; % 10-12->Residential
LH = 1:3; %Large Hotel (B1)
FSR = 4:6; %Full Service Restaurant (B2)
SH = 7:9; %Small Restaurant (B3)
RC = 10:12; %Residential Community (B4)

choice = 'RC'; % Choose which building to plot data

switch choice
  case 'LH' % Large Hotel
    building = LH;
    ylimPower = [0, 750];
    ylimHeat = [0, 1100];
    figname = 'FC1B1';
  case 'FSR'
    building = FSR;
    ylimPower = [0, 110];
    ylimHeat = [0, 140];
    figname = 'FC1B2';
  case 'SH'
    building = SH;
    ylimPower = [0, 220];
    ylimHeat = [0, 95];
    figname = 'FC1B3';
  case 'RC'
    building = RC;
    ylimPower = [0, 170];
    ylimHeat = [0, 410];
    figname = 'FC1B4';
  otherwise
    disp('Invalid choice!');
end

pos = [403, 246, 560, 420];
FS = 20;
%{
%% Code to save as 3 subfigures
% for j=building(1):building(end)
%     h1=figure(1); %Power
%     set(h1,'position', pos);
%     subplot(1,3,i);
%     plot(t_plot,power_demand(:,building(i))/1e3, t_plot, power_MGT(:,building(i))/1e3, t_plot, (power_demand(:,building(i))-power_MGT(:,building(i)))/1e3, 'LineWidth', 2); %in kW
%     xlim([0 25]); ylim(ylimPower); grid on;
%     xlabel('Time [h]', 'FontSize', FS); ylabel('Demand and commitment [kW]', 'FontSize', FS); hold off;
%
%     h2=figure(2); %Heat
%     set(h2,'position', pos);
%     subplot(1,3,i)
%     plot(t_plot,heat_demand(:,building(i))/1e3,t_plot, heat_MGT(:,building(i))/1e3, t_plot, (heat_demand(:,building(i))-heat_MGT(:,building(i)))/1e3, 'LineWidth', 2);   %in Kw
%     xlim([0 25]); ylim(ylimHeat); grid on;
%     xlabel('Time [h]', 'FontSize', FS); ylabel('Demand and commitment [kW]', 'FontSize', FS);
%
%     i=i+1;
% end
%
% if save_fig
%     set(gcf, 'PaperPositionMode', 'auto');
%     print(h1,[savepath,'Power_MGT', figname],'-depsc')
%     print(h2,[savepath, 'Heat_MGT', figname],'-depsc')
%     disp(['All data saved to folder ', savepath]);
% end
%}

%% Code to save all figures individually
i = 1;
for j = building(1):building(end)
  figure(i); %Power
  set(gcf, 'position', pos);
  set(gca, 'fontsize', 14);
  plot(t_plot, power_demand(:, building(i))/1e3, t_plot, power_MGT(:, building(i))/1e3, t_plot, (power_demand(:, building(i)) - power_MGT(:, building(i)))/1e3, 'LineWidth', 2); %in kW
  xlim([0, 25]);
  ylim(ylimPower);
  grid on; %title(['Power D', num2str(i), ' B ', num2str(j)]);
  xlabel('Time [h]', 'FontSize', FS);
  ylabel('Demand and commitment [kW]', 'FontSize', FS);
  hold off;
  
  figure(i+3); %Heat
  set(gcf, 'position', pos);
  set(gca, 'fontsize', 14);
  plot(t_plot, heat_demand(:, building(i))/1e3, t_plot, heat_MGT(:, building(i))/1e3, t_plot, (heat_demand(:, building(i)) - heat_MGT(:, building(i)))/1e3, 'LineWidth', 2); %in Kw
  xlim([0, 25]);
  ylim(ylimHeat);
  grid on; %title(['Heat D', num2str(i), ' B ', num2str(j)])
  xlabel('Time [h]', 'FontSize', FS);
  ylabel('Demand and commitment [kW]', 'FontSize', FS);
  if j == 11 | j == 12
    ylim([0, 200]);
  end
  i = i + 1;
end

if save_fig
  set(gcf, 'PaperPositionMode', 'auto');
  print('-f1', [savepath, 'Power_MGT', figname, 'D1_v2'], '-depsc')
  print('-f2', [savepath, 'Power_MGT', figname, 'D2_v2'], '-depsc')
  print('-f3', [savepath, 'Power_MGT', figname, 'D3_v2'], '-depsc')
  
  print('-f4', [savepath, 'Heat_MGT', figname, 'D1_v2'], '-depsc')
  print('-f5', [savepath, 'Heat_MGT', figname, 'D2_v2'], '-depsc')
  print('-f6', [savepath, 'Heat_MGT', figname, 'D3_v2'], '-depsc')
  
  disp(['All data saved to folder ', savepath]);
end
