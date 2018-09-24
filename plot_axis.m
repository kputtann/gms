function plot_axis
% This function sets the axes and line properties
% Puts grid and adjusts grid transparency
% Adjusts the LineWidth of the graph
% Adjusts the FontSize of axes

grid on;
h_line = findobj(gcf, 'type', 'line');
set(h_line, 'LineWidth',2);
h_axes = findobj(gcf, 'type', 'axes');
set(h_axes,'LineWidth',1,'FontSize',12,'GridAlpha',0.25);