function plot_legend(hL,hObj)
% This function sets the legend properties
% Adjusts the LineWidth and FontSize

set(hL,'FontSize',12);
hTL=findobj(hObj,'type','line'); % get the lines, not text
set(hTL,'LineWidth',2) % set linewidth
hTL=findobj(hObj,'type','Text'); % get the text
set(hTL,'FontSize',12) % set fontsize
% ax = gca; ax.LineWidth = 1; ax.FontSize = 14; ax.GridAlpha = 0.25; lgd = legend('So','To'); lgd.FontSize = 14;
% hL = gca; hObj = gca;
