function [] = apply_axis_properties(ax, linewidth)
% Apply common properties to axis
%{
-------------------------- Description ------------------------------------
Apply common properties, such as font weight and box to axis object.

---------------------------- Input ----------------------------------------
ax : Axis object
linewidth : Line width for axis (0.8 seems thick enough)

---------------------------- Output ---------------------------------------
None

----------------------------- Vesrions ------------------------------------
%{
v1 : Suraj R Pawar, 7-21-2020
    - Initialize
%}
%}
box on;
ax.FontWeight = 'bold';
ax.LineWidth = linewidth;
end

