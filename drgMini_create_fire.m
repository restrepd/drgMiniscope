% Create a figure to ensure colormap functions are accessible
% This is used to get around the problem with the lack of definition
% of fire with Matlab 2025b




hFig = figure;

set(hFig, 'units','normalized','position',[.49 .1 .05 .3])


prain=[0:0.6/99:0.6];
pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
%             colormap jet

% Apply the 'fire' colormap to the current figure
colormap(fire);

% Retrieve the colormap RGB matrix
fireMap = colormap;

shading interp
ax=gca;
set(ax,'XTickLabel','')

% Save the colormap matrix to a .mat file
save('/Users/restrepd/Documents/GitHub/drgMiniscope/fireColormap2023b.mat', 'fireMap');

% Later, to reload and apply the saved colormap:
% load('fireColormap2023b.mat', 'fireMap');
% colormap(fireMap);