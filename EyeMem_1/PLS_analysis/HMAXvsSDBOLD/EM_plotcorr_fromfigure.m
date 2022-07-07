function EM_plotcorr_fromfigure(curfig)

if nargin == 0
  curfig = gcf;
end

ax=gca(curfig); % get axes of figure that is open

X = [ax.Children(1:end-1).XData]';
Y = [ax.Children(1:end-1).YData]';

flipit=0;
if flipit
  Y = Y*-1;
end

f = figure; 
f.Position =  [1211 894 150 150];

sc = scatter(X, Y, 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', 50, 'LineWidth', 1); % 'Linewidth', 1.5
box on
lsline
% xlabel('dprime')
xlabel('criterion')
ylabel('HMAX slope')
title(sprintf('r = %1.2f', corr(X, Y)))

% PREOUT = '/Users/kloosterman/OneDrive/Documents/Eyemem/Figures/rawpics';
PREOUT = '/Users/kloosterman/Dropbox/PROJECTS/EyeMem/plots';
saveas(f, fullfile(PREOUT, 'scatterfromPLS_critOA_alone'), 'epsc')


