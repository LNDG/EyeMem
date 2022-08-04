function [outputArg1,outputArg2] = EM_plot_eyedata(timelock)
%plot pupil and various other eye-related data
%   Detailed explanation goes here

%% plot pupil time courses
baselinecorrect = 0;
usediff = 1;

age_colors = {'r' 'b'};

f=figure; hold on
clear sh
for iage = 1:2
  cfg=[];
  cfg.latency = [-2 7];
  cfg.channel = 'EYE_DIAMETER';
  timelock_plot = ft_selectdata(cfg, timelock{iage}); 
  if baselinecorrect
    cfg=[];
    cfg.latency = [-2 0];
    cfg.avgovertime = 'yes';
    timelockbaseline = ft_selectdata(cfg, timelock_plot)
    timelock_plot.individual = (timelock_plot.individual - timelockbaseline.individual) ...
      ./ timelockbaseline.individual *100;
  elseif usediff
    timelock_plot.individual = diff(timelock_plot.individual);
  end
  plotmean = squeeze(mean(timelock_plot.individual));
  plotsem = squeeze(std(timelock_plot.individual)) / sqrt(size(timelock_plot.individual,1));
    
  sh(iage) = shadedErrorBar(timelock_plot.time, plotmean, plotsem, {age_colors{iage}, 'Linewidth', 1.5}, 1);
end

%TODO stats
ax=gca


