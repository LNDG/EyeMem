load('/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/preproc/eye/eye_sub-13.mat')

smoothit = 1;

cfg=[];
cfg.latency = [0 5];
data = ft_selectdata(cfg, data)

cfg=[];
cfg.keeptrials='yes';
timelock = ft_timelockanalysis(cfg, data)

%     x_edges = np.linspace(192, 832, num=456) # 455+1 since shape of hist otherwise 454*340
%     y_edges = np.linspace(144, 624, num=342)
x_edges = linspace(192, 832, 456); % 455+1 since shape of hist otherwise 454*340
y_edges = linspace(144, 624, 342);

histdat = nan( length(x_edges)-1, length(y_edges)-1, 150);

for itrial = 1:150
  X = squeeze(timelock.trial(itrial, 2, :));
  Y = squeeze(timelock.trial(itrial, 3, :));
  if smoothit
    Image = histcounts2(X, Y, x_edges, y_edges);
    % smoothing
    sigma = 10; % set sigma to the value you need
    sz = 2*ceil(2.6 * sigma) + 1; % See note below
    mask = fspecial('gauss', sz, sigma);
    histdat(:,:,itrial) = conv2(Image, mask, 'same');
  else
    histdat(:,:,itrial) = histcounts2(X, Y, x_edges, y_edges);
  end
end

close all
figure; imagesc(histdat(:,:,1)); colorbar
figure; imagesc(mean(histdat,3)); colorbar

%% pca: make 150 by X matrix
pcadat = permute(histdat, [3 1 2]);
pcadat = reshape(pcadat, 150, []);
[coeff, score, latent, ~, explained] = pca(pcadat);
explained = (explained ./ sum(explained)) *100;

nx = size(histdat,1);
ny = size(histdat,2);
close all
figure;
for icomp=1:9
  comp = reshape(coeff(:,icomp), [nx ny] );
  subplot(3,3,icomp); 
  maxabsval = max(abs(comp(:)));
  imagesc(1:nx, 1:ny, comp, [-maxabsval maxabsval]); colorbar
  title(explained(icomp))
end

