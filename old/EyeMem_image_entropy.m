%% calculate image entropy for EyeMem stimuli

stimpath = '/Users/kloosterman/Dropbox/PROJECTS/EyeMem/fMRI_final/fMRI_exp_final/stimuli_640x480';  % /EYEMEM004/preproc2/fractals
cd(stimpath)
cond_names = {'fractals' 'landscapes' 'naturals1' 'streets1' 'streets2'};
ncond = length(cond_names);
npics = 60;

%intialize entropy matrix
im_entropy = zeros(ncond,npics); % dimord cond pic

%now loop
for icond = 1:ncond
    cd(fullfile(stimpath, cond_names{icond}))
    pic_list = dir('*.png');
    if isempty(pic_list)
        pic_list = dir('*.bmp');
    end
    
    for ipic = 1:length(pic_list)
        im_entropy(icond,ipic) = entropy(imread(pic_list(ipic).name));
    end
end

%% Plotting
close all
figure;
barweb(mean(im_entropy,2), std(im_entropy,0,2) / sqrt(npics), 0.5, 'Entropy');
hold on; box on
legend(cond_names); legend boxoff;
ax=gca;
ax.FontSize = 14;
ax.YLim = [6.5 8];
ax.YLabel.String = {'entropy'};