function Eyemem_convert_edf_prepare( )
% Convert raw data and make file structure

% overwrite = 0;

if ismac
    root_folder = '/Users/kloosterman/gridmaster2012';
else
    root_folder = '/home/mpib';
end
PREIN = fullfile(root_folder, 'LNDG/EyeMem/eye_raw' );
PREOUT = fullfile(root_folder, 'LNDG/EyeMem/data' );

cd(PREIN)
subj_list = dir(fullfile(PREIN, 'EYEMEM*'));

ijob = 0; inputfiles={};
for isub = 1:length(subj_list)
    subj_id = subj_list(isub).name;
    
    cd(fullfile(PREIN, subj_list(isub).name))
    
    subj_no = str2num(subj_list(isub).name(7:9));
    runlist = dir(sprintf('S%dp*.edf', subj_no));
    for irun = 1:length(runlist)
        ijob = ijob + 1;
        edf_name = runlist(irun).name;
        
        inputfiles{ijob} = fullfile(pwd, edf_name);
        
        exp_phase = edf_name(strfind(edf_name, 'p') + 1); % index of exp phase: 1, study; 3, resting state
        if strcmp(exp_phase, '1')
            run_no = edf_name(strfind(edf_name, 'r') + 1); % index of exp phase: 1, study; 3, resting state
            folder_name = ['run' run_no];
            output_filename = sprintf('%s_run%s', subj_id, run_no );
        elseif strcmp(exp_phase, '3')
            folder_name = 'restingstate';
            output_filename = sprintf('%s_restingstate', subj_id);            
        end
        outputfiles{ijob} = fullfile(PREOUT, subj_list(isub).name, 'eye', folder_name, output_filename );
       
    end

end
cd(PREIN)

cfg = [];
cfg.function_to_run = 'Eyemem_convert_edf_core'; %         [asc, data, event] = EM_eyelink_to_fieldtrip(edf_name);
cfg.compile = 'yes';
cfg.parallel = 'torque';
% cfg.parallel = 'local';
cfg.timreq = 4; % in minutes
cfg.memreq = 1; % in GB
cfg.stack = 5;

EM_submit_to_tardis(cfg, inputfiles, outputfiles)
%TESTING:
% EM_submit_to_tardis(cfg, inputfiles(1:10), outputfiles(1:10))



