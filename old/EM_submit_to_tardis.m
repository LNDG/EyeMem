function EM_submit_to_tardis(cfg, varargin)

% Submit jobs to tardis with ft qsub
% argin: funhandle, cfg.parallel, cfg1 ...
% EXAMPLE:
% cfg = [];
% cfg.function_to_run = 'EM_eyelink_to_fieldtrip'; %         [asc, data, event] = EM_eyelink_to_fieldtrip(edf_name);
% cfg.compile = 'yes';
% cfg.parallel = 'torque';
% cfg.timreq = 5; % in minutes
% cfg.memreq = 0.5; % in GB
% cfg.stack = 1;
% EM_submit_to_tardis(cfg, inputfiles, outputfiles)


for i=1:length(varargin)
    disp(varargin(i))    
end

switch cfg.parallel
    case 'local'
%         cellfun( str2func(cfg.function_to_run), cfg1, cfg2, inputfile, outputfile);
        cellfun( str2func(cfg.function_to_run), varargin{:});
    case 'peer'
        peercellfun(str2func(cfg.function_to_run), varargin{:});
    case {'torque'}
        setenv('TORQUEHOME', 'yes')  %    yes or ''
        mkdir('~/qsub'); cd('~/qsub');
                
%         memreq = 8; % in GB
        switch cfg.compile
            case 'no'
%                 nnodes = 24; % how many licenses?
%                 stack = round(length(varargin{1})/nnodes);
%                 stack=1;
                qsubcellfun(str2func(cfg.function_to_run), varargin{:}, 'cfg.compile', 'no', ...
                    'memreq', cfg.memreq, 'timreq', cfg.timreq*60, 'stack', cfg.stack, 'StopOnError', false, 'backend', cfg.parallel, 'options', '-l nodes=1:ppn=1');
            case 'yes'
                cfg.compiledfun = qsubcompile(str2func(cfg.function_to_run), 'toolbox', {'signal', 'stats'});
                qsubcellfun(cfg.compiledfun, varargin{:}, ...
                    'memreq', cfg.memreq, 'timreq', cfg.timreq*60, 'stack', cfg.stack, 'StopOnError', false, 'backend', cfg.parallel, 'options', '-l nodes=1:ppn=1');
        end
    otherwise
        error('Unknown backend, aborting . . .\n')
end
