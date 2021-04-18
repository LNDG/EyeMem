% EM_make_SUBJ_list


PREIN = '/Users/kloosterman/gridmaster2012/LNDG/EyeMem/data'

cd(PREIN)

subj_list = dir('EYEMEM*')

subj_str = [];
for isub = 1:length(subj_list)
    subj_str = cat(2, subj_str, sprintf(' ''%s'' ', subj_list(isub).name));
end