trlall = Eyemem_sorttrials();

%% plot ncorrect vs incorrect
size(trlall)

% trlold = trlall(trlall )

% subplot()

trl = permute(trlall, [4 3 1 2]); % condsubj trl col
trl = reshape(trl, 101, [], 14); % collapse over runs

test = sum(trl(:,10,:));
