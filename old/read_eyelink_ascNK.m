function [data, event] = EM_read_eyelink_ascNK(filename)

% READ_EYELINK_ASC reads the header information, input triggers, messages
% and all data points from an Eyelink *.asc file and converts it into a
% fieldtrip data struct and produces the event structure.
%
% NK edit: add blink parsing
% Use as
%   asc = read_eyelink_asc(filename)

% Copyright (C) 2010, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: read_eyelink_asc.m 945 2010-04-21 17:41:20Z roboos $

% Step 1: %convert edf to ascii text file
edf2asc_path = ''; % TODO find edf2asc

[~,infile] = fileparts(batch(irun).dataset); % TODO
if ~exist([infile  '.asc'], 'file')
    unix(sprintf('%s -sg %s', edf2asc_path, infile)); %-y
end

% Step 2: %convert ascii to matlab struct

fid = fopen(filename, 'rt');

asc.header  = {};
asc.msg     = {};
asc.input   = [];
asc.sfix    = {};
asc.efix    = {};
asc.ssacc   = {};
asc.esacc   = {};
asc.sblink   = {}; % NK edit: add blink parsing
asc.eblink   = {};
asc.fsample = []; % NK edit: add fsample
asc.dat     = [];
current   = 0;

while ~feof(fid)
    tline = fgetl(fid);
    
    if regexp(tline, '^[0-9]');
        tmp   = sscanf(tline, '%f');
        nchan = numel(tmp);
        current = current + 1;
        
        if size(asc.dat,1)<nchan
            % increase the allocated number of channels
            asc.dat(nchan,:) = 0;
        end
        
        if size(asc.dat, 2)<current
            % increase the allocated number of samples
            asc.dat(:,end+10000) = 0;
        end
        
        % add the current sample to the data matrix
        asc.dat(1:nchan, current) = tmp;
        
        
    elseif regexp(tline, '^INPUT')
        [val, num] = sscanf(tline, 'INPUT %d %d');
        this.timestamp = val(1);
        this.value     = val(2);
        if isempty(asc.input)
            asc.input = this;
        else
            asc.input = cat(1, asc.input, this);
        end
        
        
    elseif regexp(tline, '\*\*.*')
        asc.header = cat(1, asc.header, {tline});
        
        
    elseif regexp(tline, '^MSG')
        asc.msg = cat(1, asc.msg, {tline});
        if regexp(tline, '!MODE RECORD CR')
            tok = tokenize(tline);
            asc.fsample = str2num(tok{6});
        end
        
        
    elseif regexp(tline, '^SFIX')
        asc.sfix = cat(1, asc.sfix, {tline});
        
        
    elseif regexp(tline, '^EFIX')
        asc.efix = cat(1, asc.efix, {tline});
        
        
    elseif regexp(tline, '^SSACC')
        asc.ssacc = cat(1, asc.ssacc, {tline});
        
        
    elseif regexp(tline, '^ESACC')
        asc.esacc = cat(1, asc.esacc, {tline});
        
        
    elseif regexp(tline, '^SBLINK')
        asc.sblink = cat(1, asc.sblink, {tline});
        
        
    elseif regexp(tline, '^EBLINK')
        asc.eblink = cat(1, asc.eblink, {tline});
        
        
        
        
    else
        % all other lines are not parsed
    end
    
end

% close the file?
fclose(fid);

% remove the samples that were not filled with real data
asc.dat = asc.dat(:,1:current);

% Step 3: turn asc struct into fieldtrip data structure
data = [];
data.label ={'h'; 'v'; 'p'};
data.trial = {ascdat.dat(2:end,:)};
data.fsample = ascdat.fsample;
data.time = {0:1/data.fsample:length(ascdat.dat(1,:))/data.fsample-1/data.fsample};
data.sampleinfo = [1 length(ascdat.dat(1,:))];

% Step 4: make event structure
if ~exist([infile  '_ascdat.mat'], 'file') % && overwrite
    disp('Reading asc . . .') % read asc file into matlab
    ascdat = read_eyelink_ascNK( [infile  '.asc'] );
    save([infile '_ascdat.mat'], 'ascdat')
    
    % create event structure
    evcell = cell(length(ascdat.msg),1);
    event=struct('type', evcell, 'sample', evcell, 'value', evcell, 'offset', evcell, 'duration', evcell );
    for i=1:length(ascdat.msg)
        strtok = tokenize(ascdat.msg{i});
        event(i).type = strtok{3};
        smpstamp = find(ascdat.dat(1,:) == str2double(strtok{2})); % find sample index of trigger in ascii dat
        if ~isempty(smpstamp)
            event(i).sample = smpstamp(1);
        else
            event(i).sample = nan; % give nan if no data was recorded at msg time
        end
        event(i).value = ascdat.msg{i}; % event(i).value = [strtok{3:end}]; %trigger value: e.g. ResponseMIBOff
    end
    save([infile '_event.mat'], 'event')
else
    load([infile '_ascdat.mat'])
    load([infile '_event.mat'])
end



