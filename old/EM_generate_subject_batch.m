function batch = EM_generate_subject_batch()

% Batch variable containing all subjects

batch = [];

% batch(1).SUBJ = 'EYEMEM001';
% batch(1).dateofbirth = '1-Jan-1956';
% batch(1).agegroup = 'old';
% batch(1).responsemapping = 1;
% batch(1).scandate = '1-Jan-1956';
%
% for isub=1:101
%     if isub < 10
%         batch(isub).SUBJ = sprintf('EYEMEM00%d', isub);
%     elseif isub < 100
%         batch(isub).SUBJ = sprintf('EYEMEM0%d', isub);
%     else
%         batch(isub).SUBJ = sprintf('EYEMEM%d', isub);
%     end
%
%     batch(isub).dateofbirth = 'PLACEHOLDER';
%     batch(isub).agegroup = 'old young';
%     batch(isub).responsemapping = 12;
%     batch(isub).scandate = '-2016';
%     batch(isub).gender = 'PLACEHOLDER';
%     batch(isub).handedness = 'PLACEHOLDER';
%
% end
%
% test = gencode(batch)

%%

batch(1).SUBJ = 'EYEMEM001'; % PILOT1
batch(1).dateofbirth = '8-Apr-1944';
batch(1).agegroup = 'old';
batch(1).responsemapping = 12;
batch(1).scandate = '13-Jun-2016';
batch(1).gender = 'm';
batch(1).handedness = '?';

batch(2).SUBJ = 'EYEMEM002'; % PILOT2: SUBJ sleeping arms
batch(2).dateofbirth = '10-Mar-1944';
batch(2).agegroup = 'old';
batch(2).responsemapping = 12;
batch(2).scandate = '13-Jun-2016';
batch(2).gender = 'f';
batch(2).handedness = 'r';

batch(3).SUBJ = 'EYEMEM003'; % PILOT3
batch(3).dateofbirth = '3-Apr-1947';
batch(3).agegroup = 'old';
batch(3).responsemapping = 12;
batch(3).scandate = '13-Jun-2016';
batch(3).gender = 'f';
batch(3).handedness = 'r';

batch(4).SUBJ = 'EYEMEM004';  % PILOT4
batch(4).dateofbirth = '7-Oct-1942';
batch(4).agegroup = 'old';
batch(4).responsemapping = 12;
batch(4).scandate = '13-Jun-2016';
batch(4).gender = 'm';
batch(4).handedness = 'r';

batch(5).SUBJ = 'EYEMEM005'; %4 runs, short resp duration: drop?
batch(5).dateofbirth = '7-Aug-1947';
batch(5).agegroup = 'old';
batch(5).responsemapping = 1;
batch(5).scandate = '14-June-2016';
batch(5).gender = 'f';
batch(5).handedness = 'r';

batch(6).SUBJ = 'EYEMEM006'; %4 runs, short resp duration: drop?
batch(6).dateofbirth = '20-Jun-1943';
batch(6).agegroup = 'old';
batch(6).responsemapping = 2;
batch(6).scandate = '14-Jun-2016';
batch(6).gender = 'm';
batch(6).handedness = 'r';

batch(7).SUBJ = 'EYEMEM007'; %4 runs, short resp duration: drop?
batch(7).dateofbirth = '22-Sep-1990';
batch(7).agegroup = 'young';
batch(7).responsemapping = 1;
batch(7).scandate = '14-Jun-2016';
batch(7).gender = 'm';
batch(7).handedness = '?';

batch(8).SUBJ = 'EYEMEM008'; %4 runs, short resp duration: drop?
batch(8).dateofbirth = '15-Jan-1945';
batch(8).agegroup = 'old';
batch(8).responsemapping = 2;
batch(8).scandate = '14-Jun-2016';
batch(8).gender = 'm';
batch(8).handedness = 'l';

batch(9).SUBJ = 'EYEMEM009';  %5 runs!, short resp duration: drop?
batch(9).dateofbirth = '3-Mar-1994';
batch(9).agegroup = 'young';
batch(9).responsemapping = 2;
batch(9).scandate = '14-Jun-2016';
batch(9).gender = 'f';
batch(9).handedness = '?';

batch(10).SUBJ = 'EYEMEM010'; % 5 runs, 2.5 s resp duration!
batch(10).dateofbirth = '5-May-1987'; % good et quality
batch(10).agegroup = 'young';
batch(10).responsemapping = 1;
batch(10).scandate = '14-Jun-2016';
batch(10).gender = 'f';
batch(10).handedness = 'r';

batch(11).SUBJ = 'EYEMEM011';
batch(11).dateofbirth = '27-Jan-1992';
batch(11).agegroup = 'young';
batch(11).responsemapping = 2;
batch(11).scandate = '15-Jun-2016';
batch(11).gender = 'f';
batch(11).handedness = 'r';

batch(12).SUBJ = 'EYEMEM012';
batch(12).dateofbirth = '17-Aug-1943';
batch(12).agegroup = 'old';
batch(12).responsemapping = 2;
batch(12).scandate = '20-Jun-2016';
batch(12).gender = 'm';
batch(12).handedness = 'r';

batch(13).SUBJ = 'EYEMEM013'; % wrong button pressed in run1 and 2
batch(13).dateofbirth = '15-Nov-1946';
batch(13).agegroup = 'old';
batch(13).responsemapping = 2;
batch(13).scandate = '20-Jun-2016';
batch(13).gender = 'f';
batch(13).handedness = 'r';

batch(14).SUBJ = 'EYEMEM014';
batch(14).dateofbirth = '11-Feb-1948';
batch(14).agegroup = 'old';
batch(14).responsemapping = 2;
batch(14).scandate = '20-Jun-2016';
batch(14).gender = 'm';
batch(14).handedness = 'r';

batch(15).SUBJ = 'EYEMEM015'; % mri motion
batch(15).dateofbirth = '19-Jan-1988';
batch(15).agegroup = 'young';
batch(15).responsemapping = 2;
batch(15).scandate = '20-Jun-2016';
batch(15).gender = 'f';
batch(15).handedness = 'r';

batch(16).SUBJ = 'EYEMEM016';
batch(16).dateofbirth = '13-Jan-1943';
batch(16).agegroup = 'old';
batch(16).responsemapping = 1;
batch(16).scandate = '21-Jun-2016';
batch(16).gender = 'f';
batch(16).handedness = 'r';

batch(17).SUBJ = 'EYEMEM017';
batch(17).dateofbirth = '4-Jun-1988';
batch(17).agegroup = 'young';
batch(17).responsemapping = 1;
batch(17).scandate = '21-Jun-2016';
batch(17).gender = 'm';
batch(17).handedness = 'r';

batch(18).SUBJ = 'EYEMEM018';
batch(18).dateofbirth = '7-Oct-1990';
batch(18).agegroup = 'young';
batch(18).responsemapping = 1;
batch(18).scandate = '21-Jun-2016';
batch(18).gender = 'f';
batch(18).handedness = 'r';

batch(19).SUBJ = 'EYEMEM019';
batch(19).dateofbirth = '26-Jun-1990';
batch(19).agegroup = 'young';
batch(19).responsemapping = 1;
batch(19).scandate = '21-Jun-2016';
batch(19).gender = 'm';
batch(19).handedness = 'r';

batch(20).SUBJ = 'EYEMEM020'; %motion
batch(20).dateofbirth = '21-Jun-1993';
batch(20).agegroup = 'young';
batch(20).responsemapping = 1;
batch(20).scandate = '24-Jun-2016';
batch(20).gender = 'f';
batch(20).handedness = 'l';

batch(21).SUBJ = 'EYEMEM021';
batch(21).dateofbirth = '14-Feb-1947';
batch(21).agegroup = 'old';
batch(21).responsemapping = 1;
batch(21).scandate = '24-Jun-2016';
batch(21).gender = 'f';
batch(21).handedness = 'r';

batch(22).SUBJ = 'EYEMEM022'; % feeling sick resting state, blinks
batch(22).dateofbirth = '28-Sep-1993';
batch(22).agegroup = 'young';
batch(22).responsemapping = 1;
batch(22).scandate = '24-Jun-2016';
batch(22).gender = 'm';
batch(22).handedness = 'r';

batch(23).SUBJ = 'EYEMEM023';
batch(23).dateofbirth = '25-Dec-1942';
batch(23).agegroup = 'old';
batch(23).responsemapping = 1;
batch(23).scandate = '27-Jun-2016';
batch(23).gender = 'f';
batch(23).handedness = 'r';

batch(24).SUBJ = 'EYEMEM024';
batch(24).dateofbirth = '27-Dec-1950';
batch(24).agegroup = 'old';
batch(24).responsemapping = 2;
batch(24).scandate = '27-Jun-2016';
batch(24).gender = 'f';
batch(24).handedness = 'r';

batch(25).SUBJ = 'EYEMEM025';
batch(25).dateofbirth = '27-Dec-1950';
batch(25).agegroup = 'old';
batch(25).responsemapping = 2;
batch(25).scandate = '27-Jun-2016';
batch(25).gender = 'm';
batch(25).handedness = 'r';

batch(26).SUBJ = 'EYEMEM026';
batch(26).dateofbirth = '27-Apr-1948';
batch(26).agegroup = 'old';
batch(26).responsemapping = 2;
batch(26).scandate = '27-Jun-2016';
batch(26).gender = 'm';
batch(26).handedness = 'r';

batch(27).SUBJ = 'EYEMEM027';
batch(27).dateofbirth = '10-Nov-1944';
batch(27).agegroup = 'old';
batch(27).responsemapping = 2;
batch(27).scandate = '28-Jun-2016';
batch(27).gender = 'm';
batch(27).handedness = 'r';

batch(28).SUBJ = 'EYEMEM028';
batch(28).dateofbirth = '2-Sep-1943';
batch(28).agegroup = 'old';
batch(28).responsemapping = 2;
batch(28).scandate = '28-Jun-2016';
batch(28).gender = 'm';
batch(28).handedness = 'r';

batch(29).SUBJ = 'EYEMEM029';
batch(29).dateofbirth = '18-Mar-1946';
batch(29).agegroup = 'old';
batch(29).responsemapping = 1;
batch(29).scandate = '28-Jun-2016';
batch(29).gender = 'f';
batch(29).handedness = 'r';

batch(30).SUBJ = 'EYEMEM030';
batch(30).dateofbirth = '11-Apr-1944';
batch(30).agegroup = 'old';
batch(30).responsemapping = 2;
batch(30).scandate = '28-Jun-2016';
batch(30).gender = 'm';
batch(30).handedness = 'r';

batch(31).SUBJ = 'EYEMEM031';
batch(31).dateofbirth = '10-Apr-1987';
batch(31).agegroup = 'young';
batch(31).responsemapping = 12;
batch(31).scandate = '28-Jun-2016';
batch(31).gender = 'm';
batch(31).handedness = 'both';

batch(32).SUBJ = 'EYEMEM032';
batch(32).dateofbirth = '17-Jul-1985';
batch(32).agegroup = 'young';
batch(32).responsemapping = 1;
batch(32).scandate = '29-Jun-2016';
batch(32).gender = 'f';
batch(32).handedness = 'r';

batch(33).SUBJ = 'EYEMEM033';
batch(33).dateofbirth = '31-Jul-1988';
batch(33).agegroup = 'young';
batch(33).responsemapping = 1;
batch(33).scandate = '1-Jul-2016';
batch(33).gender = 'f';
batch(33).handedness = 'r';

batch(34).SUBJ = 'EYEMEM034';
batch(34).dateofbirth = '8-May-1944';
batch(34).agegroup = 'old';
batch(34).responsemapping = 1;
batch(34).scandate = '1-Jul-2016';
batch(34).gender = 'm';
batch(34).handedness = 'r';

batch(35).SUBJ = 'EYEMEM035';
batch(35).dateofbirth = '12-Apr-1944';
batch(35).agegroup = 'old';
batch(35).responsemapping = 2;
batch(35).scandate = '1-Jul-2016';
batch(35).gender = 'f';
batch(35).handedness = 'r';

batch(36).SUBJ = 'EYEMEM036';
batch(36).dateofbirth = '8-Dec-1946';
batch(36).agegroup = 'old';
batch(36).responsemapping = 2;
batch(36).scandate = '1-Jul-2016';
batch(36).gender = 'm';
batch(36).handedness = 'r';

batch(37).SUBJ = 'EYEMEM037';
batch(37).dateofbirth = '4-Aug-1944';
batch(37).agegroup = 'old';
batch(37).responsemapping = 1;
batch(37).scandate = '1-Jul-2016';
batch(37).gender = 'm';
batch(37).handedness = 'r';

batch(38).SUBJ = 'EYEMEM038';
batch(38).dateofbirth = '3-Oct-1944';
batch(38).agegroup = 'old';
batch(38).responsemapping = 2;
batch(38).scandate = '1-Jul-2016';
batch(38).gender = 'm';
batch(38).handedness = 'both';

batch(39).SUBJ = 'EYEMEM039';
batch(39).dateofbirth = '7-Nov-1986';
batch(39).agegroup = 'young';
batch(39).responsemapping = 2;
batch(39).scandate = '4-Jul-2016';
batch(39).gender = 'f';
batch(39).handedness = 'r';

batch(40).SUBJ = 'EYEMEM040';
batch(40).dateofbirth = '28-Apr-1990';
batch(40).agegroup = 'young';
batch(40).responsemapping = 2;
batch(40).scandate = '4-Jul-2016';
batch(40).gender = 'f';
batch(40).handedness = 'r';

batch(41).SUBJ = 'EYEMEM041';
batch(41).dateofbirth = '23-Aug-1943';
batch(41).agegroup = 'old';
batch(41).responsemapping = 1;
batch(41).scandate = '4-Jul-2016';
batch(41).gender = 'f';
batch(41).handedness = 'r';

batch(42).SUBJ = 'EYEMEM042';
batch(42).dateofbirth = '2-Dec-1940';
batch(42).agegroup = 'old';
batch(42).responsemapping = 2;
batch(42).scandate = '4-Jul-2016';
batch(42).gender = 'f';
batch(42).handedness = 'r';

batch(43).SUBJ = 'EYEMEM043';
batch(43).dateofbirth = '27-Feb-1943';
batch(43).agegroup = 'old';
batch(43).responsemapping = 2;
batch(43).scandate = '4-Jul-2016';
batch(43).gender = 'f';
batch(43).handedness = 'r';

batch(44).SUBJ = 'EYEMEM044';
batch(44).dateofbirth = '5-Dec-1987';
batch(44).agegroup = 'young';
batch(44).responsemapping = 1;
batch(44).scandate = '4-Jul-2016';
batch(44).gender = 'm';
batch(44).handedness = 'r';

batch(45).SUBJ = 'EYEMEM045';
batch(45).dateofbirth = '23-May-1942';
batch(45).agegroup = 'old';
batch(45).responsemapping = 1;
batch(45).scandate = '5-Jul-2016';
batch(45).gender = 'f';
batch(45).handedness = 'r';

batch(46).SUBJ = 'EYEMEM046';
batch(46).dateofbirth = '12-Aug-1941';
batch(46).agegroup = 'old';
batch(46).responsemapping = 2;
batch(46).scandate = '5-Jul-2016';
batch(46).gender = 'f';
batch(46).handedness = 'r';

batch(47).SUBJ = 'EYEMEM047'; % Distortion right frontal
batch(47).dateofbirth = '26-Jul-1986';
batch(47).agegroup = 'young';
batch(47).responsemapping = 1;
batch(47).scandate = '5-Jul-2016';
batch(47).gender = 'f';
batch(47).handedness = 'r';

batch(48).SUBJ = 'EYEMEM048';
batch(48).dateofbirth = '13-Aug-1941';
batch(48).agegroup = 'old';
batch(48).responsemapping = 1;
batch(48).scandate = '5-Jul-2016';
batch(48).gender = 'f';
batch(48).handedness = 'r';

batch(49).SUBJ = 'EYEMEM049';
batch(49).dateofbirth = '14-Nov-1987';
batch(49).agegroup = 'young';
batch(49).responsemapping = 1;
batch(49).scandate = '5-Jul-2016';
batch(49).gender = 'f';
batch(49).handedness = 'r';

batch(50).SUBJ = 'EYEMEM050';
batch(50).dateofbirth = '30-Jan-1996';
batch(50).agegroup = 'young';
batch(50).responsemapping = 2;
batch(50).scandate = '5-Jul-2016';
batch(50).gender = 'm';
batch(50).handedness = 'r';

batch(51).SUBJ = 'EYEMEM051';
batch(51).dateofbirth = '14-Jul-1987';
batch(51).agegroup = 'young';
batch(51).responsemapping = 1;
batch(51).scandate = '6-Jul-2016';
batch(51).gender = 'm';
batch(51).handedness = 'r';

batch(52).SUBJ = 'EYEMEM052';
batch(52).dateofbirth = '3-Aug-1991';
batch(52).agegroup = 'young';
batch(52).responsemapping = 1;
batch(52).scandate = '8-Jul-2016';
batch(52).gender = 'f';
batch(52).handedness = 'l';

batch(53).SUBJ = 'EYEMEM053';
batch(53).dateofbirth = '22-Jan-1988';
batch(53).agegroup = 'young';
batch(53).responsemapping = 1;
batch(53).scandate = '8-Jul-2016';
batch(53).gender = 'm';
batch(53).handedness = 'r';

batch(54).SUBJ = 'EYEMEM054'; % response button switched after run 1
batch(54).dateofbirth = '25-Nov-1944';
batch(54).agegroup = 'old';
batch(54).responsemapping = 1;
batch(54).scandate = '8-Jul-2016';
batch(54).gender = 'm';
batch(54).handedness = 'l';

batch(55).SUBJ = 'EYEMEM055'; % Memory task: subject did not know one could press until 3 s!
batch(55).dateofbirth = '25-Oct-1993'; % run 1/2 many blinks
batch(55).agegroup = 'young';
batch(55).responsemapping = 2;
batch(55).scandate = '8-Jul-2016';
batch(55).gender = 'f';
batch(55).handedness = 'r';

batch(56).SUBJ = 'EYEMEM056'; % runs 2-5 done first, and then run 1! Edf's renamed on LNDG
batch(56).dateofbirth = '24-Nov-1993'; % run1 done first in memtest
batch(56).agegroup = 'young';
batch(56).responsemapping = 1;
batch(56).scandate = '8-Jul-2016';
batch(56).gender = 'f';
batch(56).handedness = 'r';

batch(57).SUBJ = 'EYEMEM057'; % Many short blinks
batch(57).dateofbirth = '15-May-1988';
batch(57).agegroup = 'young';
batch(57).responsemapping = 2;
batch(57).scandate = '8-Jul-2016';
batch(57).gender = 'f';
batch(57).handedness = 'r';

batch(58).SUBJ = 'EYEMEM058'; % respiration sometimes clips
batch(58).dateofbirth = '22-Aug-1992';
batch(58).agegroup = 'young';
batch(58).responsemapping = 2;
batch(58).scandate = '9-Jul-2016';
batch(58).gender = 'f';
batch(58).handedness = 'r';

batch(59).SUBJ = 'EYEMEM059'; %KAum blinking!
batch(59).dateofbirth = '14-Dec-1989';
batch(59).agegroup = 'young';
batch(59).responsemapping = 2;
batch(59).scandate = '9-Jul-2016';
batch(59).gender = 'f';
batch(59).handedness = 'r';

batch(60).SUBJ = 'EYEMEM060'; %Drop subject! Too Sleepy, ET don't work. Aorted after run 3
batch(60).dateofbirth = '1-Jan-1941'; % resting state eye sometimes lost, shouted 3 times to keep her awake
batch(60).agegroup = 'old';
batch(60).responsemapping = 2;
batch(60).scandate = '11-Jul-2016';
batch(60).gender = 'f';
batch(60).handedness = 'r';

batch(61).SUBJ = 'EYEMEM061'; % complained of teary eyes
batch(61).dateofbirth = '14-May-1949';
batch(61).agegroup = 'old';
batch(61).responsemapping = 1;
batch(61).scandate = '11-Jul-2016';
batch(61).gender = 'm';
batch(61).handedness = 'r';

batch(62).SUBJ = 'EYEMEM062'; % Eye fixation, blink problems
batch(62).dateofbirth = '29-Aug-1988';
batch(62).agegroup = 'young';
batch(62).responsemapping = 2;
batch(62).scandate = '11-Jul-2016';
batch(62).gender = 'f';
batch(62).handedness = 'r';

batch(63).SUBJ = 'EYEMEM063'; % fan to 1 always now
batch(63).dateofbirth = '4-Jan-1943'; % big head, bad eye quality
batch(63).agegroup = 'old'; %stimulus pc hangup, exp aborted after run 2, drop subject?
batch(63).responsemapping = 1;
batch(63).scandate = '11-Jul-2016';
batch(63).gender = 'm';
batch(63).handedness = 'r';

batch(64).SUBJ = 'EYEMEM064'; % subject coughed, motion artf?
batch(64).dateofbirth = '19-May-1988';
batch(64).agegroup = 'young';
batch(64).responsemapping = 2;
batch(64).scandate = '11-Jul-2016';
batch(64).gender = 'f';
batch(64).handedness = 'r';

batch(65).SUBJ = 'EYEMEM065'; % very good subject
batch(65).dateofbirth = '1-Apr-1988';
batch(65).agegroup = 'young';
batch(65).responsemapping = 2;
batch(65).scandate = '12-Jul-2016';
batch(65).gender = 'm';
batch(65).handedness = 'r';

batch(66).SUBJ = 'EYEMEM066'; % Study phase started instead of resting state. So two resting state present, the 1st is run 1. Is trimmed to 474 volumes
batch(66).dateofbirth = '18-Apr-1948';
batch(66).agegroup = 'old';
batch(66).responsemapping = 2;
batch(66).scandate = '12-Jul-2016';
batch(66).gender = 'f';
batch(66).handedness = 'r';

batch(67).SUBJ = 'EYEMEM067'; %many blinks but not during pictures
batch(67).dateofbirth = '24-Sep-1991';
batch(67).agegroup = 'young';
batch(67).responsemapping = 2;
batch(67).scandate = '12-Jul-2016';
batch(67).gender = 'm';
batch(67).handedness = 'r';

batch(68).SUBJ = 'EYEMEM068'; % blinks a lot, mostly not during pics
batch(68).dateofbirth = '9-Jan-1990';
batch(68).agegroup = 'young';
batch(68).responsemapping = 2;
batch(68).scandate = '12-Jul-2016';
batch(68).gender = 'f';
batch(68).handedness = 'r';

batch(69).SUBJ = 'EYEMEM069';
batch(69).dateofbirth = '31-May-1991';
batch(69).agegroup = 'young';
batch(69).responsemapping = 1;
batch(69).scandate = '13-Jul-2016';
batch(69).gender = 'm';
batch(69).handedness = 'r';

batch(70).SUBJ = 'EYEMEM070'; % Asian, check for registration
batch(70).dateofbirth = '6-Sep-1989';
batch(70).agegroup = 'young';
batch(70).responsemapping = 2;
batch(70).scandate = '13-Jul-2016';
batch(70).gender = 'm';
batch(70).handedness = 'r';

batch(71).SUBJ = 'EYEMEM071';
batch(71).dateofbirth = '16-Jan-1995';
batch(71).agegroup = 'young';
batch(71).responsemapping = 2;
batch(71).scandate = '14-Jul-2016';
batch(71).gender = 'm';
batch(71).handedness = 'r';

batch(72).SUBJ = 'EYEMEM072'; % contact lenses, many small blinks, check for motion
batch(72).dateofbirth = '3-Oct-1994';
batch(72).agegroup = 'young';
batch(72).responsemapping = 2;
batch(72).scandate = '14-Jul-2016';
batch(72).gender = 'm';
batch(72).handedness = 'r';

batch(73).SUBJ = 'EYEMEM073';
batch(73).dateofbirth = '4-Apr-1995';
batch(73).agegroup = 'young';
batch(73).responsemapping = 1;
batch(73).scandate = '14-Jul-2016';
batch(73).gender = 'm';
batch(73).handedness = 'r';

batch(74).SUBJ = 'EYEMEM074'; % subj tired, some motion
batch(74).dateofbirth = '2-May-1994';
batch(74).agegroup = 'young';
batch(74).responsemapping = 1;
batch(74).scandate = '15-Jul-2016';
batch(74).gender = 'm';
batch(74).handedness = 'r';

batch(75).SUBJ = 'EYEMEM075';
batch(75).dateofbirth = '22-Jul-1946';
batch(75).agegroup = 'old';
batch(75).responsemapping = 1;
batch(75).scandate = '18-Jul-2016';
batch(75).gender = 'f';
batch(75).handedness = 'r';

batch(76).SUBJ = 'EYEMEM076';
batch(76).dateofbirth = '6-Jul-1990';
batch(76).agegroup = 'young';
batch(76).responsemapping = 2;
batch(76).scandate = '18-Jul-2016';
batch(76).gender = 'm';
batch(76).handedness = 'r';

batch(77).SUBJ = 'EYEMEM077';
batch(77).dateofbirth = '13-Sep-1991';
batch(77).agegroup = 'young';
batch(77).responsemapping = 1;
batch(77).scandate = '18-Jul-2016';
batch(77).gender = 'm';
batch(77).handedness = 'r';

batch(78).SUBJ = 'EYEMEM078'; % tired
batch(78).dateofbirth = '25-Aug-1988';
batch(78).agegroup = 'young';
batch(78).responsemapping = 1;
batch(78).scandate = '18-Jul-2016';
batch(78).gender = 'm';
batch(78).handedness = '?';

batch(79).SUBJ = 'EYEMEM079'; % few blinks, some motion
batch(79).dateofbirth = '20-May-1945';
batch(79).agegroup = 'old';
batch(79).responsemapping = 2;
batch(79).scandate = '19-Jul-2016';
batch(79).gender = 'm';
batch(79).handedness = 'r';

batch(80).SUBJ = 'EYEMEM080';
batch(80).dateofbirth = '20-Oct-1990';
batch(80).agegroup = 'young';
batch(80).responsemapping = 2;
batch(80).scandate = '19-Jul-2016';
batch(80).gender = 'm';
batch(80).handedness = 'r';

batch(81).SUBJ = 'EYEMEM081'; % ET loses pupil sometimes (rs, run 3)
batch(81).dateofbirth = '18-Jun-1948';
batch(81).agegroup = 'old';
batch(81).responsemapping = 2;
batch(81).scandate = '19-Jul-2016';
batch(81).gender = 'm';
batch(81).handedness = 'r';

batch(82).SUBJ = 'EYEMEM082'; % run2: many blinks
batch(82).dateofbirth = '29-Aug-1944';
batch(82).agegroup = 'old';
batch(82).responsemapping = 1;
batch(82).scandate = '19-Jul-2016';
batch(82).gender = 'f';
batch(82).handedness = 'r';

batch(83).SUBJ = 'EYEMEM083';
batch(83).dateofbirth = '16-Nov-1987';
batch(83).agegroup = 'young';
batch(83).responsemapping = 2;
batch(83).scandate = '19-Jul-2016';
batch(83).gender = 'f';
batch(83).handedness = 'r';

batch(84).SUBJ = 'EYEMEM084'; % some motion, blinks, fixation breaks
batch(84).dateofbirth = '14-Apr-1941'; % complains about neck problems
batch(84).agegroup = 'old';
batch(84).responsemapping = 2;
batch(84).scandate = '22-Jul-2016';
batch(84).gender = 'm';
batch(84).handedness = 'r';

batch(85).SUBJ = 'EYEMEM085'; % almost no blinks
batch(85).dateofbirth = '17-Nov-1942';
batch(85).agegroup = 'old';
batch(85).responsemapping = 1;
batch(85).scandate = '22-Jul-2016';
batch(85).gender = 'f';
batch(85).handedness = 'r';

batch(86).SUBJ = 'EYEMEM086'; % some motion
batch(86).dateofbirth = '22-Apr-1952';
batch(86).agegroup = 'old';
batch(86).responsemapping = 2;
batch(86).scandate = '25-Jul-2016';
batch(86).gender = 'm';
batch(86).handedness = 'r';

batch(87).SUBJ = 'EYEMEM087';
batch(87).dateofbirth = '3-Feb-1949';
batch(87).agegroup = 'old';
batch(87).responsemapping = 1;
batch(87).scandate = '25-Jul-2016';
batch(87).gender = 'm';
batch(87).handedness = 'r';

batch(88).SUBJ = 'EYEMEM088'; %some blinking, closes eyes sometimes
batch(88).dateofbirth = '23-Sep-1985';
batch(88).agegroup = 'young';
batch(88).responsemapping = 1;
batch(88).scandate = '25-Jul-2016';
batch(88).gender = 'f';
batch(88).handedness = 'l';

batch(89).SUBJ = 'EYEMEM089'; % big head, small pillow, blinzelt sehr oft, drop?
batch(89).dateofbirth = '3-Nov-1942';
batch(89).agegroup = 'old';
batch(89).responsemapping = 1;
batch(89).scandate = '25-Jul-2016';
batch(89).gender = 'm';
batch(89).handedness = 'r';

batch(90).SUBJ = 'EYEMEM090'; % run5 blinks more often
batch(90).dateofbirth = '4-Jan-1992';
batch(90).agegroup = 'young';
batch(90).responsemapping = 2;
batch(90).scandate = '25-Jul-2016';
batch(90).gender = 'm';
batch(90).handedness = 'r';

batch(91).SUBJ = 'EYEMEM091'; % reran localizer + t1 b/c subject's jacket zipper distortion
batch(91).dateofbirth = '4-Jan-1944';
batch(91).agegroup = 'old';
batch(91).responsemapping = 1;
batch(91).scandate = '26-Jul-2016';
batch(91).gender = 'f';
batch(91).handedness = 'r';

batch(92).SUBJ = 'EYEMEM092'; %many blinks (contacts), ET quality unstable
batch(92).dateofbirth = '16-Mar-1989';
batch(92).agegroup = 'young';
batch(92).responsemapping = 1;
batch(92).scandate = '26-Jul-2016';
batch(92).gender = 'f';
batch(92).handedness = 'r';

batch(93).SUBJ = 'EYEMEM093'; % Very good subject
batch(93).dateofbirth = '13-Sep-1991';
batch(93).agegroup = 'young';
batch(93).responsemapping = 2;
batch(93).scandate = '27-Jul-2016';
batch(93).gender = 'm';
batch(93).handedness = 'r';

batch(94).SUBJ = 'EYEMEM094';
batch(94).dateofbirth = '2-Jun-1942';
batch(94).agegroup = 'old';
batch(94).responsemapping = 1;
batch(94).scandate = '28-Jul-2016';
batch(94).gender = 'f';
batch(94).handedness = 'r';

batch(95).SUBJ = 'EYEMEM095'; % quite some blinking, some backpains
batch(95).dateofbirth = '01-Jul-1996'; % 20 years, birthdate unknown
batch(95).agegroup = 'young';
batch(95).responsemapping = 2;
batch(95).scandate = '18-Aug-2016';
batch(95).gender = 'f';
batch(95).handedness = 'r';

batch(96).SUBJ = 'EYEMEM096'; % bad et rs, runs1-3, 4-5 much better; some motion
batch(96).dateofbirth = '16-Dec-1943';
batch(96).agegroup = 'old';
batch(96).responsemapping = 1;
batch(96).scandate = '18-Aug-2016';
batch(96).gender = 'f';
batch(96).handedness = 'r';

batch(97).SUBJ = 'EYEMEM097'; % tired, many blinks, run5 missed a few trials
batch(97).dateofbirth = '20-Jul-1990';
batch(97).agegroup = 'young';
batch(97).responsemapping = 1;
batch(97).scandate = '18-Aug-2016';
batch(97).gender = 'f';
batch(97).handedness = '?'; % unknown

batch(98).SUBJ = 'EYEMEM098'; % all supi
batch(98).dateofbirth = '6-Apr-1987';
batch(98).agegroup = 'young';
batch(98).responsemapping = 2;
batch(98).scandate = '18-Aug-2016';
batch(98).gender = 'm';
batch(98).handedness = 'r';


    batch(99).SUBJ = 'EYEMEM099'; % Small eyes, ET impossible, drop?
    batch(99).dateofbirth = '22-Sep-1943';
    batch(99).agegroup = 'old';
    batch(99).responsemapping = 1;
    batch(99).scandate = '19-Aug-2016';
    batch(99).gender = 'm';
    batch(99).handedness = 'r';


batch(100).SUBJ = 'EYEMEM100'; % many blinks, tired at the end
batch(100).dateofbirth = '12-Aug-1948';
batch(100).agegroup = 'old';
batch(100).responsemapping = 1;
batch(100).scandate = '19-Aug-2016';
batch(100).gender = 'm';
batch(100).handedness = 'r';

batch(101).SUBJ = 'EYEMEM101'; % run3: taped off button pressed, score only 28 %. Drop run?
batch(101).dateofbirth = '19-Feb-1948';
batch(101).agegroup = 'old';
batch(101).responsemapping = 1;
batch(101).scandate = '19-Aug-2016';
batch(101).gender = 'm';
batch(101).handedness = 'r';

fprintf('\nN dropped subjects: %d\n\n', length(find(cellfun(@isempty, {batch.SUBJ}) == 1)) ) 
