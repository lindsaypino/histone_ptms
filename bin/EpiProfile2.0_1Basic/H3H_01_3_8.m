function H3H_01_3_8(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%

% check
out_filename = 'H3_01_3_8';
fprintf(1,'%s..',out_filename);
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    return;
end;

% init
His = init_histone(special);

% relocate
His = relocate(cur_outpath,out_filename,His);

% calculate
unitdiff = 1.0032;
[pep_rts,pep_intens,mono_isointens] = calculate_layout(MS1_index,MS1_peaks,ptol,unitdiff,His,special);

% output
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);

% draw
num_MS1 = size(MS1_index,1);
isorts = MS1_index(1:num_MS1,2);
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts,mono_isointens,MS2_index,MS2_peaks,special);

% Get PSM
if 1==special.nDAmode
    GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end;

function His = init_histone(special)
%%

His.pep_seq = 'TKQTAR';
His.mod_short = {'unmod';
    'K4me1';
    'K4me2';
    'K4me3';
    'K4ac'};
if 5==special.nsubtype
    His.mod_type = {'0,pr;2,pr;';
        '0,pr;2,hme1;';
        '0,pr;2,hme2;';
        '0,pr;2,hme3;';
        '0,pr;2,hac;'};
else
    His.mod_type = {'0,pr;2,pr;';
        '0,pr;2,me1;';
        '0,pr;2,me2;';
        '0,pr;2,me3;';
        '0,pr;2,ac;'};
end;

His.pep_ch = repmat([1 2],length(His.mod_type),1);
%{
His.pep_mz = [816.4574	408.7323
    830.4730	415.7402
    788.4625	394.7349
    802.4781	401.7427
    802.4417	401.7245];
%}
His.pep_mz = calculate_pepmzH(His);
His.rt_ref = [18.34
    21.58
    10.82
    10.8
    16.64];
His.display = ones(length(His.mod_type),1);

% main ch
main_ch = His.pep_ch(1,2);
if main_ch~=His.pep_ch(1,1)
    [npep,ncharge] = size(His.pep_mz);%#ok
    new_ch = [main_ch,setdiff(His.pep_ch(1,:),main_ch)];
    x = zeros([1,ncharge]);
    for ino=1:ncharge
        x(ino) = find(His.pep_ch(1,:)==new_ch(ino));
    end;
    tune = [3 4];%1:npep;
    His.pep_mz(tune,:) = His.pep_mz(tune,x);
    His.pep_ch(tune,:) = His.pep_ch(tune,x);
end;

function His0 = relocate(cur_outpath,out_filename,His0)
%%

[path1,name1] = fileparts(cur_outpath);
[path2,name2] = fileparts(path1);
out_file1 = fullfile(fullfile(fileparts(path2),name2,name1),[out_filename,'.mat']);
if 0==exist(out_file1,'file')
    fprintf(1,'%s: not exist.\n',out_file1);
    return;
end;
load(out_file1);
for ino=1:length(His0.rt_ref)
    His0.rt_ref(ino) = auc(ino,1);
end;

function [pep_rts,pep_intens,mono_isointens] = calculate_layout(MS1_index,MS1_peaks,ptol,unitdiff,His,special)%#ok
%%

[npep,ncharge] = size(His.pep_mz);
num_MS1 = size(MS1_index,1);
pep_rts = zeros([npep,ncharge]);
pep_intens = zeros([npep,ncharge]);
mono_isointens = zeros([num_MS1,npep]);

% unmod
% K4me1
% K4me2
% K4me3
% K4ac
for hno=1:5
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge) = cur_rts;
        pep_intens(hno,1:ncharge) = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end;
end;