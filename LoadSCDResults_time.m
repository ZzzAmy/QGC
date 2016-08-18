
total_f1score = 0;
total_running_time = 0;
count = 0;

file=dir('/Users/iankuoli/Downloads/TT_SCD/*');

for i = 3:length(file)
    SCD_path = strcat('/Users/iankuoli/Downloads/code_nips2012_results2/', file(i).name);

    fid = fopen(SCD_path, 'r');
    cellData = textscan(fid, '%s%s', 'delimiter', '=');
    fclose(fid);

    str_cluster = cellData{1,2}{6};
    C = regexp(str_cluster,'],[','split');
    C{1,1} = C{1,1}(1,3:end);
    C{1,end} = C{1,end}(1,1:end-2);

    f1score = str2num(cellData{1,2}{5});
    running_time = str2num(cellData{1,2}{9});

    total_f1score = total_f1score + f1score;
    total_running_time = total_running_time + running_time;

    count = count + 1;
end

avg_f1score = total_f1score / count;
avg_running_time = total_running_time / count;
