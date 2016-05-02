mag_file_path = '/Volumes/KINGSTON/5blindmice/Bangerter_Fatmouse/Research_Grayson/gre_te_661_deceased_27';
phase_file_path = '/Volumes/KINGSTON/5blindmice/Bangerter_Fatmouse/Research_Grayson/gre_te_661_deceased_28';

save_fileName = 'mouseDead661';
meta_fileName = 'mouseDeadMeta';



%% Load mag files (stored in a separate folder)

mag_files = dir(fullfile(mag_file_path,'*.dcm'));

Mfile_names = {};
mags={};
for k=1:length(mag_files)
    Mfile_names{k} = fullfile(mag_file_path,mag_files(k).name);
    mags{k} = im2double(dicomread(dicominfo(Mfile_names{k})));
    
end

%% Load phase files (stored in separate folder in the same order as mags)

phase_files = dir(fullfile(phase_file_path,'*.dcm'));

Pfile_names = {};
phase={};
for k=1:length(phase_files)
    Pfile_names{k} = fullfile(phase_file_path,phase_files(k).name);
    phase{k} = im2double(dicomread(dicominfo(Pfile_names{k})));%must use im2double so that we can use complex math
 
end



%% Put mag with phase

data = cell(length(Mfile_names),1);
for m=1:length(Mfile_names)
    %= zeros(size(s{2},1),size(s{2},2));
    data{m} = mags{m}.*exp(1i*phase{m});
end


save(save_fileName,'data');
meta = dicominfo(Mfile_names{1});
save(meta_fileName,'meta');
%{
for r=1:length(Mfile_names)
    figure();
    subplot(1,2,1);
    imshow(abs(data{r}),[]);
     subplot(1,2,2);
    imshow(angle(data{r}),[]);
end
%}