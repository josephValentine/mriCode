%% Load Raw Data From Scanner
%{
reads raw data files from scanner into a cell array which is saved as a .mat
file. data{n,k} where n is the image # and k is the coil #. for example: to
extract all images from 5th coil you could just use data{:,5}.

coil data is not averaged here as this could mess up reconstructions which
rely on phase data.

the readMeasDataVB15 file must be located in the same folder as this file
%}
%--------------------------------------------------------------------------
%       file path for .dat files
%--------------------------------------------------------------------------
file_path = '/Users/uu/Desktop/rawdata/in_vivo/knee4';

%--------------------------------------------------------------------------
%       file path to save .mat file
%--------------------------------------------------------------------------
save_name = 'SSFP_knee4.mat';

%--------------------------------------------------------------------------
%       other options
%--------------------------------------------------------------------------
show_progress_in_console = true;

%display slices after saving file
show_slices = true;

%if showing slices, do we want to view 1st image from each coil?
%false will show each image with 1st coil
show_by_coil = false;

%combine averages with sum of squares? false is arithmetic avg
%FALSE FOR PHASE SENSITIVE TECHNIQUES!!
slice_avg_by_SoS = false;

%% Read Filenames

files = dir(fullfile(file_path,'*.dat'));

file_names = {};
for k=1:length(files)
    file_names{k} = fullfile(file_path,files(k).name);
end

%% Load Data
%initialize cell array with data for each image as row and each coil as col
img = readMeasDataVB15(file_names{1});
if(size(size(img),2)==4)%only 3 cols if no averages
    
    avg_num = size(img,3);
    coil_num = size(img,4);
else
    avg_num = 1;
    coil_num = size(img,3);
end

data = cell(length(file_names),coil_num);

for n=1:length(file_names)
    
    if show_progress_in_console
        text='Reading file #: %d';
        message=sprintf(text,n)
    end
    
    img = readMeasDataVB15(file_names{n});
    img = ifftshift(ifft2(ifftshift(img)));%transform from kspace
    
    for k=1:coil_num
        
        % get amount to add to avg_img for each coil(this is the avg val
        % from all averages)
        avg_img = zeros(size(img,1),size(img,2));
        if avg_num ~= 1
            if slice_avg_by_SoS
                
                for p=1:avg_num
                    avg_img = avg_img + img(:,:,p,k).^2;
                end
                avg_img = sqrt(avg_img);
                
            else
                for p=1:avg_num
                    avg_img = avg_img + img(:,:,p,k);
                end
                avg_img = avg_img/avg_num;
            end
        else
            avg_img = img(:,:,k);
        end
        
        
        %save in list in correct location
        data{(k-1)*length(file_names) + n} = avg_img;
    end
    
    
    
end

%save as mat file
save(save_name,'data');

%% Display Slices
if show_slices
    if show_by_coil
        for k=1:coil_num
            figure();
            subplot('121');
            imshow(abs(data{1,k}),[]);
            subplot('122');
            imshow(angle(data{1,k}),[]);
            
        end
    else
        for n=1:coil_num
            for k=1:length(file_names)
                figure();
                subplot('121');
                imshow(abs(data{k,n}),[]);
                string = sprintf('Magnitude of image %d coil %d',k,n);
                title(string);
                subplot('122');
                imshow(angle(data{k,n}),[]);
                string = sprintf('Phase of image %d coil %d',k,n);
                title(string);

            end
        end
    end
end
