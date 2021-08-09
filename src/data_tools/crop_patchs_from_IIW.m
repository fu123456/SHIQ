% This script is to crop patches from input images in the original
% dataset and images ( diffuse images and specular highlight intensity
% images) produced by the adopted multi-image specular highlight
% removal method RPCA (see our paper). All cropped patches will be
% further screened by human subjects to produce desired patches and
% then we produce associated ground truth data. Note that, the pixel
% sizes of input, diffuse and highlight intensity images are
% necessarily equal for running this script.

clear all;
close all;
clc;

% You need to accordingly modify these variables "DATA_DIR", "result_dir"

DATA_DIR='<your data dir>'; % Input your dir for input data

% All subdirs of DATA_DIR
names{1}='<your data subdir>'; % Input your subdir

%% Parameters
% Width and height for cropped regions
m=200; % width
n=200; % height
s=100; % step size

result_dir='<your result dir>'; % input your dir for saving produced results

for i=1:numel(names)
    dataDIR=fullfile(DATA_DIR,names{i});
    dataFiles=dir(fullfile(dataDIR,'*.jpg'));
    for j=1:numel(dataFiles)
        [~,name,~]=fileparts(fullfile(dataDIR,dataFiles(j).name));
        if ~exist(fullfile(result_dir,names{i}))
            mkdir(fullfile(result_dir,names{i}));
        end
        img=imread(fullfile(dataDIR,dataFiles(j).name));
        disp(name);

        % Mkdir subdir by input image name
        if ~exist(fullfile(result_dir,names{i},name))
            mkdir(fullfile(result_dir,names{i},name));
        end

        [h,w,~]=size(img);
        kk=floor((w-m)/s)+1; % number of regions along y
        ll=floor((h-n)/s)+1; % number of regions along x

        % Begin cropping and save results to subdir
        for iii=1:kk % height
            for jjj=1:ll % width
                %% The top left x and y coordinates
                y=1+s*(iii-1);
                x=1+s*(jjj-1);
                crop_region=img(x:x+m-1,y:y+n-1,:);
                % Save patches with the png format or the format you need
                imwrite(crop_region,[result_dir '/' names{i} '/' name '/' sprintf('%02d_%02d_I.png',iii,jjj)]);
                fprintf('save crop_region to %02d_%02d.png\n',iii,jjj);
            end
        end
    end
end
