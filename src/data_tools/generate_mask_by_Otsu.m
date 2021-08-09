% This script is to produce ground truth mask images using the
% modified Otsu algorithm (see "rolling_L_Otsu.m") on specular
% highlight layers.

clear all;
close all;
clc;

% Please modify two variables: (a) DATA_DIR; (b) output_dir
DATA_DIR='<your dir>'; % input your dir for input data
output_dir='<your dir>' % input your dir for saving results
if ~exist(output_dir)
    mkdir(output_dir)
end

% The number of Otsu iteration
num_otsu=2;

dataDIR=DATA_DIR;
dataFiles=dir(fullfile(dataDIR,'*_S.png')); % "_S": highlight images
for j=1:numel(dataFiles)
    [~,name,~]=fileparts(fullfile(dataDIR,dataFiles(j).name));
    disp(name);
    img=im2double(imread(fullfile(dataDIR,dataFiles(j).name)));
    % Rolling Otsu
    mask=Rolling_L_Otsu(img,num_otsu);
    result_path=output_dir;
    % Save binary mask image
    mask_name=strrep(name,'_S','_T'); % "_T" denotes mask images
    imwrite(mask,[result_path '/' mask_name '.png']);
end
