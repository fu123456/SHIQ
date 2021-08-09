% This script is to perform white balance operation for our dataset
% due to the color distortion of ground truth specular-free image. We
% found that the network trained on the data (input and associated
% ground truth highlight-free images) without color balance generally
% produces perceptually unsatisfactory results. This operation can be
% effectively alleviate the this issue.

clear all;
close all;
clc;

% Please modify two variables: (a) DATA_DIR; (b) output_dir
DATA_DIR='<your dir>'; % Input your dir for input data
output_dir='<your dir>'; % input your dir for saving results
if ~exist(output_dir)
    mkdir(output_dir)
end

% "_A" denotes input highlight images; "_D" denotes highlight-free
% images; "_S" denotes highlight intensity images. Please modify the
% all names of the data with this standard
dataDIR=DATA_DIR;
dataFiles=dir(fullfile(dataDIR,'*_S.png')); % highlight images
for j=1:numel(dataFiles)
    [~,name,~]=fileparts(fullfile(dataDIR,dataFiles(j).name));
    diffuse_name=strrep(name,'_S','_D'); % name of specular-free image
    input_name=strrep(name,'_S','_A'); % name of input image
    new_generate_name=strrep(name,'_S','_A'); % name of new generated image
    disp(name);
    highlight=im2double(imread(fullfile(dataDIR,dataFiles(j).name)));
    diffuse=im2double(imread(fullfile(dataDIR,[diffuse_name '.png'])));
    new_generate_image=diffuse+repmat(rgb2gray(highlight),[1 1 3]);
    result_path=output_dir;
    % Save new generated images
    imwrite(diffuse,[result_path '/' diffuse_name '.png']);
    imwrite(new_generate_image,[result_path '/' new_generate_name '.png']);
    imwrite(repmat(rgb2gray(highlight),[1 1 3]),[result_path '/' name '.png']);
end
