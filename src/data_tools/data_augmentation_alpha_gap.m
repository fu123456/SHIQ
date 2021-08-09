% Data augmentation for our dataset using two operations including (a)
% specular highlight boosting; (b) specular highlight attenuation.
% This scheme can produce more diverse training images with wider
% range of highlight intensity, making our network more robust and
% general to real-world highlights

clear all;
close all;
clc;

% Please accordingly modify variables "DATA_DIR", "beta", "theta", "k" and
% "alpha". Note "output_dir" is automatically generated based on "DATA_DIR"
DATA_DIR='<you dir>'; % Input your dir for data

%% Parameters
beta=10;
theta=0.2;
% alpha=0.8;
k=0;

%% Boosting specular highlights
for alpha=1.5:0.5:5
    output_dir=[DATA_DIR num2str(alpha)];
    if ~exist(output_dir)
        mkdir(output_dir)
    end
    disp(output_dir)
    % _A: input; _S: specular; _D: diffuse
    dataDIR=DATA_DIR;
    dataFiles=dir(fullfile(dataDIR,'*_S.png'));
    for j=1:numel(dataFiles)
        [~,name,~]=fileparts(fullfile(dataDIR,dataFiles(j).name));
        disp(name);
        diffuse_name=strrep(name,'_S','_D'); % name of specular-free image
        input_name=strrep(name,'_S','_A'); % name of input image
        new_generate_name=strrep(name,'_S','_A'); % name of new generted image
        highlight=im2double(imread(fullfile(dataDIR,dataFiles(j).name)));
        if size(highlight,3)==3
            highlight=rgb2gray(highlight);
        end
        diffuse=im2double(imread(fullfile(dataDIR,[diffuse_name '.png'])));
        % Boost for highlight
        g=1./(1+exp(-beta.*(highlight-theta)));
        highlight_boost=highlight.*(k+alpha.*g);
        new=diffuse+repmat(highlight_boost,[1 1 3]);
        % Save results
        result_path=output_dir;
        imwrite(diffuse,[result_path '/' diffuse_name '.png']);
        imwrite(new,[result_path '/' new_generate_name '.png']);
        imwrite(highlight_boost,[result_path '/' name '.png']);
    end
end

%% Attenuating specular highlights
for alpha=0.2:0.2:0.8
    output_dir=[DATA_DIR num2str(alpha)];
    if ~exist(output_dir)
        mkdir(output_dir)
    end
    disp(output_dir)
    % _A: input; _S: specular; _D: diffuse
    dataDIR=DATA_DIR;
    dataFiles=dir(fullfile(dataDIR,'*_S.png'));
    for j=1:numel(dataFiles)
        [~,name,~]=fileparts(fullfile(dataDIR,dataFiles(j).name));
        disp(name);
        diffuse_name=strrep(name,'_S','_D'); % name of specular-free image
        input_name=strrep(name,'_S','_A'); % name of input image
        new_generate_name=strrep(name,'_S','_A'); % name of new generated image
        highlight=im2double(imread(fullfile(dataDIR,dataFiles(j).name)));
        if size(highlight,3)==3
            highlight=rgb2gray(highlight);
        end
        diffuse=im2double(imread(fullfile(dataDIR,[diffuse_name '.png'])));
        % Boost for highlight
        g=1./(1+exp(-beta.*(highlight-theta)));
        highlight_boost=highlight.*(k+alpha.*g);
        new=diffuse+repmat(highlight_boost,[1 1 3]);
        % Save results
        result_path=output_dir;
        imwrite(diffuse,[result_path '/' diffuse_name '.png']);
        imwrite(new,[result_path '/' new_generate_name '.png']);
        imwrite(highlight_boost,[result_path '/' name '.png']);
    end
end
