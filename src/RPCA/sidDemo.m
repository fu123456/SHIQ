% Xiaojie Guo, Oct 2013. 
% Questions? xj.max.guo@gmail.com
%
% Reference: Robust Separation of Reflection Using Multiple Images,
% Xiaojie Guo, Xiaochun Cao and Yi Ma, CVPR 2014

close all;clear all;clc;

% Load data
load('./data/dataFace.mat');

tmpT = transfm;
tmpI = I0R;
canonSize = floor([max(ROI{1}(2,:))-min(ROI{1}(2,:)),...
             max(ROI{1}(1,:))-min(ROI{1}(1,:)),]+1);
         
% For R channel
[Fotr, Tr, Rr, Nr, tran] = SID(tmpI, tmpT, canonSize);
tmpT = tran;

% For G
mode = 0;
tmpI = I0G;
[Fotg, Tg, Rg, Ng] = SID(tmpI, tmpT, canonSize, mode);

% For B
tmpI = I0B;
[Fotb, Tb, Rb, Nb] = SID(tmpI, tmpT, canonSize, mode);
close all
%% Show results
for i = 1 : size(Fotr,2)
   figure;hold on;
   subplot(1,3,1);
   tmp = cat(3, reshape(Fotr(:,i),canonSize),...
       reshape(Fotg(:,i),canonSize),reshape(Fotb(:,i),canonSize));
   imshow(tmp,[]);
   
   subplot(1,3,2);
   tmp = cat(3, reshape(Tr(:,i),canonSize),...
       reshape(Tg(:,i),canonSize),reshape(Tb(:,i),canonSize));
   imshow(tmp,[]);
   subplot(1,3,3);
   tmp = cat(3, reshape(Rr(:,i),canonSize),...
       reshape(Rg(:,i),canonSize),reshape(Rb(:,i),canonSize));
   imshow(tmp,[])
   hold off
end