% Kelvin Zhang, Arvind Ganesh, February 2011. 
% Questions? zhangzdfaint@gmail.com, abalasu2@illinois.edu
%
% Copyright: Perception and Decision Laboratory, University of Illinois, Urbana-Champaign
%            Microsoft Research Asia, Beijing
%
% Reference: TILT: Transform Invariant Low-rank Textures  
%            Zhengdong Zhang, Xiao Liang, Arvind Ganesh, and Yi Ma. Proc. of ACCV, 2010.
%

function H=compute_homography(X, U)
% compute_homography(X, U) will compute the homography from pts X to pts U.
% ---------------------------input----------------------------------------
% X:        2*4 matrix.
% U:        2*4 matrix.
% ---------------------------output---------------------------------------
% H:        3*3 matrix.
L=zeros(8, 8);
b=zeros(8, 1);
for i=1:4
    L(2*i-1, :)=[0 0 0 -X(1, i) -X(2, i) -1 U(2, i)*X(1, i) U(2, i)*X(2, i)];
    b(2*i-1)=-U(2, i);
    L(2*i, :)=[X(1, i) X(2, i) 1 0 0 0 -U(1, i)*X(1, i) -U(1, i)*X(2, i)];
    b(2*i)=U(1, i);
end
H=reshape([L\b;1], 3, 3)';

