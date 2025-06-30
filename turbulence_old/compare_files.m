close all; clear all;

mfc_dir = "/p/global/hyeoksu/MFC/reynolds/case/";
dir1 = mfc_dir + "p013-ml3-002/single-bubblescale-at-t0/restart_data/lustre_0.dat";
dir2 = mfc_dir + "p013-ml3-002/single-bubblescale-init/restart_data/lustre_0.dat";

% Parameters
mp = 1024; np = 1024; pp = 512; sys_size = 6;

% Read data
fileID = fopen(dir1,'r');
A = fread(fileID,'double');
fclose(fileID);

fileID = fopen(dir2,'r');
B = fread(fileID,'double');
fclose(fileID);

% Reassign density & velocity components
qcA = permute(reshape(A, mp, np, pp, sys_size),[4 1 2 3]);
qcB = permute(reshape(B, mp, np, pp, sys_size),[4 1 2 3]);

% Compare data
dq = qcA - qcB;
dq_abs = abs(dq);

