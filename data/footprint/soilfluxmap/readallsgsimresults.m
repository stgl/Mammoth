% readallsgsimresults.m
%
% Matlab script to read all sgsim results and save them:

load grid_data

i = 1;

datamat = readsgsim(strcat('sgsim',num2str(i),'.out'),nx,ny,ns);
   
save(strcat('Simulation',num2str(i)),'datamat');
