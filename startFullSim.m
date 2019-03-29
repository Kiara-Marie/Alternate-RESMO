function [time,nden,eden,deac_n_min,deac_dr,deac_pd,Te,rx,ry,rz,vx,vy,vz,vol,y0] = startFullSim(density,n,N,t_max,steps,NP,dirname,name)
addToPath = pwd;
addpath(addToPath);
%STARTFULLSIM Run this to see the bifurcation model                                           

        %     density = %peak density in um-3 (usually around 0.04)
        %     N = number of shells
        %     n = initial PQN
        %     sigma_z=1*1000; %Gaussian width
        %     sigma_x=0.70*1000;
        %     t_max = max time in ns
        %     steps = how many time steps
        %     NP = number of particles (for particle model)
              sigma_z=1*1000; %Gaussian width in um
             sigma_x=0.70*1000; %um
states = 100;

[time,nden,eden,deac_n_min,deac_dr,deac_pd,Te,rx,ry,rz,vx,vy,vz,vol,y0] = start_sim_fun(density,N,n,sigma_z,sigma_x,t_max,steps,dirname,name);
%particlesFromShells(rx,ry,rz,vx,vy,vz,nden,eden,NP, states, t_max,steps);


end

