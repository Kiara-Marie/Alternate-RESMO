function [time,nden,eden,deac_n_min,deac_dr,deac_pd,Te,rx,ry,rz,vx,vy,vz,vol,y0] = start_sim_fun(density,N,n,sigma_z,sigma_x,t_max,steps,dirname,name,noHydro)
%START_SIM_FUN This part solves the rate equations
        %     density = %peak density in um-3 (usually around 0.04)
        %     N = number of shells
        %     sigma_z=1*1000; %Gaussian width
        %     sigma_x=0.70*1000;
        %     t_max = max time in ns
        %     steps = how many time steps
        %     n = initial PQN
        
    tic
   
    sigma_y=sigma_z;
    
    d_p=density; 

    sigma_env=5;%consider number of sigma environments

    pos=linspace(0,sigma_env*sigma_z-0.5*sigma_env*sigma_z/(N-0.5),N);
    
    pos_x=linspace(0.5*sigma_env*sigma_x/(N-0.5),sigma_env*sigma_x,N)';
    pos_y=linspace(0.5*sigma_env*sigma_y/(N-0.5),sigma_env*sigma_y,N)';
    pos_z=linspace(0.5*sigma_env*sigma_z/(N-0.5),sigma_env*sigma_z,N)';
    
    
    d=arrayfun(@(z) d_p*exp(-(z^2)/(2*sigma_z^2)),pos);

    mkdir(dirname);
   % filename= [dirname , '\' , name];
    if noHydro
        single = noHydro;
    else
        single = false;
    end
    %solve rate equations
    [time,nden,eden,deac_n_min,deac_dr,deac_pd,Te,rx,ry,rz,vx,vy,vz,vol,y0]=shell_rate_eqn_sim(d, pos_x, pos_y, pos_z, n, steps, t_max, single);
    %save workspace
   % save(filename +"_den0_"+string(nden(1,1))+  ".mat");
  % save("withTSpan.mat");
    
    toc
end

