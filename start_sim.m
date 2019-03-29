%This part solves the rate equations and saves results in workspace
for density=[0.08]
    tic
    N=20;%number of shells
    t_max=200;
    steps=100;


    sigma_z=0.42*1000;%Gaussian width
    sigma_y=0.42*1000;
    sigma_x=1.00*1000;
%     sigma_z=1*1000;%Gaussian width
%     sigma_y=0.55*1000;
%     sigma_x=0.70*1000;

    n=49; %PQN
    
   
    d_p=density; %peak density in um-3

    sigma_env=5;%consider number of sigma environments

    pos=linspace(0,sigma_env*sigma_z-0.5*sigma_env*sigma_z/(N-0.5),N);
    
    pos_x=linspace(0.5*sigma_env*sigma_x/(N-0.5),sigma_env*sigma_x,N)';
    pos_y=linspace(0.5*sigma_env*sigma_y/(N-0.5),sigma_env*sigma_y,N)';
    pos_z=linspace(0.5*sigma_env*sigma_z/(N-0.5),sigma_env*sigma_z,N)';
    
    
    d=arrayfun(@(z) d_p*exp(-(z^2)/(2*sigma_z^2)),pos);

    foldername=['R_PQN_',num2str(n),'othertry'];
    mkdir(foldername);
    filename=[foldername,'\','l=',num2str(N),'_shell_sim_n=',num2str(n),'_d0=', num2str(d_p), '_sigma_', num2str(sigma_z), 'mm_', ...
        num2str(sigma_x),'mm','_tfinal',num2str(t_max),'ns_',num2str(sigma_env),'sigmaenv'];
    
    %solve rate equations
    [time,nden,eden,deac_n_min,deac_dr,deac_pd,Te,rx,ry,rz,vx,vy,vz,vol,y0]=shell_rate_eqn_sim(d, pos_x, pos_y, pos_z, n, t_max/steps, t_max, false);
    %save workspace
    save(strcat([filename, '.mat']))
    
    %reduce file size by deleting shell specific information
    eden=EDEN(eden,vol);
    nden=NDEN(nden, vol);
    deac_n_min=NDEN(deac_n_min, vol);
    deac_pd=NDEN(deac_pd, vol);
    deac_dr=EDEN(deac_dr, vol);
    clearvars y0
    
    save(strcat([filename, '_small', '.mat']))
    
    toc
end


fprintf("Done everything\n");




