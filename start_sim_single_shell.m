%This part solves the rate equations and saves results in workspace
for density=0.04
    tic
    
    N=1000;%number of shells
    t_max=1000;
    steps=1000;


    sigma_z=0.42*1000;%Gaussian width
    sigma_y=0.42*1000;
    sigma_x=0.75*1000;
    n=49; %PQN
    
    d_p=density; %peak density in um-3
    
    d_eff=density;

    sigma_env=5;%consider number of sigma environments

    pos=linspace(0,sigma_env*sigma_z-0.5*sigma_env*sigma_z/(N-0.5),N);
    
    pos_x=linspace(0.5*sigma_env*sigma_x/(N-0.5),sigma_env*sigma_x,N)';
    pos_y=linspace(0.5*sigma_env*sigma_y/(N-0.5),sigma_env*sigma_y,N)';
    pos_z=linspace(0.5*sigma_env*sigma_z/(N-0.5),sigma_env*sigma_z,N)';
    
    
    d=arrayfun(@(z) d_p*exp(-(z^2)/(2*sigma_z^2)),pos);
    
    volume= 4/3*pi* (pos_x.*pos_y.*pos_z - [0; pos_x(1:end-1)].*[0; pos_y(1:end-1)].*[0; pos_z(1:end-1)]);
    
    [m,i]=max(volume.*d');
    d_p=d(i);
    d=d_p;
    
    pos_x=1;
    pos_y=1;
    pos_z=1;
    

    foldername=['Version2\Uncoupled_shells_PQN_',num2str(n),'_single_av'];
    mkdir(foldername);
    filename=[foldername,'\','l=',num2str(N),'_shell_sim_n=',num2str(n),'_d0=',...
        num2str(d_p),'_d_eff=',num2str(d_eff),'_single',...
        '_tfinal',num2str(t_max),'ns_'];
    
    %solve rate equations
    [time,nden,eden,deac_n_min,deac_dr,deac_pd,Te,rx,ry,rz,vx,vy,vz,vol,y0]=shell_rate_eqn_sim(d, pos_x, pos_y, pos_z, n, t_max/steps, t_max,true);
    %save workspace
    save(strcat([filename, '.mat']))
    
  
    
    toc
end


