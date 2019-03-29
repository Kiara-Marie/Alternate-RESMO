

function kion=kION(n,T) % ionizing collisions
    
    % output unit is um^3
    kB = 1.3806504e-23;% #m2 kg s-2 K-1;J K-1; 8.62e-5 #eV K-1,
    Rydhc = 2.179872e-18;% #J: Ryd [cm-1] --> [J]
    epsi=Rydhc./(power(n,2)*kB*T); % find reduced initial energy
    kion=11*sqrt(Rydhc./(kB*T)).*kNOT(T).*exp(-epsi)./(power(epsi,2.33)+4.38*power(epsi,1.72)+1.32*epsi);

end