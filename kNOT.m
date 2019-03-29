function knot=kNOT(T)
    
    el = 1.6021765e-19;% #C
    kB = 1.3806504e-23;% #m2 kg s-2 K-1;J K-1; 8.62e-5 #eV K-1,
    epsilon = 8.854187817e-12;% #C2 J-1 m-1
    Rydhc = 2.179872e-18;% #J: Ryd [cm-1] --> [J]
    emass = 9.1093822e-31;% #kg
    
    knot=1e9*power(el,4.0)./(kB*T*sqrt(emass*Rydhc)*power(4*pi*epsilon,2)); % PRL 2008 PVS 
    % orginal units of knot = m^3/s
    % 1e9 converts it into um^3 / ns
    % dividing by (4\pi\epsilon)^2 converts a.u. --> s.i.
end