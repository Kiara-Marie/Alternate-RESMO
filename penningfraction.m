%%%%%%%%%%%%%%%%%%%%
% Penning fraction %
%%%%%%%%%%%%%%%%%%%%

    % This function calculates the inital seed electrons from Penning fraction

    % n are initial PQN-levels before Penning ionization
    % den is Ry-atom density in 10^12 pcc

    % eden is the # of electrons produced
    % rden are the remaining PQN-levels after Penning ionization

function [PenningFraction eden rden]=penningfraction(n,den)     % den in 10^12 pcc

    a0=5.2917721092e-5;     % bohr radius in um
    Rn0=n.^2*a0;            % radius of Rydb. atom by bohr model using semi-classical method
    Rmax=1.8*(Rn0*2);       % Robicheaux paper, within this distance, 90% penning ionize (~10ns)

    %PenningFraction=zeros(length(n),1);
    PenningFraction_analytic=zeros(length(n),1);

    % Calculates number of Penning partners based on Erlang distribution:
    for i=1:length(n)
    %PenningFraction(i)=quadv(@(r)4*pi*den*r.^2.*exp(-4*pi*den*r.^3/3),0,Rmax(i),1e-10); % proportion between 0 and Rmax
    PenningFraction_analytic(i)=1-exp(-4*pi*den*Rmax.^3/3); % integral is solved analytically
    end

    PenningFraction = PenningFraction_analytic;         % 90% ionization within certain time (assume rest non-interactive)
    eden=PenningFraction/2*den*.9;      % the dens. of electrons produce is half the proportion (1e- per partner)
    rden=(1-PenningFraction*0.9)*den;    % this is remaining density of rydbergs

end