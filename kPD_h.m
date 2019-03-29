function kpd=kPD_h(n) %helpfunction, only takes single value
    kpd=4.13e7/(2*pi*n^5)*(0.3054+(3e-5)*(n^2-16));
end