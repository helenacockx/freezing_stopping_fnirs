function DPF = age2DPF(age, wl)
% This function calculates the DPF based on age and wavelength(wl) as given by the formula of
% Scholkmann & Wolf, 2013 (https://doi.org/10.1117/1.JBO.18.10.105004)

a = 223.3;
b = 0.05624;
g = 0.8493;
d = -5.723 * 10^(-7);
e = 0.001245;
z = -0.9025;

DPF = a + b*age^g + d*wl^3 + e*wl^2 + z*wl;
