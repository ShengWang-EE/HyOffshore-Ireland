function [windTurbine] = windTurbinePara()
%WINDTURBINEPARA Summary of this function goes here
%   Detailed explanation goes here
windTurbine.ratedPower = 2; windTurbine.cutinSpeed = 3; windTurbine.ratedSpeed = 11; windTurbine.cutoutSpeed = 22;
windTurbine.hubHeight = 100; windTurbine.rotorRadius = 50; windTurbine.roughness = 0.0002;
windTurbine.C_T = 0.88 / 2; windTurbine.airDensity = 1.25;
end

