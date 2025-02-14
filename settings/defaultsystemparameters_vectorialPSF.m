function [sys,opt,out]=defaultsystemparameters_vectorialPSF
%units are in m

sys = [];
out = [];    
opt.Et = 0; % Display Et during the calculation.
opt.Ef = 0;    % Display Ef during the calculation.
opt.calbar = 0;    % Display calculation progress bar.
opt.mem = 50000;
opt.pixSize = 10e-9; % A pixel size of 5e-9 is recommended for fine calculation.
opt.radiusCanvas = 0.7e-6; % Lateral range.
opt.depthCanvas = 0.8e-6; % Axial range. 
opt.polAngle = 0.0*pi; % Linear polarization angle.
opt.phaseImage = false; % Show the phase image at the back focal plane.
opt.intImage = false;

% Setting specific parameters
% [sys,out]=effInit_oil_exc(sys,out,opt);        % Assigning initial parameters. 
sys.NA=1.35;
sys.nm=1.406;
sys.ns=1.406;
sys.Mt=100;
sys.wa=7e-3;
sys.Na=200; 
sys.phaseImage = opt.phaseImage;
sys.intImage = opt.intImage;
sys.wa=7e-3;
sys.rz = 0.8e-6;% Axial range.
sys.pl = 0; % Angle of the linear polarization.
sys.beadradius = 0; %convolution with bead
sys.maskshift = [0 0];