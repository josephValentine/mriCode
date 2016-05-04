%%-------------------------------------------------------------------------------------
% test bench for simulating compressed sensing with varying parameters.

clear;
%whitebg('k');%change background to black (it looks nice!)

% Inital Parameters
alpha = pi/3; phi = 0; dphi = pi; Nr = 200;%makes iterative return exact phase result of steadystate in SSFP.m
TR = 10/1000; TE = TR/2; f0 = 0;

FoSD             = 0.05;%.05       	% Off-Resonance Std. Dev. in Hz
MaxGradFreq      = 500; %50; 	% Max Frequency for Off-Resonance Linear Gradient
NoiseSD          = 0.004; %0.004;       % Signal Noise Std. Dev.

% Generate Tissue Model
addpath('/home/josephwv/Desktop/code/matlab_code/fatWaterSeparation')
addpath('/home/josephwv/Desktop/code/matlab_code/ssfp')
addpath('/home/josephwv/Desktop/code/matlab_code/tissueGeneration')
lookup_table = {'brain128.png', 'knee128.png', 'simpleFatMuscle.png'};
filename = lookup_table{2};
tissues = TissueElements();
tissue = TissueModel(filename, tissues, FoSD, MaxGradFreq);

phis = [0 0 0 0];
dphis = [0 1/4 2/4 3/4] * 2 * pi; 
df = abs(tissues.Fat.F0);
dt = 1 / (2 * df);
%TEs = [TE - 3*dt/2, TE - dt/2, TE + dt/2, TE + 3*dt/2];
TEs = [TE , TE + dt, TE , TE + dt];
rimg = SSFPImage(tissue, alpha, phis(1), dphis(1), TR, TEs(1), Nr, NoiseSD, f0);
rimg2 = SSFPImage(tissue, alpha, phis(2), dphis(2), TR, TEs(2), Nr, NoiseSD, f0);
rimg3 = SSFPImage(tissue, alpha, phis(3), dphis(3), TR, TEs(3), Nr, NoiseSD, f0);
rimg4 = SSFPImage(tissue, alpha, phis(4), dphis(4), TR, TEs(4), Nr, NoiseSD, f0);

theta = -[0 1/4 1/2 3/4] * pi;
[ img, img2, img3, img4 ] = AdjustSSFPConstantPhase(rimg, rimg4, rimg3, rimg2, theta);%makes steady states work

PlotMRImage(img, img2, img3, img4);
