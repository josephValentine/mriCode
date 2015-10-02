
% Algorithm Bench
clear;
whitebg('k');
% Inital Parameters
alpha = pi/8; phi = 0; dphi = pi; Nr = 200;%makes iterative return exact phase result of steadystate
TR = 10/1000; TE = TR/2; f0 = 0;

FoSD             = 0.05;%.05       	% Off-Resonance Std. Dev. in Hz
MaxGradFreq      = 500; %50; 	% Max Frequency for Off-Resonance Linear Gradient
NoiseSD          = 0.004; %0.004;       % Signal Noise Std. Dev.

% Generate Tissue Model
lookup_table = {'TissueMasks/PhantomDeathStar.png', 'TissueMasks/Phantom64Tissues.png', 'TissueMasks/MuscleCartFat64.png'};
filename = lookup_table{1};
tissues = TissueElements();
tissue = TissueModel(filename, tissues, FoSD, MaxGradFreq);

%% ------------------------------------------------------------------------
% Plot Spectrum of 4 Phase Cycled Images
if 0

    phis = -[0 1/4 2/4 3/4] * pi;
    dphis = [0 1/4 2/4 3/4] * 2 * pi; 
    df = abs(tissues.Fat.F0);
    dt = 1 / (2 * df);
    TEs = [TE , TE + dt, TE , TE + dt];

    PlotTissueSpectrum(alpha, phis(1), dphis(1), TEs(1), TR);
    PlotTissueSpectrum(alpha, phis(2), dphis(2), TEs(2), TR);
    PlotTissueSpectrum(alpha, phis(3), dphis(3), TEs(3), TR);
    PlotTissueSpectrum(alpha, phis(4), dphis(4), TEs(4), TR);
end

%% ------------------------------------------------------------------------
% Generate 4 Phase Cycled Images
if 0
    %phis = -[0 1/4 2/4 3/4] * pi;
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
    %PlotSpectrum(img,20);
    
    data = {img,img2,img3,img4};
    save('sim_data.mat','data');
end

%% ------------------------------------------------------------------------
% Generate 6 Phase Cycled Images to compare dixon and ideal
if 0
    %phis = -[0 1/4 2/4 3/4] * pi;
    phis = [0 0 0 0];
    dphis = [0 0 0 0] * 2 * pi; 
    df = abs(tissues.Fat.F0);
    dt = 1 / (2 * df);
%     TEs = [TE - 3*dt/2, TE - dt/2, TE + dt/2, TE + 3*dt/2, TE - dt/2 - 2*dt, TE + 3*dt/2 - 2*dt];
    TEs = [TE, TE + dt, TE + 2*dt,TE + dt/6, TE + 5*dt/6, TE + 3*dt/2];

    img = SSFPImage(tissue, alpha, phis(1), dphis(1), TR, TEs(1), Nr, NoiseSD, f0);
    img2 = SSFPImage(tissue, alpha, phis(2), dphis(2), TR, TEs(2), Nr, NoiseSD, f0);
    img3 = SSFPImage(tissue, alpha, phis(3), dphis(3), TR, TEs(3), Nr, NoiseSD, f0);
    img4 = SSFPImage(tissue, alpha, phis(4), dphis(4), TR, TEs(4), Nr, NoiseSD, f0);
    img5 = [];
    img6 = SSFPImage(tissue, alpha, phis(2), dphis(2), TR, TEs(5), Nr, NoiseSD, f0);
    img7 = [];
    img8 = SSFPImage(tissue, alpha, phis(1), dphis(1), TR, TEs(6), Nr, NoiseSD, f0);%inphase
    
    theta = -[0 0 0 0 0 0] *pi;%was neg and 2pi
     
    img = img * exp(-1i*theta(1));
    img2 = img2 * exp(-1i*theta(2));
    img3 = img3 * exp(-1i*theta(3));
    img4 = img4 * exp(-1i*theta(4));
    img6 = img6 * exp(-1i*theta(5));
    img8 = img8 * exp(-1i*theta(1));%fix
   
    
    %PlotMRImage(img, img2, img3, img4, img6, img8);
%     PlotSpectrum(img,27);
%     PlotSpectrum(img,39);
%     PlotSpectrum(img2,27);
%     PlotSpectrum(img2,39);
%     PlotSpectrum(img3,27);
%     PlotSpectrum(img3,39);
%     PlotSpectrum(img4,27);
%     PlotSpectrum(img4,39);
    
end

%% ------------------------------------------------------------------------
% Generate 6 Phase Cycled Images for ssfp separation
if 1
    
    phis = [0 0 0 0];
    dphis = [0 1/4 1/2 3/4] * 2 * pi; 
    df = abs(tissues.Fat.F0);
    dt = 1 / (2 * df);
    TEs = [TE , TE + dt, TE, TE + dt, TE - dt, TE - dt];
   
    img = SSFPImage(tissue, alpha, phis(1), dphis(1), TR, TEs(1), Nr, NoiseSD, f0);
    img2 = SSFPImage(tissue, alpha, phis(2), dphis(2), TR, TEs(2), Nr, NoiseSD, f0);
    img3 = SSFPImage(tissue, alpha, phis(3), dphis(3), TR, TEs(3), Nr, NoiseSD, f0);
    img4 = SSFPImage(tissue, alpha, phis(4), dphis(4), TR, TEs(4), Nr, NoiseSD, f0);
    img5 = [];
    img6 = SSFPImage(tissue, alpha, phis(2), dphis(2), TR, TEs(5), Nr, NoiseSD, f0);
    img7 = [];
    img8 = SSFPImage(tissue, alpha, phis(1), dphis(4), TR, TEs(6), Nr, NoiseSD, f0);
    
    theta = [0 1/4 1/2 3/4 1/4 3/4] *pi;
     
    img = img * exp(-1i*theta(1));
    img2 = img2 * exp(-1i*theta(2));
    img3 = img3 * exp(-1i*theta(3));
    img4 = img4 * exp(-1i*theta(4));
    img6 = img6 * exp(-1i*theta(5));
    img8 = img8 * exp(-1i*theta(6));
   

    
    PlotMRImage(img);
    %PlotMRImage(img, img2, img3, img4, img6, img8);
%     PlotSpectrum(img,27);
%     PlotSpectrum(img2,27);
%     PlotSpectrum(img3,27);
%     PlotSpectrum(img4,27);
    
end

%% ------------------------------------------------------------------------
% Plot phase gradient
if 0
    [fx, fy] = gradient(angle(img));
    [fx2, fy2] = gradient(angle(img2));
    [fx3, fy3] = gradient(angle(img3));
    [fx4, fy4] = gradient(angle(img4));
  grad = sqrt(fx.^2 + fy.^2);
  grad2 = sqrt(fx2.^2 + fy2.^2);
  grad3 = sqrt(fx3.^2 + fy3.^2);
  grad4 = sqrt(fx4.^2 + fy4.^2);
  
  [valmap] = CalculateMinDistToNull(abs(fx),abs(fy));
  invmap = ones(64,64);
  figure();
  
  surf(valmap);
  figure();
  surf(abs(img));
  %figure();
  %surf((invmap - grad));
%     data = grad(:,27);%col num
%     figure();
%     plot(grad(:,27));
%     figure();
%     plot(grad2(:,27));
%     figure();
%     plot(grad3(:,27));
%     figure();
%     plot(grad4(:,27));
   
end
%% ------------------------------------------------------------------------
% test to see where swap error occurs 
if 0
    range = 1:64;
    data = img(:,27);
    data2 = img2(:,27);
    data3 = img3(:,27);
    data4 = img4(:,27);
    figure();
    subplot(2,1,1); plot(range,abs(data),range,abs(data2),':',range,abs(data3),':',range,abs(data4),':');
    subplot(2,1,2); plot(range,angle(data),range,angle(data2),range,angle(data3),range,angle(data4));
end

%% ------------------------------------------------------------------------
% Separate water and fat images with Dixon Separation
if 0
    [ water, fat ] = DixonSeparation(img, img2);
    PlotMRImage(water, fat);
end

%% ------------------------------------------------------------------------
% compare 3pt and ideal
if 0
    [water,fat] = DixonSeparation(img,img2);
    [ water2, fat2 ] = DixonSeparation(img, img2, img3);
    %[ water3, fat3 ] = IdealSSFP(img, img2, img3,TE, TE + dt, TE + 2*dt);
    [ water3, fat3 ] = IdealSSFP(img4, img6, img8,TE + dt/6, TE + 5*dt/6, TE + 3*dt/2);
PlotMRImage(img,img2,img3);
    
    %     PlotMRImage(water, fat);
%     PlotMRImage(water2, fat2);
%     PlotMRImage(water3,fat3);
%     
    %[CNR] = CalculateCNR(Signal, TissueModel, Tissue1, Tissue2)
    %[CNR] = CalculateCNR(Signal, TissueModel, Tissue1, Tissue2)
    %[CNR] = CalculateCNR(Signal, TissueModel, Tissue1, Tissue2)
    
    [SNRfatDixon] = CalculateSNR(fat2, tissue, 'Fat')
    [SNRfatIdeal] = CalculateSNR(fat3, tissue, 'Fat')
    
    [SNRcartDixon] = CalculateSNR(water2, tissue, 'Cartilage')
    [SNRcartIdeal] = CalculateSNR(water3, tissue, 'Cartilage')
    
    [SNRmusDixon] = CalculateSNR(water2, tissue, 'Muscle')
    [SNRmustIdeal] = CalculateSNR(water3, tissue, 'Muscle')
    
    
end

%% ------------------------------------------------------------------------
% Separate water and fat images with Dixon Separation with 4 Images all combinations
if 0
    [ water, fat ] = DixonSeparation(img, img2);
    [ water2, fat2 ] = DixonSeparation(img2, img3);
    [ water3, fat3 ] = DixonSeparation(img3, img4);
    [ water4, fat4 ] = DixonSeparation(img4, img);
    
    PlotMRImage(water, fat);
    PlotMRImage(water2, fat2);
    PlotMRImage(water3, fat3);
    PlotMRImage(water4, fat4);    
end

%% ------------------------------------------------------------------------
% Separate water and fat images with Dixon Separation with 5 Images all combinations
if 0    
    [ water, fat ] = DixonSeparation(img, img2);
    [ water2, fat2 ] = DixonSeparation(img2, img3);
    [ water3, fat3 ] = DixonSeparation(img3, img4);
    [ water4, fat4 ] = DixonSeparation(img4, img5);
    
    PlotMRImage(water, fat);
    PlotMRImage(water2, fat2);
    PlotMRImage(water3, fat3);
    PlotMRImage(water4, fat4);
end

%% ------------------------------------------------------------------------
% Separate water and fat images with SSFP Separation
if 0
    [water, fat ] = SSFPSeparation(img, img2, img3, img4, 1);
    PlotMRImage(water, fat);
    
    figure();
    subplot(1,6,1); imshow(abs(img),[]);
    subplot(1,6,2); imshow(abs(img2),[]);
    subplot(1,6,3); imshow(abs(img3),[]);
    subplot(1,6,4); imshow(abs(img4),[]);
    subplot(1,6,5); imshow(abs(water),[]);
    subplot(1,6,6); imshow(abs(fat),[]);
    
    
%     [waterIwaterEst, waterIfatE, waterIsigmaN] = CalculateTissueSignalAverage(water, tissue, tissues.WhiteMatter, tissues.Fat);
%     [fatIwaterE, fatIfatE, fatIsigmaN] = CalculateTissueSignalAverage(fat, tissue, tissues.WhiteMatter, tissues.Fat);
%     [waterIwaterEst, waterIfatE, waterIsigmaN fatIwaterE, fatIfatE, fatIsigmaN]
%     
%     [SNRW, sigmaN] = CalculateSNR(water, tissue, tissues.WhiteMatter);
%     SNRF = CalculateSNR(fat, tissue, tissues.Fat);
%     CNRW = CalculateCNR(water, tissue, tissues.WhiteMatter, tissues.Fat);
%     CNRF = CalculateCNR(fat, tissue, tissues.WhiteMatter, tissues.Fat);
%     [SNRW SNRF CNRW CNRF sigmaN]
end

%% ------------------------------------------------------------------------
% Separate water and fat images with SSFP Separation
if 1
    [water, fat ] = SSFPSeparationWithFieldMap(img, img2, img3, img4, img5, img6, img7, img8, 1);
    PlotMRImage(water, fat);
    
    %[jimg1,jimg2,jimg3] = SSFPGradientRecon(img,img2,img3,img4,img6,img8);
    %[ water2, fat2 ] = DixonSeparation(jimg1, jimg2,jimg3);
    %PlotMRImage(jimg1, jimg2, jimg3);
    [water2,fat2]=SSFPGradientReconV2(img,img2,img3,img4,img6,img8,'gradient');
    PlotMRImage(water2,fat2);
    
    %[water3, fat3 ] = SSFPGradientReconM(img, img2, img3, img4, img5, img6, img7, img8, 1);
    %PlotMRImage(water3, fat3);
%     figure();
%     subplot(2,6,1); imshow(abs(img),[]);
%     subplot(2,6,2); imshow(abs(img2),[]);
%     subplot(2,6,3); imshow(abs(img3),[]);
%     subplot(2,6,4); imshow(abs(img4),[]);
%     subplot(2,6,5); imshow(abs(water),[]);
%     subplot(2,6,6); imshow(abs(fat),[]);
%     subplot(2,6,7); imshow(angle(img),[]);
%     subplot(2,6,8); imshow(angle(img2),[]);
%     subplot(2,6,9); imshow(angle(img3),[]);
%     subplot(2,6,10); imshow(angle(img4),[]);
%     subplot(2,6,11); imshow(angle(water),[]);
%     subplot(2,6,12); imshow(angle(fat),[]);
    
    %figure();
    %subplot(1,2,1); imshow(angle(conj(img6).*img2) / 2, []);
    %subplot(1,2,2); imshow(angle(conj(img8).*img4) / 2, []);
    
%     [waterIwaterEst, waterIfatE, waterIsigmaN] = CalculateTissueSignalAverage(water, tissue, tissues.WhiteMatter, tissues.Fat);
%     [fatIwaterE, fatIfatE, fatIsigmaN] = CalculateTissueSignalAverage(fat, tissue, tissues.WhiteMatter, tissues.Fat);
%     [waterIwaterEst, waterIfatE, waterIsigmaN fatIwaterE, fatIfatE, fatIsigmaN]
%     
%     [SNRW, sigmaN] = CalculateSNR(water, tissue, tissues.WhiteMatter);
%     SNRF = CalculateSNR(fat, tissue, tissues.Fat);
%     CNRW = CalculateCNR(water, tissue, tissues.WhiteMatter, tissues.Fat);
%     CNRF = CalculateCNR(fat, tissue, tissues.WhiteMatter, tissues.Fat);
%     [SNRW SNRF CNRW CNRF sigmaN]
end

%% ------------------------------------------------------------------------
% Separate water and fat images with SSFP Separation with 5 Images
if 0
	[ water, fat ] = SSFPSeparationWith5Images(img, img2, img3, img4, img5, 0);
    PlotMRImage(water, fat);
    
    [waterIwaterEst, waterIfatE, waterIsigmaN] = CalculateTissueSignalAverage(water, tissue, tissues.WhiteMatter, tissues.Fat);
    [fatIwaterE, fatIfatE, fatIsigmaN] = CalculateTissueSignalAverage(fat, tissue, tissues.WhiteMatter, tissues.Fat);
    [waterIwaterEst, waterIfatE, waterIsigmaN fatIwaterE, fatIfatE, fatIsigmaN]
    
%     [SNRW, sigmaN] = CalculateSNR(water, tissue, tissues.WhiteMatter);
%     SNRF = CalculateSNR(fat, tissue, tissues.Fat);
%     CNRW = CalculateCNR(water, tissue, tissues.WhiteMatter, tissues.Fat);
%     CNRF = CalculateCNR(fat, tissue, tissues.WhiteMatter, tissues.Fat);
%     [SNRW SNRF CNRW CNRF sigmaN]
end

%% ------------------------------------------------------------------------
% Noise Analysis
if 0
    NoiseAnalysis();
end

