function [ water, fat ] = SSFPSeparationIdeal( img1, img2, img3, img4, img5, img6, te1,te2,te3 )
% SSFPSeparation separates supplied images into water and fat
% components.
%
% [water,fat] = SSFPSeparation(img1,img2,img3,img4,region) uses a
% SSFPDixon separation method to separate four images into water and
% fat components.
%
% img1   - First image
% img2   - Second image
% img3   - Third image
% img4   - Fourth image
% region - a boundary box which contains water in all images

% Init
rows = length(img1(:,1));
cols = length(img1(1,:));
s = {img1, img2, img3, img4};

% Select Top Signal Magnitude
[top_s, indices] = SelectTop3Pixels(s);



[water, fat] = IdealSeparationInOrderDebug( s{1}, s{2}, s{3}, s{4}, img5, img6,indices{1}, indices{2}, indices{3},te1,te2,te3 );


end

function [SelectPixels, SelectIndices] = SelectTop3Pixels(Images)

% Initialize
N = length(Images);
[Ny, Nx] = size(Images{1});
SelectPixels  = {zeros(Ny,Nx), zeros(Ny,Nx), zeros(Ny,Nx)};
SelectIndices = {zeros(Ny,Nx), zeros(Ny,Nx), zeros(Ny,Nx)};

% Select top 3 pixels from Phase Cycled Images
for r = 1:Ny
    for c = 1:Nx
        
        % Select top 3 pixels
        signal = zeros(1, N);
        for n = 1:length(Images)
            signal(n) = Images{n}(r,c);
        end
        
        [mag, indices] = InsertionSort(abs(signal));%abs(signal)
        
        % Save results
        for n = 1:3
            SelectPixels{n}(r,c) = signal(indices(N+1-n));
            SelectIndices{n}(r,c) = indices(N+1-n);
        end
        
    end
end
end

function [Array, Indices] = InsertionSort(Array)

Indices = 1:numel(Array);
for n = (2:numel(Array))
    
    value = Array(n);
    index = Indices(n);
    m = n - 1;
    
    while (m >= 1) && (Array(m) > value)
        Array(m+1) = Array(m);
        Indices(m+1) = Indices(m);
        m = m-1;
    end
    
    Array(m+1) = value;
    Indices(m+1) = index;
end
end

function [water, fat] = IdealSeparationInOrderDebug( Image1, Image2, Image3, Image4, Image5, Image6, Index1, Index2, Index3,te1,te2,te3 )
N = 3;
[Ny, Nx] = size(Image1);
water = zeros(Ny, Nx);
fat = zeros(Ny, Nx);
top1 = zeros(Ny, Nx);   % First In-Out
top2 = zeros(Ny, Nx);   % Second Out-In
top3 = zeros(Ny, Nx);   % left over
phase = zeros(Ny, Nx);
for n = 1:N
    for r = 1:Ny
        for c =1:Nx
            % Get Index Positions into Index Array
            Index = [0 0 0 0];
            Index = Index + (Index1(r,c) == 1) * [1 0 0 0] + (Index1(r,c) == 2) * [0 1 0 0] + (Index1(r,c) == 3) * [0 0 1 0] + (Index1(r,c) == 4) * [0 0 0 1];
            Index = Index + (Index2(r,c) == 1) * [1 0 0 0] + (Index2(r,c) == 2) * [0 1 0 0] + (Index2(r,c) == 3) * [0 0 1 0] + (Index2(r,c) == 4) * [0 0 0 1];
            Index = Index + (Index3(r,c) == 1) * [1 0 0 0] + (Index3(r,c) == 2) * [0 1 0 0] + (Index3(r,c) == 3) * [0 0 1 0] + (Index3(r,c) == 4) * [0 0 0 1];
            
            % Go through combinations
            if( Index(1) && Index(2) && Index(3) )
                if( (abs(Image1(r,c)) + abs(Image2(r,c))) > (abs(Image2(r,c)) + abs(Image3(r,c))) )
                    top1(r,c) = Image1(r,c);
                    top2(r,c) = Image2(r,c);
                    top3(r,c) = Image3(r,c);
                    phase(r,c) = angle(conj(Image5(r,c)).*Image2(r,c)) / 2;
                else
                    top1(r,c) = Image3(r,c);
                    top2(r,c) = Image2(r,c);
                    top3(r,c) = Image1(r,c);
                    phase(r,c) = angle(conj(Image5(r,c)).*Image2(r,c)) / 2;
                end
            elseif( Index(1) && Index(2) && Index(4))
                top1(r,c) = Image1(r,c);
                top2(r,c) = Image2(r,c);
                top3(r,c) = Image4(r,c);
                phase(r,c) = angle(conj(Image5(r,c)).*Image2(r,c)) / 2;
            elseif( Index(1) && Index(3) && Index(4))
                top1(r,c) = Image3(r,c);
                top2(r,c) = Image4(r,c);
                top3(r,c) = Image1(r,c);
                phase(r,c) = angle(conj(Image6(r,c)).*Image4(r,c)) / 2;
            elseif( Index(2) && Index(3) && Index(4))
                if( (abs(Image2(r,c)) + abs(Image3(r,c))) > (abs(Image3(r,c)) + abs(Image4(r,c))) )
                    top1(r,c) = Image3(r,c);
                    top2(r,c) = Image2(r,c);
                    top3(r,c) = Image4(r,c);
                    phase(r,c) = angle(conj(Image5(r,c)).*Image2(r,c)) / 2;
                else
                    top1(r,c) = Image3(r,c);
                    top2(r,c) = Image4(r,c);
                    top3(r,c) = Image2(r,c);
                    phase(r,c) = angle(conj(Image6(r,c)).*Image4(r,c)) / 2;
                end
            else
                'ideal Error!'
            end
            
        end
    end
end

img1 = zeros(Ny,Nx);
img2 = zeros(Ny,Nx);
img3 = zeros(Ny,Nx);
for r = 1:Ny
    for c = 1:Nx
        norm = abs(top1(r,c)) + abs(top2(r,c)) + abs(top3(r,c)) / 3;
        img1(r,c) = norm * exp(1i*angle(top1(r,c)));
        img2(r,c) = norm * exp(1i*angle(top2(r,c)));
        img3(r,c) = norm * exp(1i*angle(top3(r,c)));
    end
end

[water, fat] = IdealSSFP(img1,img2,img3,te1,te2,te3);
end
