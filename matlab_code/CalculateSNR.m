function [SNR] = CalculateSNR(Signal, TissueModel, Tissue)
    Rows = size(TissueModel,1);
    Cols = size(TissueModel,2);
    
    s = [];
    n = [];
    for r = 1:Rows
        for c = 1:Cols
            if strcmp(Tissue, TissueModel(r,c).Name)
                s = [s abs(Signal(r,c))];
            elseif strcmp('Empty', TissueModel(r,c).Name)
                n = [n abs(Signal(r,c))];
            end
        end
    end
    
    mu = mean(s);
    sigmaN = std(n);
    
    SNR = abs(mu)/sigmaN;
end