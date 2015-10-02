function PlotSpectrum(image,colnum)
    data = image(:,colnum);
    figure();
    subplot(2,1,1); plot(abs(data));
    subplot(2,1,2); plot(angle(data));
    
       
    

end


