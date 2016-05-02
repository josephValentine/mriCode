img1=zeros(50,50);%this represents object1
for i=5:45
    for j=5:45
        img1(i,j)=i^2;%coil oriented near top of object
    end
end

img2=zeros(50,50);%object 2
for i=5:45
    for j=5:45
        img2(i,j)=2025-i^2;%coil oriented near bottom of object
    end
end

sampled1=fftshift(fft2(fftshift(img1)));
sampled2=fftshift(fft2(fftshift(img2)));
temp=sampled1;
%now lets decimate some lines of our kspace data
for k=1:2:50
    if(k < 20 || k>30)
       sampled1(:,k)=0;
        sampled2(:,k)=0;
    end
end
figure
imshow(abs(ifftshift(ifft2(ifftshift(sampled1)))),[]);
figure
imshow(abs(ifftshift(ifft2(ifftshift(sampled2)))),[]);

%save a kernel (3*3)
kernel=sampled1(25:27,25:27);%should be the same for each sampled signal right?

%try to fix coil1

%first get coefficients


%maybe try with 9 eq even though you know some will be zeroed
Ma=zeros(9,6);
calStartX=25;
calStartY=20;%this doesnt matter since we decimated vertically
maRow=1;
for c=1:3%right shift
    for r=1:3%down shift
      % coil 1 
       xshift=r-1;
       yshift=c-1;
       
      Ma(maRow,1)=sampled1(calStartX+xshift,calStartY+yshift);
      Ma(maRow,2)=sampled1(calStartX+2+xshift,calStartY+yshift);
      Ma(maRow,3)=sampled1(calStartX+xshift,calStartY+1+yshift);
      Ma(maRow,4)=sampled1(calStartX+2+xshift,calStartY+1+yshift);
      Ma(maRow,5)=sampled1(calStartX+xshift,calStartY+2+yshift);
      Ma(maRow,6)=sampled1(calStartX+2+xshift,calStartY+2+yshift);
      %coil 2
      Ma(maRow,7)=sampled2(calStartX+xshift,calStartY+yshift);
      Ma(maRow,8)=sampled2(calStartX+2+xshift,calStartY+yshift);
      Ma(maRow,9)=sampled2(calStartX+xshift,calStartY+1+yshift);
      Ma(maRow,10)=sampled2(calStartX+2+xshift,calStartY+1+yshift);
      Ma(maRow,11)=sampled2(calStartX+xshift,calStartY+2+yshift);
      Ma(maRow,12)=sampled2(calStartX+2+xshift,calStartY+2+yshift);

      maRow=maRow+1;
    end

end
    
      
Mk=[kernel(1,1)
    kernel(2,1)
    kernel(3,1)
    kernel(1,2)
    kernel(2,2)
    kernel(3,2)
    kernel(1,3)
    kernel(2,3)
    kernel(3,3)];%this is almost certanly wrong. i dont really get the whole kernel thing yet

%ak=inv(Ma' * Ma)* Ma' * Mk;%solve the system
ak=Ma\Mk;
%ak = (conj(Ma)*Ma + eye(9,9)) \ (conj(Ma)*Mk);

%fill in the missing rows on left
%take a 3by3 region and move down each row
fixed1=sampled1;
for x=3:2:19%move to each missing col (skipping first for now)
    for y=2:49%go down, filling each point. (skipping first for now)
       
        coil1=ak(1)*sampled1(x-1,y-1)+ ak(2)*sampled1(x+1,y-1)+ ak(3)*sampled1(x-1,y)+ ak(4)sampled1*(x+1,y)+ak(5)sampled1*(x-1,y+1)+ ak(6)*sampled1(x+1,y+1);
        coil2=ak(7)*sampled1(x-1,y-1)+ ak(8)*sampled1(x+1,y-1)+ ak(9)*sampled1(x-1,y)+ ak(10)sampled1*(x+1,y)+ak(11)sampled1*(x-1,y+1)+ ak(12)*sampled1(x+1,y+1);
        fixed1(x,y)=ak*[coil1 coil2];
    
    end
end
%fill right half
for x=31:2:50%move to each missing col (skipping first for now)
    for y=2:49%go down, filling each point. (skipping first for now)
        
        coil1=[sampled1(x-1,y-1), sampled1(x+1,y-1), sampled1(x-1,y), sampled1(x+1,y),sampled1(x-1,y+1), sampled1(x+1,y+1)];
        coil2=[sampled2(x-1,y-1), sampled2(x+1,y-1), sampled2(x-1,y), sampled2(x+1,y),sampled2(x-1,y+1), sampled2(x+1,y+1)];
        fixed1(x,y)=ak.*coil1
    
    
    end
end

figure(1);
subplot(2,2,1)
imshow(abs(ifftshift(ifft2(ifftshift(temp)))),[]);%transform from kspace,[]);
xlabel('coil1 image space');

subplot(2,2,2);
imshow(angle(temp),[]);
xlabel('coil1 k space');

subplot(2,2,3)
imshow(abs(ifftshift(ifft2(ifftshift(sampled1)))),[]);
xlabel('coil2 image space');

subplot(2,2,4);
imshow(angle(sampled1),[]);
xlabel('coil2 k space');

figure(2);
subplot(2,2,1)
imshow(abs(ifftshift(ifft2(ifftshift(fixed1)))),[]);%transform from kspace,[]);
xlabel('coil1 image space');

subplot(2,2,2);
imshow(angle(fixed1),[]);
xlabel('coil1 k space');

subplot(2,2,3)
imshow(abs(ifftshift(ifft2(ifftshift(sampled2)))),[]);
xlabel('coil2 image space');

subplot(2,2,4);
imshow(angle(sampled2),[]);
xlabel('coil2 k space');

        
        