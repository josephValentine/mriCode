function [minDist] = CalculateMinDistToNull(Fx,Fy)
% return an int matrix representing the number of pixels to the nearest
% likely signal null for each pixel
% this is done by setting a threshold that most likely represents a gradient value near a phase shift
% (usually about 1-1.5) and traveling in four directions out from each point untill we reach that threshold

%this can be improved by letting it travel in more that just our four straight cardinal lines
Fx = abs(Fx);
Fy = abs(Fy);

rows = size(Fx,1);
cols = size(Fx,2);
minDist = zeros(rows,cols);
count=[];
keepGoing = false;

for r = 1:rows
    for c = 1:cols
        %go up
        rowtemp = r;
        count(1) = 0;
        if(rowtemp > 1 && Fy(rowtemp,c) < .6)
            keepGoing = true;
        end
        while(keepGoing)
            rowtemp = rowtemp - 1;
            if(rowtemp <= 1 || Fy(rowtemp,c) > .6)
                keepGoing = false;
            end
            
            count(1)=count(1)+1;
            
        end
        
        %go down
        rowtemp = r;
        count(2) = 0;
        if(rowtemp < rows && Fy(rowtemp,c) < .6)
            keepGoing = true;
        end
        while(keepGoing)
           rowtemp = rowtemp + 1;
            if(rowtemp >= rows || Fy(rowtemp,c) > .6)
                keepGoing = false;
            end
            
            count(2)=count(2)+1;
            
        end
        
        %go left
        coltemp = c;
       count(3) = 0;
        if(coltemp > 1 && Fx(r,coltemp) < .6)
            keepGoing = true;
        end
        while(keepGoing)
            coltemp = coltemp - 1;
            if(coltemp <= 1 || Fx(r,coltemp) > .6)
                keepGoing = false;
            end
            
            count(3)=count(3)+1;
           
        end
        
        %go right
        coltemp = c;
        count(4) = 0;
        if(coltemp < cols && Fx(r,coltemp) < .6)
            keepGoing = true;
        end
        while(keepGoing)
            coltemp = coltemp + 1;
            if(coltemp >= cols || Fx(r,coltemp) > .6)
                keepGoing = false;
            end
            
            count(4)=count(4)+1;
            
        end
        
        minDist(r,c) = min(count);
    end
end


