function [water,fat] = SSFPGradientRecon(in1, out1, in2, out2, fix1, fix2, mode)
%function [inPhase,outPhase1,fixPhase] = SSFPGradientRecon(in1, out1, in2, out2, fix1, fix2)
rows = length(in1(:,1));
cols = length(in1(1,:));

inPhase = zeros(rows,cols);
outPhase1 = zeros(rows,cols);
fixPhase = zeros(rows,cols);
water = zeros(rows,cols);
fat = zeros(rows,cols);
% correct mag for all images
for r=1:rows
    for c=1:cols
        s(1)=abs(in1(r,c));
        s(2)=abs(out1(r,c));
        s(3)=abs(in2(r,c));
        s(4)=abs(out2(r,c));
        
        [val,ind] = min(s);
        
        if ind ==1
            inPhase(r,c) = (s(2)+s(3)+s(4))/3;
            outPhase1(r,c) = (s(2)+s(3)+s(4))/3;
            %outPhase2(r,c) = (s(2)+s(3)+s(4))/3;
        elseif ind ==2
            inPhase(r,c) = (s(1)+s(3)+s(4))/3;
            outPhase1(r,c) = (s(1)+s(3)+s(4))/3;
            %outPhase2(r,c) = (s(1)+s(3)+s(4))/3;
        elseif ind ==3
            inPhase(r,c) = (s(2)+s(1)+s(4))/3;
            outPhase1(r,c) = (s(2)+s(1)+s(4))/3;
            %outPhase2(r,c) = (s(2)+s(1)+s(4))/3;
        else
            inPhase(r,c) = (s(2)+s(3)+s(1))/3;
            outPhase1(r,c) = (s(2)+s(3)+s(1))/3;
            %outPhase2(r,c) = (s(2)+s(3)+s(1))/3;
        end
    end
end

%find phase for inphase image
[fxin1,fyin1]=gradient(angle(in1));
[fxin2,fyin2]=gradient(angle(in2));
[fxout1,fyout1]=gradient(angle(out1));
[fxout2,fyout2]=gradient(angle(out2));

% [fxin1,fyin1]=gradient(abs(in1));
% [fxin2,fyin2]=gradient(abs(in2));
% [fxout1,fyout1]=gradient(abs(out1));
% [fxout2,fyout2]=gradient(abs(out2));

%new
[valmap1] = CalculateMinDistToNull(fxin1,fyin1);
[valmap2] = CalculateMinDistToNull(fxout1,fyout1);
[valmap3] = CalculateMinDistToNull(fxin2,fyin2);
[valmap4] = CalculateMinDistToNull(fxout2,fyout2);


for r=1:rows
    for c=1:cols
        s=[];
        if strcmp(mode, 'gradient')
            
            s(1)=sqrt(fxin1(r,c)^2 + fyin1(r,c)^2);
            s(2)=sqrt(fxout1(r,c)^2 + fyout1(r,c)^2);
            s(3)=sqrt(fxin2(r,c)^2 + fyin2(r,c)^2);
            s(4)=sqrt(fxout2(r,c)^2 + fyout2(r,c)^2);
            
            [val,ind] = max(s);
            
            %             g(1)=valmap1(r,c);
            %             g(2)=valmap2(r,c);
            %             g(3)=valmap3(r,c);
            %             g(4)=valmap4(r,c);
            %             [val,ind] = min(s);
            
            
        else
           
            s(1)=abs(in1(r,c));
            s(2)=abs(out1(r,c));
            s(3)=abs(in2(r,c));
            s(4)=abs(out2(r,c));
            
            [val,ind] = min(s);
            
            
        end
        
        if ind == 1%3 and 2 or 3 and 4
            if ((abs(in2(r,c)) + abs(out1(r,c))) > (abs(in2(r,c)) + abs(out2(r,c))))
                %if s(2) < s(4)
                
                inPhase(r,c) = inPhase(r,c) * exp(1i * angle(in2(r,c)));
                outPhase1(r,c) = outPhase1(r,c) * exp(1i * angle(out1(r,c)));
                fixPhase(r,c) =  angle(conj(fix1(r,c)).*out1(r,c)) / 2;
            else
                inPhase(r,c) = inPhase(r,c) * exp(1i * angle(in2(r,c)));
                outPhase1(r,c) = outPhase1(r,c) * exp(1i * angle(out2(r,c)));
                fixPhase(r,c) = angle(conj(fix2(r,c)).*out2(r,c)) / 2;
            end
            %             else
            %                 if g(2) > g(4)
            %                   inPhase(r,c) = inPhase(r,c) * exp(1i * angle(in1(r,c)));
            %                   outPhase1(r,c) = outPhase1(r,c) * exp(1i * angle(out1(r,c)));
            %                   fixPhase(r,c) =  angle(conj(fix1(r,c)).*out1(r,c)) / 2;
            %                 else
            %                     inPhase(r,c) = inPhase(r,c) * exp(1i * angle(in1(r,c)));
            %                     outPhase1(r,c) = outPhase1(r,c) * exp(1i * angle(out2(r,c)));
            %                     fixPhase(r,c) =  angle(conj(fix1(r,c)).*out2(r,c)) / 2;
            %                 end
            % end
        elseif ind == 2%1 and 4 or 3 and 4
            %if ((s(1) + s(4)) < (s(3) + s(4)))
            %if ((abs(in1(r,c)) + abs(out2(r,c))) > (abs(in2(r,c)) + abs(out2(r,c))))
            %if s(1) < s(3)
            %   inPhase(r,c) = inPhase(r,c) * exp(1i * angle(in1(r,c)));
            %  outPhase1(r,c) = outPhase1(r,c) * exp(1i * angle(out2(r,c)));
            % fixPhase(r,c) = angle(conj(fix2(r,c)).*out2(r,c)) / 2;
            % else
            %if gind ~= 2
            inPhase(r,c) = inPhase(r,c) * exp(1i * angle(in2(r,c)));
            outPhase1(r,c) = outPhase1(r,c) * exp(1i * angle(out2(r,c)));
            fixPhase(r,c) = angle(conj(fix2(r,c)).*out2(r,c)) / 2;
            %             else
            %                 inPhase(r,c) = inPhase(r,c) * exp(1i * angle(in2(r,c)));
            %                 outPhase1(r,c) = outPhase1(r,c) * exp(1i * angle(out1(r,c)));
            %                 fixPhase(r,c) = angle(conj(fix2(r,c)).*out1(r,c)) / 2;
            %             end
            %end
        elseif ind ==3%1 and 2 or 1 and 4
            %if ((s(1) + s(2)) < (s(1) + s(4)))
            %if ((abs(in1(r,c)) + abs(out1(r,c))) > (abs(in1(r,c)) + abs(out2(r,c))))
            %if s(2) < s(4)
            %if gind ~= 3
            inPhase(r,c) = inPhase(r,c) * exp(1i * angle(in1(r,c)));
            outPhase1(r,c) = outPhase1(r,c) * exp(1i * angle(out1(r,c)));
            fixPhase(r,c) = angle(conj(fix1(r,c)).*out1(r,c)) / 2;
            %             else
            %                 inPhase(r,c) = inPhase(r,c) * exp(1i * angle(in2(r,c)));
            %                 outPhase1(r,c) = outPhase1(r,c) * exp(1i * angle(out1(r,c)));
            %                 fixPhase(r,c) = angle(conj(fix1(r,c)).*out1(r,c)) / 2;
            %             end
            %else
            %   inPhase(r,c) = inPhase(r,c) * exp(1i * angle(in1(r,c)));
            %  outPhase1(r,c) = outPhase1(r,c) * exp(1i * angle(out2(r,c)));
            % fixPhase(r,c) = angle(conj(fix2(r,c)).*out2(r,c)) / 2;
            % end
        elseif ind ==4%1 and 2 or 3 and 2
            if ((abs(in1(r,c)) + abs(out1(r,c))) > (abs(in2(r,c)) + abs(out1(r,c))))
                %if s(1) < s(3)
                %if gind ~= 4
                % if valmap1(r,c) > valmap3(r,c)
                inPhase(r,c) = inPhase(r,c) * exp(1i * angle(in1(r,c)));
                outPhase1(r,c) = outPhase1(r,c) * exp(1i * angle(out1(r,c)));
                fixPhase(r,c) = angle(conj(fix1(r,c)).*out1(r,c)) / 2;
            else
                inPhase(r,c) = inPhase(r,c) * exp(1i * angle(in2(r,c)));
                outPhase1(r,c) = outPhase1(r,c) * exp(1i * angle(out1(r,c)));
                fixPhase(r,c) = angle(conj(fix1(r,c)).*out1(r,c)) / 2;
            end
            %             else
            %                 if valmap1(r,c) > valmap3(r,c)
            %                     inPhase(r,c) = inPhase(r,c) * exp(1i * angle(in1(r,c)));
            %                     outPhase1(r,c) = outPhase1(r,c) * exp(1i * angle(out2(r,c)));
            %                     fixPhase(r,c) = angle(conj(fix1(r,c)).*out2(r,c)) / 2;
            %                 else
            %                     inPhase(r,c) = inPhase(r,c) * exp(1i * angle(in2(r,c)));
            %                     outPhase1(r,c) = outPhase1(r,c) * exp(1i * angle(out2(r,c)));
            %                     fixPhase(r,c) = angle(conj(fix1(r,c)).*out2(r,c)) / 2;
            %                 end
            %             end
        end
    end
end
for r=1:rows
    for c=1:cols
        %final reconstruct
        [water(r,c), fat(r,c)] = DixonSeparationWithFieldMapPhase( inPhase(r,c) , outPhase1(r,c), fixPhase(r,c) );
    end
end


%{
   % correct mag for all images
for r=1:rows
    for c=1:cols
        s(1)=abs(in1(r,c));
        s(2)=abs(out1(r,c));
        s(3)=abs(in2(r,c));
        s(4)=abs(out2(r,c));
        
        [val,ind] = min(s);
        
        if ind ==1
            inPhase(r,c) = (s(2)+s(3)+s(4))/3;
            outPhase1(r,c) = (s(2)+s(3)+s(4))/3;
            %outPhase2(r,c) = (s(2)+s(3)+s(4))/3;
        elseif ind ==2
            inPhase(r,c) = (s(1)+s(3)+s(4))/3;
            outPhase1(r,c) = (s(1)+s(3)+s(4))/3;
            %outPhase2(r,c) = (s(1)+s(3)+s(4))/3;
        elseif ind ==3
            inPhase(r,c) = (s(2)+s(1)+s(4))/3;
            outPhase1(r,c) = (s(2)+s(1)+s(4))/3;
            %outPhase2(r,c) = (s(2)+s(1)+s(4))/3;
        else
            inPhase(r,c) = (s(2)+s(3)+s(1))/3;
            outPhase1(r,c) = (s(2)+s(3)+s(1))/3;
            %outPhase2(r,c) = (s(2)+s(3)+s(1))/3;
        end
    end
end

%find phase for inphase image
[fxin1,fyin1]=gradient(angle(in1));
[fxin2,fyin2]=gradient(angle(in2));
[fxout1,fyout1]=gradient(angle(in1));
[fxout2,fyout2]=gradient(angle(in2));

for r=1:rows
    for c=1:cols
        s=[];
        s(1)=sqrt(fxin1(r,c)^2 + fyin1(r,c)^2);
        s(2)=sqrt(fxout1(r,c)^2 + fyout1(r,c)^2);
        s(3)=sqrt(fxin2(r,c)^2 + fyin2(r,c)^2);
        s(4)=sqrt(fxout2(r,c)^2 + fyout2(r,c)^2);
        
        [val,ind] = max(s);
        
        if ind == 1%3 and 2 or 3 and 4
            if ((s(3) + s(2)) < (s(3) + s(4)))
                inPhase(r,c) = inPhase(r,c) * exp(1i * angle(in2(r,c)));
                outPhase1(r,c) = outPhase1(r,c) * exp(1i * angle(out1(r,c)));
                fixPhase(r,c) =  angle(conj(fix1(r,c)).*out1(r,c)) / 2;
            else
                inPhase(r,c) = inPhase(r,c) * exp(1i * angle(in2(r,c)));
                outPhase1(r,c) = outPhase1(r,c) * exp(1i * angle(out2(r,c)));
                fixPhase(r,c) = angle(conj(fix2(r,c)).*out2(r,c)) / 2;
            end
        elseif ind == 2%1 and 4 or 3 and 4
            if ((s(1) + s(4)) < (s(3) + s(4)))
                inPhase(r,c) = inPhase(r,c) * exp(1i * angle(in1(r,c)));
                outPhase1(r,c) = outPhase1(r,c) * exp(1i * angle(out2(r,c)));
                fixPhase(r,c) = angle(conj(fix2(r,c)).*out2(r,c)) / 2;
            else
                inPhase(r,c) = inPhase(r,c) * exp(1i * angle(in2(r,c)));
                outPhase1(r,c) = outPhase1(r,c) * exp(1i * angle(out2(r,c)));
                fixPhase(r,c) = angle(conj(fix2(r,c)).*out2(r,c)) / 2;
            end
        elseif ind ==3%1 and 2 or 1 and 4
            if ((s(1) + s(2)) < (s(1) + s(4)))
                inPhase(r,c) = inPhase(r,c) * exp(1i * angle(in1(r,c)));
                outPhase1(r,c) = outPhase1(r,c) * exp(1i * angle(out1(r,c)));
                fixPhase(r,c) = angle(conj(fix1(r,c)).*out1(r,c)) / 2;
            else
                inPhase(r,c) = inPhase(r,c) * exp(1i * angle(in1(r,c)));
                outPhase1(r,c) = outPhase1(r,c) * exp(1i * angle(out2(r,c)));
                fixPhase(r,c) = angle(conj(fix2(r,c)).*out2(r,c)) / 2;
            end
        elseif ind ==4%1 and 2 or 3 and 2
            if ((s(1) + s(2)) < (s(3) + s(2)))
                inPhase(r,c) = inPhase(r,c) * exp(1i * angle(in1(r,c)));
                outPhase1(r,c) = outPhase1(r,c) * exp(1i * angle(out1(r,c)));
                fixPhase(r,c) = angle(conj(fix1(r,c)).*out1(r,c)) / 2;
            else
                inPhase(r,c) = inPhase(r,c) * exp(1i * angle(in2(r,c)));
                outPhase1(r,c) = outPhase1(r,c) * exp(1i * angle(out1(r,c)));
                fixPhase(r,c) = angle(conj(fix1(r,c)).*out1(r,c)) / 2;
            end
        end
    end
end
for r=1:rows
    for c=1:cols
%final reconstruct
    [water(r,c), fat(r,c)] = DixonSeparationWithFieldMapPhase( inPhase(r,c) , outPhase1(r,c), fixPhase(r,c) );
    end
end

%}