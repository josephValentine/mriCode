% not complete
% modifingy an old ssfp to make a gre
% will need to create a GREImage file as well once this is complete.
function [ Mc ] = GRE(beta, M0, alpha, phi, dphi, TR, TE, T1, T2, Nr)
Tr=5*10^(-3);
Te=2.5*10^(-3);
T1=300*10^(-3);
T2=85*10^(-3);
ddphi = 180;

alphastop = 20; %where do we want so sample for first graph
ssfpavg=[];
avggre=[];
savemagxy=[];
savephase=[];
savemagxygre=[];
savephasegre=[];
phi=0;


   

    %--------------------------------------------------------------------------
    %                     GRE SEQUENCE (as a modified ssfp sequence)
    %--------------------------------------------------------------------------

     rxalpha=[1 0 0;0 cos(alpha) sin(alpha);0 -sin(alpha) cos(alpha)];

     Mgre = zeros(3,200);

          %first flip
          phi=0;
          rzdphi=[cos(phi) sin(phi) 0;-sin(phi) cos(phi) 0;0 0 1];
          rznegdphi=[cos(-phi) sin(-phi) 0;-sin(-phi) cos(-phi) 0;0 0 1];
          Mgre(:,1) = [0 0 1];

          Mgre(:,1) = rzdphi*M(:,1);
          Mgre(:,1) = rxalpha*M(:,1);
          Mgre(:,1) = rznegdphi*M(:,1);


            for n=2:200 %assume steady state after 200 flips

                %relaxation
                Mgre(1,n) = Mgre(1,n-1)*exp(-Tr/T2);%x
                Mgre(2,n) = Mgre(2,n-1)*exp(-Tr/T2);%y
                Mgre(3,n) = 1+(Mgre(3,n-1)-1)*exp(-Tr/T1);%z

                %convert offres into angle to tip
                cycles=offres*Tr;
                rotation_angle=rem(cycles,1)*2*pi;
                rzoffres =[cos(rotation_angle) sin(rotation_angle) 0;-sin(rotation_angle) cos(rotation_angle) 0;0 0 1];
                 Mgre(:,n) = rzoffres*Mgre(:,n);%dephase by offres

                 %next tip
                 %delete phase information! to make it gre
                Mgre(1,n)=0;
                Mgre(2,n)=0;


                phi=phi+dphi;%will got over 2pi but shouldnt matter
                
                rzdphi=[cos(phi) sin(phi) 0;-sin(phi) cos(phi) 0;0 0 1];
                rznegdphi=[cos(-phi) sin(-phi) 0;-sin(-phi) cos(-phi) 0;0 0 1];

                Mgre(:,n)=rzdphi*Mgre(:,n);
                Mgre(:,n)=rxalpha*Mgre(:,n);
                Mgre(:,n)=rznegdphi*Mgre(:,n);



            end

            %take steady state sample and relax in Te
            Mgre(1,200)=Mgre(1,200)*exp(-Te/T2);
            Mgre(2,200)=Mgre(2,200)*exp(-Te/T2);


            cycles=offres*Te;
            rotation_angle=rem(cycles,1)*2*pi;
            rzoffres =[cos(rotation_angle) sin(rotation_angle) 0;-sin(rotation_angle) cos(rotation_angle) 0;0 0 1];
            Mgre(:,200) = rzoffres*Mgre(:,200);


            magxygre= [magxygre,sqrt(Mgre(1,200).^2 + Mgre(2,200).^2)];
            phasegre = [phasegre,atan2(Mgre(2,200),Mgre(1,200))];

            if (alpha*(180/pi)==alphastop)
                savemagxygre= [magxygre,sqrt(Mgre(1,200).^2 + Mgre(2,200).^2)];
                savephasegre = [phasegre,atan2(Mgre(2,200),Mgre(1,200))];
            end
