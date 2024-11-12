%set up time and space parameters
timestep = 0.01;
totaltime=5;
time = (0:timestep:totaltime);
numtsteps = totaltime/timestep;
distancestep = 0.02;
totaldistance = 3;
distance = (distancestep:distancestep:totaldistance);
numdsteps = totaldistance/distancestep;

%initialize protein and RNA matrices: rows are distance, columns are time
Xrna = zeros(numdsteps, numtsteps);
Xprot = zeros(numdsteps, numtsteps);
Yrna = zeros(numdsteps, numtsteps);
Yprot = zeros(numdsteps, numtsteps);

%initial conditions
Xrna(1,1) = 1;
Xprot(1,1) = 1;
Yrna(numdsteps,1) = 1;
Yprot(numdsteps, 1) = 1;



for (t=1:numtsteps)
    for (i=1:numdsteps)
        
        
       
        %find laplacian for each species with absorbing boundary conditions
        switch i
            case 1
                lap_Xrna = (Xrna(2,t) - 2*Xrna(1,t))/ distancestep^2;
                lap_Xprot = (Xprot(2,t) - 2*Xprot(1,t))/ distancestep^2;
                lap_Yrna = (Yrna(2,t) - 2*Yrna(1,t))/ distancestep^2;
                lap_Yprot = (Yprot(2,t) - 2*Yprot(1,t))/ distancestep^2;
            case numdsteps
                lap_Xrna = (Xrna(numdsteps-1,t)-2*Xrna(numdsteps,t)) / distancestep^2;
                lap_Xprot = (Xprot(numdsteps-1,t)-2*Xprot(numdsteps,t)) / distancestep^2;
                lap_Yrna = (Yrna(numdsteps-1,t)-2*Yrna(numdsteps,t)) / distancestep^2;
                lap_Yprot = (Yprot(numdsteps-1,t)-2*Yprot(numdsteps,t)) / distancestep^2;
            otherwise               
                lap_Xrna = (Xrna(i+1, t) + Xrna(i-1, t) - 2*Xrna(i, t)) / distancestep^2;
                lap_Xprot = (Xprot(i+1, t) + Xprot(i-1, t) - 2*Xprot(i, t)) / distancestep^2;
                lap_Yrna = (Yrna(i+1, t) + Yrna(i-1, t) - 2*Yrna(i, t)) / distancestep^2;
                lap_Yprot = (Yprot(i+1, t) + Yprot(i-1, t) - 2*Yprot(i, t)) / distancestep^2;
        end
        
        
        %find rates for each species
        dXprot_dt = Xrna(i, t) - Xprot(i, t) + 0.0001*lap_Xprot;
        dXrna_dt = (1 - ( (Yprot(i,t)^2)/((0.33^2) + (Yprot(i,t)^2)) )) - Xrna(i,t) + 0.0001*lap_Xrna;
        dYprot_dt = Yrna(i, t) - Yprot(i, t) + 0.0001*lap_Yprot;
        dYrna_dt = (1 - ( (Xprot(i,t)^2)/((0.33^2) + (Xprot(i,t)^2)) )) - Yrna(i,t) + 0.0001*lap_Yrna;
               
        %forward euler
        Xprot(i, t+1) = Xprot(i,t) + dXprot_dt*timestep;
        Xrna(i, t+1) = Xrna(i,t) + dXrna_dt*timestep;
        Yprot(i, t+1) = Yprot(i,t) + dYprot_dt*timestep;
        Yrna(i, t+1) = Yrna(i,t) + dYrna_dt*timestep;
        
    end    
    
end


figure(1)
plot(time, Xprot(1,:), 'red','DisplayName', 'X Protein')
hold on
plot(time, Xrna(1,:),'red','LineStyle','--', 'DisplayName', 'X RNA')
plot(time, Yprot(1,:),'b', 'DisplayName','Y Protein')
plot(time, Yrna(1,:),'b','LineStyle','--','DisplayName','Y RNA')
title('Concentration of X and Y protein and RNA vs time for first position in space')
xlabel('Time (s)')
ylabel('Concentration (mM)')
legend()
hold off

figure(2)
plot(time, Xprot((numdsteps/2),:), 'red','DisplayName', 'X Protein')
hold on
plot(time, Xrna((numdsteps/2),:),'red','LineStyle','--', 'DisplayName', 'X RNA')
plot(time, Yprot((numdsteps/2),:),'b', 'DisplayName','Y Protein')
plot(time, Yrna((numdsteps/2),:),'b','LineStyle','--','DisplayName','Y RNA')
title('Concentration of X and Y protein and RNA vs time for middle position in space')
xlabel('Time (s)')
ylabel('Concentration (mM)')
legend()
hold off

figure(3)
plot(time, Xprot(numdsteps,:), 'red','DisplayName', 'X Protein')
hold on
plot(time, Xrna(numdsteps,:),'red','LineStyle','--', 'DisplayName', 'X RNA')
plot(time, Yprot(numdsteps,:),'b', 'DisplayName','Y Protein')
plot(time, Yrna(numdsteps,:),'b','LineStyle','--','DisplayName','Y RNA')
title('Concentration of X and Y protein and RNA vs time for last position in space')
xlabel('Time (s)')
ylabel('Concentration (mM)')
legend()
hold off

figure(4)
plot(distance,Xprot(:,1),'red','DisplayName', 'X Protein')
hold on
plot(distance, Xrna(:,1),'red','LineStyle','--', 'DisplayName', 'X RNA')
plot(distance, Yprot(:,1),'b', 'DisplayName','Y Protein')
plot(distance, Yrna(:,1),'b','LineStyle','--','DisplayName','Y RNA')
title('Concentration of X and Y protein and RNA vs distance for first position in time')
xlabel('Distance (um)')
ylabel('Concentration (mM)')
legend()
hold off

figure(5)
plot(distance,Xprot(:,(numtsteps/2)),'red','DisplayName', 'X Protein')
hold on
plot(distance, Xrna(:,(numtsteps/2)),'red','LineStyle','--', 'DisplayName', 'X RNA')
plot(distance, Yprot(:,(numtsteps/2)),'b', 'DisplayName','Y Protein')
plot(distance, Yrna(:,(numtsteps/2)),'b','LineStyle','--','DisplayName','Y RNA')
title('Concentration of X and Y protein and RNA vs distance for middle position in time')
xlabel('Distance (um)')
ylabel('Concentration (mM)')
legend()
hold off

figure(6)
plot(distance,Xprot(:,numtsteps),'red','DisplayName', 'X Protein')
hold on
plot(distance, Xrna(:,numtsteps),'red','LineStyle','--', 'DisplayName', 'X RNA')
plot(distance, Yprot(:,numtsteps),'b', 'DisplayName','Y Protein')
plot(distance, Yrna(:,numtsteps),'b','LineStyle','--','DisplayName','Y RNA')
title('Concentration of X and Y protein and RNA vs distance for last position in time')
xlabel('Distance (um)')
ylabel('Concentration (mM)')
legend()
hold off


