%set up time and space parameters
timestep = 0.01;
totaltime = 30;
time = (0:timestep:totaltime);
numtsteps = totaltime/timestep;
distancestep = 0.02;
totaldistance = 3;
distance = (distancestep:distancestep:totaldistance);
numdsteps = totaldistance/distancestep;

%initialize protein and RNA matrices: rows are distance, columns are time
Xrna = zeros(numdsteps, numtsteps);
Xprot = zeros(numdsteps, numtsteps);

%initial conditions
Xrna(1,1) = 0.5;
Xprot(1,1) = 0.5;

for (t=1:numtsteps)
    for (i=1:numdsteps)
        
        
       
        %find laplacian for each species with absorbing boundary conditions
        switch i
            case 1
                lap_rna = (Xrna(2,t) - 2*Xrna(1,t))/ distancestep^2;
                lap_prot = (Xprot(2,t) - 2*Xprot(1,t))/ distancestep^2;
            case numdsteps
                lap_rna = (Xrna(numdsteps-1,t)-2*Xrna(numdsteps,t)) / distancestep^2;
                lap_prot = (Xprot(numdsteps-1,t)-2*Xprot(numdsteps,t)) / distancestep^2;
            otherwise               
                lap_rna = (Xrna(i+1, t) + Xrna(i-1, t) - 2*Xrna(i, t)) / distancestep^2;
                lap_prot = (Xprot(i+1, t) + Xprot(i-1, t) - 2*Xprot(i, t)) / distancestep^2;
        end
        
        
        %find rates for each species
        dXprot_dt = Xrna(i, t) - Xprot(i, t) + 0.0001*lap_prot;
        dXrna_dt = ( (Xprot(i,t)^2)/((0.33^2) + (Xprot(i,t)^2)) ) - Xrna(i,t) + 0.0001*lap_rna;
        
        %forward euler
        Xprot(i, t+1) = Xprot(i,t) + dXprot_dt*timestep;
        Xrna(i, t+1) = Xrna(i,t) + dXrna_dt*timestep;
        
    end    
    
end
hold on
figure(1)
plot(time, Xprot(1,:), 'DisplayName', 'Protein')
plot(time, Xrna(1,:), 'DisplayName', 'RNA')
title('Concentration of protein and RNA vs time for first position in space')
xlabel('Time (s)')
ylabel('Concentration (mM)')
legend()
hold off

figure(2)
plot(time, Xprot(2,:), 'DisplayName', 'Protein')
hold on
plot(time, Xrna(2,:), 'DisplayName', 'RNA')
title('Concentration of protein and RNA vs time for second position in space')
xlabel('Time (s)')
ylabel('Concentration (mM)')
legend()
hold off

figure(3)
plot(time, Xprot(3,:), 'DisplayName', 'Protein')
hold on
plot(time, Xrna(3,:), 'DisplayName', 'RNA')
title('Concentration of protein and RNA vs time for third position in space')
xlabel('Time (s)')
ylabel('Concentration (mM)')
legend()
hold off

figure(4)
plot(distance, Xprot(:,1), 'DisplayName', 'Protein')
hold on
plot(distance, Xrna(:,1), 'DisplayName', 'RNA')
title('Concentration of protein and RNA vs distance for first point in time')
xlabel('Distance (um)')
ylabel('Concentration (mM)')
legend()
hold off

figure(5)
plot(distance, Xprot(:,numtsteps), 'DisplayName', 'Protein')
hold on
plot(distance, Xrna(:,numtsteps), 'DisplayName', 'RNA')
title('Concentration of protein and RNA vs distance for last point in time')
xlabel('Distance (um)')
ylabel('Concentration (mM)')
hold off
legend()
