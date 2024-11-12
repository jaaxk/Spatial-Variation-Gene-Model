%set up time and space parameters
timestep = 0.01;
times = [50, 100, 200, 400];
fig=1;

for (k=1:4)
totaltime=times(k);
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
Xrna(1,1) = 1;
Xprot(1,1) = 1;

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




figure(fig)
plot(time, Xprot(1,:), 'r', 'DisplayName', 'Position 1 Protein')
hold on
plot(time, Xrna(1,:),'r','LineStyle','--','DisplayName', 'Position 1 RNA')

plot(time, Xprot(2,:), 'b','DisplayName', 'Position 2 Protein')
plot(time, Xrna(2,:), 'b','LineStyle','--','DisplayName', 'Position 2 RNA')

plot(time, Xprot(3,:), 'g','DisplayName', 'Position 3 Protein')
plot(time, Xrna(3,:), 'g','LineStyle','--', 'DisplayName', 'Position 3 RNA')
title(sprintf('Concentration vs time for 1st, 2nd and 3rd positions in space, total time = %ds', totaltime))
xlabel('Time (s)')
ylabel('Concentration (mM)')
legend()
hold off

fig = fig + 1;
figure(fig)
plot(distance, Xprot(:,1), 'r','DisplayName', 'First time point: Protein')
hold on
plot(distance, Xrna(:,1),'r','LineStyle','--', 'DisplayName', 'First time point: RNA')

plot(distance, Xprot(:,numtsteps), 'b','DisplayName', 'Last time point: Protein')
plot(distance, Xrna(:,numtsteps), 'b','LineStyle','--', 'DisplayName', 'Last time point: RNA')
title(sprintf('Concentration vs distance for first and last point in time, total time = %ds' ,totaltime))
xlabel('Distance (um)')
ylabel('Concentration (mM)')
hold off
legend()

fig = fig + 1;

end
