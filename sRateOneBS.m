clear all;
close all;


global freq;
global lambda;
global c;
global Pt;
global Gt;
global Gr;

Gr = 1;
Gt = 1;
Pt = 0.1;% 0.1W (puissance de sigfox)
freq = 868*10^6;
c = 3.0*10^8;
lambda = c/freq;
NbBS = 1;
rate = 100;% 100bps   
W = 100;% 100Hz
N = W*1.38*10^(-23)*298;%puissance du bruit

MonteCarlo = 10000;
NbNoeudsMax = 20;

% matrice de simulations
% les colonnes représentent les numéros d'itérations
% les lignes représentent le nombre de noeuds

simulations = zeros(NbNoeudsMax,1);

% matrice de résultats 
results = zeros(NbNoeudsMax,1);

% création d'une matrice avec les informations sur les BS
% 1ère ligne => postion en x
% 2ème ligne => postion en y
BS = zeros(2,NbBS);

% dimensions de l'espace considéré en 
dimX = 200000; %10 km
dimY = 200000;

% positionnement de la BS

for j =1:NbBS
    BS(1,j) = dimX/2;
    BS(2,j) = dimY/2;
end


for n = 1:NbNoeudsMax
    success = 0;
    for it = 1:MonteCarlo

        % positionnement des neuds, on considère que tous sont actifs
        
        % création d'une matrice regroupant les positions des noeuds
        % dimensions de cette matrice (2,NbNoeuds)
        % 1ère ligne => postion en x
        % 2ème ligne => postion en y
        Noeuds = zeros(2,n);

        % matrice contenant les distance des noeuds par rapport aux NbBS BS
        % dimensions de cette matrice (NbBS,NbNoeuds)
        distances = zeros(NbBS,n);

        % matrice contenant les puissances reçues par les BS en fonction des noeuds
        % dimensions de cette matrice (NbBS,NbNoeuds)
        powerReceived = zeros(NbBS,n);

        for i=1:n
            Noeuds(1,i) = rand*dimX;
            Noeuds(2,i) = rand*dimY;
            for j = 1:NbBS
                distances(j,i) = sqrt((Noeuds(1,i)-BS(1,j))^2+(Noeuds(2,i)-BS(2,j))^2);
                powerReceived(j,i) = PowerReceived(distances(j,i));
            end
        end
        % figure
        % hold on
        % plot(BS(1,:),BS(2,:),'+');
        % plot(Noeuds(1,:),Noeuds(2,:),'+');

        % On regarde si chaque station de base arrive à performer du SIC

        % matrice ordonant les puissances de plus élevée à la plus faible
        % conserve aussi les indices des noeuds liés à ces puissances
        % NbBS premières lignes pour les puissances ordonnées 
        % NbBS secondes lignes pour les indices des noeuds
        orderedPowers = zeros(2*NbBS,n);

        % instanciation des indices des noeuds
        for i =1:NbBS
            for j=1:n
                orderedPowers(NbBS+i,j)=j;
            end
        end

        SIC =0;
        for i =1:NbBS
            % création d'une matrice avec les puissances ordonnées
            orderedPowers(i,:) = powerReceived(i,:);
            change = 1;
            while(change ~= 0)
                change = 0;
                for j=1:(n-1)
                    if (orderedPowers(i,j+1) > orderedPowers(i,j))
                        tempPower = orderedPowers(i,j);
                        tempIndex = orderedPowers(NbBS+i,j);
                        orderedPowers(i,j)= orderedPowers(i,j+1);
                        orderedPowers(NbBS+i,j) = orderedPowers(NbBS + i,j+1);
                        orderedPowers(i,j+1) = tempPower;
                        orderedPowers(NbBS + i,j+1) = tempIndex;
                        change = change+1;
                    end
                end
            end

            % verification de la première condition pour n noeuds classés de la
            % plus grande puissance à la plus faible
            % reçue à la plus faible, pour un noueud i Ri <= log2(1 + Pi/N + SUM(P(i+1) à P(NbNoeuds))  

            for j=1:n
                interference = 0;
                for k = j+1:n
                    interference = interference + orderedPowers(i,k);
                end
                if ( 6.8 > (10*(log10(orderedPowers(i,j)/(interference+N)))))
                    SIC = 1;
                else 
                    if(SIC ~= 1)
                        success = success +1;
                    end
                end
            end
        end
    end
    simulations(n,1) = success;
end

for n = 1:NbNoeudsMax
    % calcul du pourcentage de réussite
    results(n) = simulations(n,1)/(MonteCarlo*n);
end

semilogy(results,'-o');
% title('rate of successfully decoded transmissions vs number of active nodes for 1 BS');
legend('NbBS = 1');
xlabel('number of active nodes');
ylabel('rate of successfully decoded transmissions');

function Pr = PowerReceived(distance)
    global Pt;
    global Gt;
    global Gr;
    global lambda;
    
    Pr = (Pt * Gt * Gr * (lambda / (4 * pi * distance))^2);
end 
