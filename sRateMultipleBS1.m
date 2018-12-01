
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
rate = 100;% 100bps   
W = 100;% 100Hz
N = W*1.38*10^(-23)*298;%puissance du bruit

MonteCarlo = 1000;
MonteCarlo2 = 1000;
NbNoeudsMax = 20;
NbBSMax = 4;

% matrice de simulations
% les colonnes représentent les numéros d'itérations
% les lignes sont par blocks de NbNoeudsMax
% il y a NbBSMax blocks

simulations = zeros(NbNoeudsMax*NbBSMax,MonteCarlo);

% matrice de résultats issus d'un grand nombre de simulations de
% disposition de BS
results = zeros(NbBSMax,NbNoeudsMax,MonteCarlo2);

% matrice de résultats issus d'un grand nombre de simulations de
% dispositions de BS et de Noeuds
finalResults = zeros(NbBSMax,NbNoeudsMax);


% dimensions de l'espace considéré en 
dimX = 200000; 
dimY = 200000;

for ite=1:MonteCarlo2
    for n = 1:NbNoeudsMax

        % création d'une matrice regroupant les positions des noeuds
        % dimensions de cette matrice (2,n)
        % 1ère ligne => postion en x
        % 2ème ligne => postion en y
        Noeuds = zeros(2,n);

        % positionnement des neuds, 

        for i=1:n
            Noeuds(1,i) = rand*dimX;
            Noeuds(2,i) = rand*dimY;
        end

        for b = 1:NbBSMax
            for it = 1:MonteCarlo

                % création d'une matrice avec les informations sur les BS
                % 1ère ligne => postion en x
                % 2ème ligne => postion en y
                BS = zeros(2,b);

                % positionnement des BS

                for j =1:b
                    BS(1,j) = rand*dimX;
                    BS(2,j) = rand*dimY;
                end

                % matrice contenant les distance des noeuds par rapport aux NbBS BS
                % dimensions de cette matrice (NbBS,NbNoeuds)
                distances = zeros(b,n);

                % matrice contenant les puissances reçues par les BS en fonction des noeuds
                % dimensions de cette matrice (NbBS,NbNoeuds)
                powerReceived = zeros(b,n);

                % instanciation des puissances reçues pour chaque noeud par
                % rapport à chaque BS  

                for i=1:n
                    for j = 1:b
                        distances(j,i) = sqrt((Noeuds(1,i)-BS(1,j))^2+(Noeuds(2,i)-BS(2,j))^2);
                        powerReceived(j,i) = PowerReceived(distances(j,i));
                    end
                end

                % matrice ordonant les puissances de la plus élevée à la plus faible
                % conserve aussi les indices des noeuds liés à ces puissances
                % NbBS premières lignes pour les puissances ordonnées 
                % NbBS secondes lignes pour les indices des noeuds

                orderedPowers = zeros(2*b,n);

                % vecteur indiquant si le signal transmis par le
                % noeud a été décodé par une base donnée (0 si décodé)
                % on tient compte de l'indexe du tableau matlab pour les noeuds

                decoded = zeros(b,n);
                for i =1:b
                    decoded(i,:)=1;
                end

                % instanciation des indices des noeuds
                for i =1:b
                    for j=1:n
                        orderedPowers(b+i,j)=j;
                    end
                end

                for i =1:b
                    % création d'une matrice avec les puissances ordonnées
                    orderedPowers(i,:) = powerReceived(i,:);
                    change = 1;
                    while(change ~= 0)
                        change = 0;
                        for j=1:(n-1)
                            if (orderedPowers(i,j+1) > orderedPowers(i,j))
                                tempPower = orderedPowers(i,j);
                                tempIndex = orderedPowers(b+i,j);
                                orderedPowers(i,j)= orderedPowers(i,j+1);
                                orderedPowers(b+i,j) = orderedPowers(b+i,j+1);
                                orderedPowers(i,j+1) = tempPower;
                                orderedPowers(b+i,j+1) = tempIndex;
                                change = change+1;
                            end
                        end
                    end
                end

                % verification de la première condition pour n noeuds classés de la
                % plus grande puissance à la plus faible
                % reçue à la plus faible, pour un noueud i Ri <= log2(1 + Pi/N + SUM(P(i+1) à P(NbNoeuds))  
                for i =1:b
                    SIC =0;
                    for j=1:n
                        interference = 0;
                        for k = j+1:n
                            interference = interference + orderedPowers(i,k);
                        end
                        if ( 6.8 <= (10*(log10(orderedPowers(i,j)/(interference+N)))))
                            if(SIC~=1)
                                decoded(i,(orderedPowers(i+b,j))) = 0;
                            end
                        else
                            SIC =1;
                        end

                    end
                end


                % remplissage du vecteur simulations pour un nombe de noeuds
                % avec un nombre de BS donné, on a un nombre de bonnes
                % transmissions

                finalDecoded = zeros(1,n);
                finalDecoded(1,:) = 1;

                for i=1:b
                    for j=1:n
                        if(decoded(i,j) == 0)
                            finalDecoded(1,j) = 0;
                        end
                    end
                end
                for j=1:n
                    if(finalDecoded(1,j) == 0)
                        simulations(((b-1)*NbNoeudsMax+n),it) = simulations(((b-1)*NbNoeudsMax+n),it)+1;
                    end
                end
            end
        end
    end

    for b = 1:NbBSMax
        for n = 1:NbNoeudsMax
            % calcul du taux de réussite 
            results(b,n,ite) = ((sum(simulations(((b-1)*NbNoeudsMax+n),:)))/(MonteCarlo*n));
        end
    end
    simulations(:,:)=0;
end

for b = 1:NbBSMax
    for n = 1:NbNoeudsMax
        % calcul du taux de réussite 
        finalResults(b,n) = (sum(results(b,n,:))/(MonteCarlo2));
    end
end



figure
hold on;
for i=1:NbBSMax
    semilogy(finalResults(i,:),'-o');
end
legend('NbBS = 1','NbBS = 2','NbBS = 3','NbBS = 4');
%title('rate of successfully decoded transmissions vs number of active nodes for multiple BS');
xlabel('number of active nodes');
ylabel('rate of successfully decoded transmissions');

function Pr = PowerReceived(distance)
    global Pt;
    global Gt;
    global Gr;
    global lambda;
    
    Pr = (Pt * Gt * Gr * (lambda / (4 * pi * distance))^2);
end 
