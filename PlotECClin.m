
clear; clc; close all;

%% Parametres
T = 100;
Nx = 100;
x = linspace(0,1,Nx+1);
dx = (x(2)-x(1));

lambda = [1; sqrt(2)];
H = [0.25, 0; 
       -1, 1.01];
B = eye(2);
K = -(1/4)*eye(2);

CFL = 0.8;
dt = CFL*dx/max(lambda);
Nt = ceil(T/dt);
dt = T/Nt;
t = linspace(0,T,Nt+1);

%% Condition initiale
Y = zeros(2,Nx+1);
Y(1,:) = cos(4*pi*x)-1;
Y(2,:) = cos(2*pi*x)-1;

normL2(1) = sqrt(trapz(x, sum(Y.^2,1)));

%% Simulation en boucle ouverte : U(t) = 0

Y_ol = zeros(2,Nx+1);
Y_ol(1,:) = cos(4*pi*x)-1;
Y_ol(2,:) = cos(2*pi*x)-1;

normL2_ol = zeros(1,Nt+1);
normL2_ol(1) = sqrt(trapz(x, sum(Y_ol.^2,1)));

for k = 1:Nt

    tk = t(k);

    % Terme source
    if tk <= 10
        d = sin(tk*x);
    else
        d = zeros(size(x));
    end

    Ynewol = Y_ol;
    
    % Schema de Lax-Friedrichs aux points interieurs
    for j = 2:Nx
        Ynewol(:,j) = 0.5*(Y_ol(:,j+1)+Y_ol(:,j-1)) ...
            - dt/(2*dx)*(lambda.* ...
            (Y_ol(:,j+1)-Y_ol(:,j-1))) ...
            + dt*[1;1]*d(j);
    end

    % Bord sortant x=1
    Ynewol(:,end) = Y_ol(:,end) ...
        - dt/dx*(lambda.* ...
        (Y_ol(:,end)-Y_ol(:,end-1))) ...
        + dt*[1;1]*d(end);

    % Bord entrant x=0, avec U(t)=0
    Ynewol(:,1) = H*Y_ol(:,end);

    Y_ol = Ynewol;

    % Norme L2 en boucle ouverte
    normL2_ol(k+1) = sqrt(trapz(x, sum(Y_ol.^2,1)));

end

%% Graphe
figure('Color','w');

plot(t,normL2_ol,'b-','LineWidth',2);

xlabel('$t$','Interpreter','Latex','FontSize',20);
ylabel('$\|Y(t,.)\|_{L^2(0,1)}$','Interpreter','Latex','FontSize',20);
%title('Co$mparaison boucle fermee / boucle ouverte');
legend('Open-loop system : U(t) = 0', ...
       'Location','best','Fontsize',20);
set(gca,'FontSize',18)

grid on;
xlim([0,T]);



% Graphe de la norme L2
%figure('Color','w');
%plot(t,normL2,'b-','LineWidth',2);
%xlabel('Temps t','FontSize',12);
%ylabel('||Y(t,.)||_{L^2(0,1)}','FontSize',12);
%title('Evolution de la norme L^2 de la solution','FontSize',13);
%grid on;
%xlim([0,T]);

% Affichage de quelques valeurs
%fprintf('Norme L2 initiale : %.7f\n',normL2(1));
%fprintf('Norme L2 finale   : %.7f\n',normL2(end));
%fprintf('Norme L2 maximale : %.7f\n',max(normL2));
