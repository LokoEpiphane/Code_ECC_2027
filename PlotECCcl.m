
clear; clc; close all;

%% Parametres
T = 100;
%% 
%% 
Nx = 100;
x = linspace(0,1,Nx+1);
dx = (x(2)-x(1));

lambda = [1; sqrt(2)];
H = [0.25, 0; 
       -1, 1.01];
B = eye(2);

CFL = 0.8;
dt = CFL*dx/max(lambda);
Nt = ceil(T/dt);
dt = T/Nt;
t = linspace(0,T,Nt+1);
normU = zeros(1,Nt+1);

%% Condition initiale
Y = zeros(2,Nx+1);
Y(1,:) = cos(4*pi*x)-1;
Y(2,:) = cos(2*pi*x)-1;

normL2 = zeros(1,Nt+1);
normL2(1) = sqrt(trapz(x, sum(Y.^2,1)));


function u = tanh_vector(x)
    u=tanh(x);
end

%% Boucle temporelle
for n = 1:Nt
    tn = t(n);

    % Commande a partir de la trace Y(t,1)
    y = Y(:,end);
    z = -y/36;
    U = 1.5*tanh_vector(z);
    

    % Source
    if tn <= 10
        d = sin(tn*x);
    else
        d = zeros(size(x));
    end

    Ynew = Y;

    % Lax-Friedrichs aux points interieurs
    % j = 2,...,Nx
    for j = 2:Nx
        Ynew(:,j) = 0.5*(Y(:,j+1)+Y(:,j-1)) ...
            - dt/(2*dx) * ...
            (lambda.*(Y(:,j+1)-Y(:,j-1))) ...
            + dt*[1;1]*d(j);
    end

    % Traitement sortant a x=1 :
     %derivee amont (les deux vitesses sont positives)
    Ynew(:,end) = Y(:,end) ...
        - dt/dx * (lambda.*(Y(:,end)-Y(:,end-1))) ...
        + dt*[1;1]*d(end);
     
    %yt=Ynew(:,end);
   % zt=-yt/18;
   % tanhyt = zt.*(exp(2*zt)-1)./(exp(2*zt)+1);
   % Ut = 1.5*tanhyt;
   % normU(n) = norm(Ut, 2);

    % Condition d'entree a x=0
    Ynew(:,1) = H*Y(:,end) + B*U;

    % Mise a jour
    Y = Ynew;

    % Norme L2 : integration numerique en espace
    normL2(n+1) = sqrt(trapz(x, sum(Y.^2,1)));

end




%% Graphe
figure('Color','w');

plot(t,normL2,'b-','LineWidth',2);
%hold on;
%plot(t,normL2_ol,'r--','LineWidth',2);

xlabel('$t$','Interpreter','latex','FontSize',22);
ylabel('$\|Y(t,.)\|_{L^2(0,1)}$','Interpreter','latex','FontSize',22);
%title('Closed-loop');
legend('  Closed-loop system', ...
      'Location','best','Fontsize',20);
set(gca,'FontSize',18)
grid on;
xlim([0,T]);

%% Graphe de la norme de U(t)
%figure('Color','w');

%plot(t(1:Nt), normU(1:Nt), 'r-', 'LineWidth', 2.5);

%xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 22);
%ylabel('$\|U(t)\|_2$', 'Interpreter', 'latex', 'FontSize', 22);

%grid on;
%xlim([0,T]);



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
