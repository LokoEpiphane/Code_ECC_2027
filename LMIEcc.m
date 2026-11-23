format long
clear; clc; close all;

%% Parametres
theta0 = 0;
theta1 = 0;
alpha0  = -1/2;
alpha1 = -1/6;
mu     = 1;

I2 = eye(2);
Lambda = diag([1, sqrt(2)]);

H = [ 0.25, 0;
     -1,    1.01];

B = eye(2);

K = -(1/2)*I2;
calK = (1/6)*K;

A = H + B*calK;

%% Parametres supplementaires
eta = 1/2;
chi = 1;
Theta = [1; 1];

%% Variables de decision
% Necessaire pour assurer la symetrie de Q*Lambda
q = sdpvar(2,1);
Q = diag(q);
P = sdpvar(2,2,'symmetric');

% Gamma0 diagonale
g = sdpvar(2,1);
Gamma0 = diag(g);
Gamma1 = diag(g);

%% Construction de la matrice bloc
M11 = A'*Q*Lambda*A - exp(-mu)*Q*Lambda;

M12 = A'*Q*Lambda*B + theta0*alpha1*Gamma1*K;

M13 =alpha0* A'*Q*Lambda*B + theta1*Gamma1*K;

M23 =alpha0*B'*Q*Lambda*B+theta0*Gamma0;

M22 = B'*Q*Lambda*B - 2*Gamma0;

M33 = alpha0*alpha0*B'*Q*Lambda*B - 2*Gamma1;


M = [M11,   M12,  M13;
     M12',  M22,  M23;
     M13',  M23', M33];

%% Contraintes
eps0 = 1e-7;

Constraints = [
    Q >= eps0*I2, ...
    Gamma0 >= eps0*I2, ...
    M <= 0, ...
    P >= eps0*I2, ...
    Q*(eta*I2 - mu*Lambda) + P <= 0, ...
    [P, Q*Theta; Theta'*Q, chi^2] >= 0
];

%% Choix du solveur
% Installer YALMIP et un solveur SDP (SDPT3 ou SeDuMi)
options = sdpsettings('verbose', 1, 'solver', 'sdpt3');

%% Resolution
sol = optimize(Constraints, [], options);

%% Resultats
if sol.problem == 0
    Qsol = value(Q);
    Gsol = value(Gamma0);
    Gsol0 = value(Gamma1);
    Msol = value(M);
     Psol = value(P);

    disp('Q ='); disp(Qsol);
    disp('Gamma0 ='); disp(Gsol);
    disp('Gamma1 ='); disp(Gsol0);
    disp('Valeurs propres de Q :');
    disp(eig(Qsol));
    disp('Valeurs propres de Gamma0 :');
    disp(eig(Gsol));
    disp('Valeurs propres de M :');
    disp(eig((Msol+Msol')/2));
    fprintf('\nMatrice P :\n');
    fprintf('%.7f  %.7f\n', Psol');

    fprintf('\nValeurs propres de P :\n');
    disp(eig(Psol));

else
    disp('Aucune solution certifiee par le solveur.');
    disp(yalmiperror(sol.problem));
end
