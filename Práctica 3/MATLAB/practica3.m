clear; clc;

global matSym m n;

m = 6; n = 6;

matSym = [43, 59,  51,  50,  39, 54;
          59, 118, 103, 101, 61, 97;
          51, 103, 121, 93,  56, 102;
          50, 101, 93,  111, 69, 99;
          39, 61,  56,  69,  77, 86;
          54, 97,  102, 99,  86, 134]; 

disp("Matriz a operar:");
disp(matSym);
disp(" ");

disp("Ejercicio 3_A (CHOLESKY) --> L es traspuesta de U:");
fprintf("\tLower (Triangulo inferior):\n");
cholesky('lower');
fprintf("\tUpper (Triangulo superior):\n");
cholesky('upper');

disp("Ejercicio 3_B (QR):");
QR();

disp("Ejercicio 3_C (SVD):");
SVD();

disp("Ejercicio 3_D (Autovalores y Autovectores):");
autos();

function cholesky(uplo)
    global matSym;    
    res = chol(matSym, uplo); 
    disp(res);
end

function QR()
    global matSym;
    res = qr(matSym);
    disp(res);
end

function SVD()
    global matSym;
    [U, S, VT] = svd(matSym);
    disp("Valores U");
    disp(U);
    disp("Valores S");
    disp(S);
    disp("Valores VT");
    disp(VT);
end

function autos()
    global matSym;
    [V, D] = eig(matSym);
    disp("Autovalores");
    disp(D);
    disp("Autovectores");
    disp(V);
end

%{
Función para generar matrices cuadradas 
de enteros, simétricas y definidas positivas.
function A = randomMatrix(n)
    A = randi(n,n);
    A = A*A';
    A = A + n*eye(n);
end
%}








