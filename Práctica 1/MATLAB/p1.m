clc;
clear;

A = randomMatrix(6, 6);
detA = det(A);
disp("Matriz A de 6x6: ");
disp(A);
disp("Determinante de A: " + detA);
disp(" ");
disp("Factorización LU de A: ")
luA = lu(A);
disp(luA);

function X = randomMatrix(m, n)
    X = randi([1,10], m,n);
end

%{
9     3     5     3     7     2
6    10     2     6     2     1
10    5     9     5     9     8
5     3     6     8     2     8
4     5     4     1     7     5
5     8     3     9     6     1
%}
