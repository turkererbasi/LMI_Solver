clear   % Clear Workspace
clc     % Clear comand window

setlmis([])

A0 = [-2.1429, 1142.9, 1.4286, 192.80; 0, 0, 1, 0; 1.8095, -3809.5, -5.0952, -192.80; 0, 0, 0, 0];
A1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, -1];

p_min = 285.7143;
p_max = 400;



for p = p_min : 1 : p_max
    disp(p)
end

n = size(A1, 1);        % Dimension of the matrix A1

K0 = lmivar(1, [n 1]);  % Full symmetric matrix (4x4)
K1 = lmivar(1, [n 1]);  % Full symmetric matrix (4x4)

% Second derivative of the LMI: 2(A1*K1 + K1*A1) < 0
lmiterm([1 1 1 K1], 2*A1, 1, 's');  % LMI #1: 2*A1'*K1
lmiterm([1 1 1 K1], 2, A1, 's');    % LMI #1: 2*K1*A1

lmiterm([-2 1 1 K1], 1, 1, 's');    % LMI #2: K1 
lmiterm([2 1 1 0], 1);              % LMI #2: I 

lmis = getlmis;

options = [0,0,10,0,0];
%[tmin,xfeas] = feasp(lmis,options,-1); 
[tmin,xfeas] = feasp(lmis);

K1_sol = dec2mat(lmis,xfeas,K1);

K1_eig = eig(K1_sol);

if all(K1_eig > 0)
    disp('K1 is positive definite.');
end

Q = A0' * K1 + K1 * A0 + 2 * p_max * (A1 * K1 + K1 * A1);

K0_sol = lyap(A1,Q);