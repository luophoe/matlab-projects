format shortg;
clc;

% Ask User to Input R value
prompt = 'Input the resistance value: ';
R = input(prompt);

% Bisection Method
Tm_l = 0; % keep track of the last Tm for error calculation
error_1 = 1; % absolute relative approximate error percentage
count_1 = 0; % keep track of the number of iterations

if R > 100
    min = 0; % the left range for the iteration if inputed R is higher than 100
    max = 850; % the right range for the iteration if inputed R is higher than 100
    fmin = 100 * (1 + 3.9083 * 10^-3 * min - 5.775 * 10^-7 * min^2) - R; % f of left range
    fmax = 100 * (1 + 3.9083 * 10^-3 * max - 5.775 * 10^-7 * max^2) - R; % f of left range
    while error_1 > 0.05
        Tm = (min + max)/2;  % calculated temperature Tm for the iteration
        fTm = 100 * (1 + 3.9083 * 10^-3 * Tm - 5.775 * 10^-7 * Tm^2) - R;
        if fmin * fTm < 0 % if fmin * fTm < 0, change range to be [min, Tm]
            max = Tm;
        else % else, change range to be [Tm, max]
            min = Tm;
        end
        error_1 = abs((Tm - Tm_l) / Tm) * 100;
        Tm_l = Tm; % update Tm_l to keep track of Tm
        count_1 = count_1 + 1;
    end
else
    min = -200; % the left range for the iteration if inputed R is lower than 100
    max = 0; % the right range for the iteration if inputed R is lower than 100
    fmin = 100 * (1 + 3.9083 * 10^-3 * min - 5.775 * 10^-7 * min^2 - 4.183 * 10^-12 * (min - 100) * min^3) - R; % f of left range
    fmax = 100 * (1 + 3.9083 * 10^-3 * max - 5.775 * 10^-7 * max^2 - 4.183 * 10^-12 * (max - 100) * max^3) - R; % f of left range
    while error_1 > 0.05
        Tm = (min + max)/2;  % calculated temperature Tm for the iteration
        fTm = 100 * (1 + 3.9083 * 10^-3 * Tm - 5.775 * 10^-7 * Tm^2 - 4.183 * 10^-12 * (Tm - 100) * Tm^3) - R;
        if fmin * fTm < 0 % if fmin * fTm < 0, change range to be [min, Tm]
            max = Tm;
        else % else, change range to be [Tm, max]
            min = Tm;
        end
        error_1 = abs((Tm - Tm_l) / Tm) * 100;
        Tm_l = Tm; % update Tm_l to keep track of Tm
        count_1 = count_1 + 1;
    end
end

% Newton Raphson Method
error_2 = 1; % absolute relative approximate error percentage
count_2 = 1; % keep track of the number of iterations (starts at 1 for the 1st iteration outside of the while loop)

if R > 100
    T_i = 450; % the initial T of range [0, 850)
    T = T_i - (5.775 * 10^-7 * T_i^2 - 3.9083 * 10^-3 * T_i + (R/100) - 1) / (2 * 5.775 * 10^-7 * T_i - 3.9083 * 10^-3); 
    % get T of first iteration using T = T_i - f(T_i)/f'(T_i)
    T_l = T; % update T_l to keep track of T
    while error_2 > 0.05
        T = T_l - (5.775 * 10^-7 * T_l^2 - 3.9083 * 10^-3 * T_l + (R/100) - 1)/(2 * 5.775 * 10^-7 * T_l - 3.9083 * 10^-3);
        % get T using T = T_l - f(T_l)/f'(T_l)
        error_2 = abs((T - T_l) / T) * 100;
        T_l = T; % update T_l to keep track of T
        count_2 = count_2 + 1;
    end
else
    T_i = -100; % the initial T of range (-200, 0)
    T = T_i - (4.183 * 10^-12 * T_i^4 - 4.183 * 10^-10 * T_i^3 + 5.775 * 10^-7 * T_i^2 - 3.9083 * 10^-3 * T_i + R/100 -1) / (4 * 4.183 * 10^-12 * T_i^3 - 3 * 4.183 * 10^-10 * T_i^2 + 2 * 5.775 * 10^-7 * T_i - 3.9083 * 10^-3);
    % get T of first iteration using T = T_i - f(T_i)/f'(T_i)
    T_l = T; % update T_l to keep track of T
    while error_2 > 0.05
        T = T_l - (4.183 * 10^-12 * T_l^4 - 4.183 * 10^-10 * T_l^3 + 5.775 * 10^-7 * T_l^2 - 3.9083 * 10^-3 * T_l + R/100 -1) / (4 * 4.183 * 10^-12 * T_l^3 - 3 * 4.183 * 10^-10 * T_l^2 + 2 * 5.775 * 10^-7 * T_l - 3.9083 * 10^-3);
        % get T using T = T_l - f(T_l)/f'(T_l)
        error_2 = abs((T - T_l) / T) * 100;
        T_l = T; % update T_l to keep track of T
        count_2 = count_2 + 1;
    end
end

DISP1 = ['The temperature obtained by bisection is ', num2str(Tm), ' C.'];
DISP2 = ['The temperature obtained by NR is ', num2str(T), ' C.'];
DISP3 = ['The number of required iterations for bisection is ', num2str(count_1), '.'];
DISP4 = ['The number of required iterations for NR is ', num2str(count_2), '.'];
DISP5 = ['The absolute relative approximate error % for bisection is ', num2str(error_1), '%.'];
DISP6 = ['The absolute relative approximate error % for NR is ', num2str(error_2), '%.'];

disp(DISP1);
disp(DISP2);
disp(DISP3);
disp(DISP4);
disp(DISP5);
disp(DISP6);