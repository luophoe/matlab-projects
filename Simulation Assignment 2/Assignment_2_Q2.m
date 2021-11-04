clc;

% Symbolic Expressions
syms x y z
f1 = x^3 - 10*x + y - z + 3;
f2 = y^3 + 10*y - 2*x - 2*z - 5;
f3 = x + y - 10*z + 2*sin(z) +5;

x_c = [0, 0, 0]; % the current x in the iteration
x_l = [1, 1, 1]; % keep track of last x with a initial value of [x, y, z]
count = 0; % keep track of the number of iterations initialized as 1
error = [1, 1, 1]; % absolute relative approximate error percentage initialized as 1 for each row

% The Jacobian matrix 
Fx = [f1; f2; f3];
J = jacobian(Fx);

while (error(1,1) >= 0.1) || (error(1,2) >= 0.1) || (error(1,3) >= 0.1)
    % get Fx output for the iteration
    a = double(subs(f1, [x, y, z], x_l));
    b = double(subs(f2, [x, y, z], x_l));
    c = double(subs(f3, [x, y, z], x_l));
    Fx_o = [a; b; c];
    
    % get J output for the iteration
    d = double(subs(J(1,:), [x, y, z], x_l));
    e = double(subs(J(2,:), [x, y, z], x_l));   
    f = double(subs(J(3,:), [x, y, z], x_l));
    J_o = [d; e; f];
    
    % Perform Gauss Elimination to get delta x
    % forward elimination
    A = [J_o Fx_o]; % combine J_o, Fx_o together to do operation
    row_num = 3; % the number of rows
    row_swap = 1; % keep track of the row number that will be swapped with the cur_row
    ope_num = row_num - 1; % the number of operations in total
    cur_row = 2; % keep track of the current row of subtraction for this operation (initially at row 2)
    cur_col = 1; % keep track of the current column of subtraction (initially at column 1)
    count_t = cur_row; % keep track of the current row being manipulated

    for a = 1:ope_num
        count_t = cur_row;
        while count_t < (row_num+1)
            A(count_t,:) = A(count_t,:)-A(a,:).*(A(count_t,cur_col)/A(a,cur_col));
            count_t = count_t + 1;
        end
        % row swap
        [~,row_swap] = max(abs(A(cur_row:row_num,cur_col+1))); % find the row number with the max 
        temp = A(cur_row,:);
        A(cur_row,:) = A(row_swap+cur_row-1, :);
        A(row_swap+cur_row-1,:) = temp;
        cur_row = cur_row + 1;
        cur_col = cur_col + 1;
    end

    % back substitution
    delta_x = zeros (row_num,1); % solution matrix of delta x
    x_s = row_num; % x for solving from the bottom right diagonally
    num_sub = 0; % number of subtraction required for A(x, row_sum)
    cur_sub = num_sub; % keep track of needed subtraction left for A(x, row_sum)

    while x_s > 0
        while cur_sub > 0
            n = row_num - (cur_sub - 1);
            A(x_s,row_num+1) = A(x_s, row_num+1) -  delta_x(n, 1)*A(x_s, n);
            cur_sub = cur_sub - 1;
        end
        delta_x(x_s,1) = A(x_s,row_num+1)/A(x_s,x_s);
        num_sub = num_sub + 1;
        cur_sub = num_sub;
        x_s = x_s - 1;
    end
    % convert delta_x to its transpose for subs format
    delta_x = transpose(delta_x);
    x_c = x_l - delta_x;
    error(1,1) = abs((x_c(1,1) - x_l(1,1)) / x_c(1,1)) * 100;
    error(1,2) = abs((x_c(1,2) - x_l(1,2)) / x_c(1,2)) * 100;
    error(1,3) = abs((x_c(1,3) - x_l(1,3)) / x_c(1,3)) * 100;
    x_l = x_c;
    count = count + 1;
end

% convert x_c, error to their transpose for display format
x_c = transpose(x_c);
error = transpose(error);

disp("The final values obtained is");
disp(x_c);
disp("The number of required iterations is");
disp(count);
disp("The absolute relative approximate error % is");
disp(error);