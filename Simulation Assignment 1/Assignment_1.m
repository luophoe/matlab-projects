format shortg
clc;
A = readmatrix("A.txt");
B = readmatrix("B.txt");
row_num = size(A, 1); % total number of rows


% Simple Matrix Inversion
disp("****** Simple Matrix Inversion ******")
S1 = inv(A)*B; % solution matrix
S1 = round(S1,5,'significant');
disp(S1);


% Gaussian Elimination with Partial Pivoting
disp("****** Gaussian Elimination with Partial Pivoting ******")
% forward elimination
A2 = [A B]; % combine A, B together to do operation
row_swap = 1; % keep track of the row number that will be swapped with the cur_row
ope_num = row_num - 1; % the number of operations in total
cur_row = 2; % keep track of the current row of subtraction for this operation (initially at row 2)
cur_col = 1; % keep track of the current column of subtraction (initially at column 1)
count = cur_row; % keep track of the current row being manipulated

for a = 1:ope_num
    count = cur_row;
    while count < (row_num+1)
          A2(count,:) = A2(count,:)-A2(a,:).*(A2(count,cur_col)/A2(a,cur_col));
          count = count + 1;
    end
    % row swap
    [~,row_swap] = max(abs(A2(cur_row:row_num,cur_col+1)));
    temp = A2(cur_row,:);
    A2(cur_row,:) = A2(row_swap+cur_row-1, :);
    A2(row_swap+cur_row-1,:) = temp;
    cur_row = cur_row + 1;
    cur_col = cur_col + 1;
end

% back substitution
S2 = zeros (row_num,1); % solution matrix
x = row_num; % x for solving from the bottom right diagonally
num_sub = 0; % number of subtraction required for A(x, row_sum)
cur_sub = num_sub; % keep track of needed subtraction left for A(x, row_sum)

while x > 0
    while cur_sub > 0
        n = row_num - (cur_sub - 1);
        A2(x,row_num+1) = A2(x, row_num+1) -  S2(n, 1)*A2(x, n);
        cur_sub = cur_sub - 1;
    end
    S2(x,1) = A2(x,row_num+1)/A2(x,x);
    num_sub = num_sub + 1;
    cur_sub = num_sub;
    x = x - 1;
end
S2 = round(S2,5,'significant');
disp(S2);


% Gauss-Seidel Iteration
disp("****** Gauss-Seidel Iteration ******")
% 1% error
S3_cur = zeros(row_num,1); % store values in the current operation
S3_last = zeros(row_num,1); % store values in the last operation
ite_req = 0; % keep track of number of iterations required
error_req = 1; % error required
stop = 0; % change to 1 when the loop can be stopped

while stop == 0
    ite_req = ite_req + 1;
    for b = 1:row_num
        x = B(b,1);
        for c = 1:row_num
            if c ~= b
                x = x - A(b,c)*S3_cur(c,1);
            end
        end
        x = x/A(b,b);
        x = round(x,5,'significant');
        S3_cur(b,1) = x;
    end
    % check for error
    count_2 = 0; % count if all the elements fits error required
    for d = 1:row_num
        error = abs((S3_cur(d,1) - S3_last(d,1))/S3_cur(d,1)*100);
        if error < error_req
            count_2 = count_2 + 1;
        end
        if count_2 == row_num
            stop = 1;
        end
    end
    S3_last = S3_cur;
   
end

disp(S3_cur);
disp(ite_req);