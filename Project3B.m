% File Name: Mark Ira Galang
% Date: 5/6/25
% Final Version: 5/15/25
% This is my way of implementing a rref solver

% Clean up stuff
clc;
clear;

% syms and stuff
syms x;

%% Part 1: Creating the rref method

function matrix = rref(matrix)
    % we need to go through some conditions
    % first we need a way to keep track of the row we are using to
    % eliminate, then use that row to eliminate other rows
    % then also check to see if the diagonal is non-zero
    % then also deal with zero rows
    % get the rows and cols
    [rows, cols] = size(matrix);
    % in the event that rows != cols, then we can only focus on the
    % diagonal by getting the minimum of the two

    % before anything, we need to deal with the zero pivots
    
    % run through each diagonal
    % we want to run through each diagonal so that we can do ref
    % check for the minimum as we want to stop before reaching out
    % of bounds
    k = 1;

    while k < min(rows, cols) + 1
        % get the current diagonal as this is the row we are using
        % to eliminate
        curr = matrix(k, k);
        % Skip if pivot is zero
        if abs(curr) < 1e-12
            % if it is zero, we need to swap the row with a pivot with
            % a non-zero row
            for row = k:rows
                % skip the row itself
                if k ~= row
                    % check if the row pivot is not zero
                    if matrix(row, k) ~= 0
                        % swap
                        temp = matrix(row, :);
                        matrix(row, :) = matrix(k, :);
                        matrix(k, :) = temp;
                        break;
                    end
                end
            end
            % restart operation by skipping everything below
            continue;
        end
        
        % Normalize the pivot row a.k.a make the pivot element equal to 1
        % it just makes life easier if we get the identity matrix
        matrix(k, :) = matrix(k, :) / curr;

        % Eliminate other rows 
        for row = 1:rows
            % make sure not to use the row we are using to eliminate
            if row ~= k
                % grab the coefficient needed to multiply and eliminate
                multiplier = matrix(row, k);
                % grab the row used to do an elimination
                r_eli = matrix(k, :);
                % grab the row we are usig to get a new row
                r_old = matrix(row, :);
                
                % Operation: R_old - Multi*R_eli = R_new 
                % rather than doing each number from each row
                % matlab has the ability to do subtraction with matricies 
                % making life a million times easier
                matrix(row, :) = round(r_old - multiplier * r_eli, 3);
            end
        end

        % increment k to move to the next diagonal
        k = k + 1;
    end
end

% function to make matrices look pretty :3
function printMatrix(matrix)
    [rows, cols] = size(matrix);
    
    % Define fixed width per entry
    % we are rounding anyway so by setting the size, everything will
    % hopefully be spaced evenly

    entry_width = 8; % e.g., 7.2f = width 7 + spacing
    sep_width = 3;   % space for " | "

    % Calculate total width: entries + separator + outer pipes
    total_width = cols * entry_width + sep_width + 1;

    % Print top line
    fprintf(' %s\n', repmat('-', 1, total_width));

    % then print matrix
    for row = 1:rows
        fprintf("|");
        for col = 1:cols - 1
            fprintf(" %7.2f", matrix(row, col));
        end
        fprintf(" | %8.2f |\n", matrix(row, end));
    end

    % Print bottom line
    fprintf(' %s\n', repmat('-', 1, total_width));
end

%% Part 2: Check if it works with given
A1_Q2 = [14, 14, -9, 3, -5; 
      14, 52, -15, 2, -32;
      -9, -15, 36, -5, 16;
      3, 2, -5, 47, 49;
      -5, -32, 16, 49, 79];

B1_Q2 = [-15; -100; 106; 329; 463];

disp("Original Matrix: ")
printMatrix([A1_Q2, B1_Q2])
disp("Solve Matrix: ")
printMatrix(rref([A1_Q2, B1_Q2])) 

% DO NOT COPY PASE LOL
% luckily we just needed to transpose
A2_Q2 = transpose([9, 7, 2, 0, 7;
      3, 6, 7, 9, 3;
      2, 9, 7, 7, 6; 
      0, 6, 8, 2, 4; 
      7, 4, 2, 2, 3]);
B2_Q2 = [35; 58; 53; 37; 39];

disp("Original Matrix: ")
printMatrix([A2_Q2, B2_Q2])
disp("Solved Matrix: ")
printMatrix(rref([A2_Q2, B2_Q2]))

%% Part 3: 10 random matricies... well sort of random

% The first three are taken from the linear algebra book
A1 = [4.5, -0.5;
        1, -3.5];
B1 = [1; -1];

disp("Original Matrix: ")
printMatrix([A1, B1])
disp("Solve Matrix: ")
printMatrix(rref([A1, B1]))

A2 = [20,   1, -1;
       1, -10,  1;
      -1,   1, 10];
B2 = [17; 13; 18];

disp("Original Matrix: ")
printMatrix([A2, B2])
disp("Solve Matrix: ")
printMatrix(rref([A2, B2]))

A3 = [3, -1,  0,  0;
     -1,  3, -1,  0;
      0, -1,  3, -1;
      0,  0, -1,  3];
B3 = [1; 0; 1; 1];

disp("Original Matrix: ")
printMatrix([A3, B3])
disp("Solve Matrix: ")
printMatrix(rref([A3, B3]))

% See if we can find the inverse matrix
A4 = [2,  3,  0; 
      1, -2, -1;
      2,  0, -1];

B4 = [1, 0, 0; 0, 1, 0; 0, 0, 1];

disp("Original Matrix: ")
printMatrix([A4, B4])
disp("Solve Matrix: ")
printMatrix(rref([A4, B4]))

% test for zero rows
A5 = [0, 0, 1; 
      2, 2, 4;
      0, 1, 3];

B5 = [1; 1; 1];

disp("Original Matrix: ")
printMatrix([A5, B5])
disp("Solve Matrix: ")
printMatrix(rref([A5, B5]))

%% The rest are randomized
% GUI menu with checkboxes to choose matrix systems
% Define the list of options
% I only added this because the anything greater than 50x50 is just too big
options = {
    '5x5', 
    '10x10', 
    '20x20', 
    '50x50', 
    '75x75', 
    '100x100'
};

% Display checkbox dialog to get user input
[selectedIndices, ok] = listdlg( ...
    'PromptString', 'Select which matrix sizes to solve:', ...
    'SelectionMode', 'multiple', ...
    'ListString', options);

% If user cancels or doesn't select anything
if ~ok || isempty(selectedIndices)
    disp('No selection made or dialog cancelled.');
else
    % Matrix size corresponding to each option
    matrix_sizes = [5, 10, 20, 50, 75, 100];

    % Loop through selected sizes and solve
    for idx = selectedIndices
        n = matrix_sizes(idx);
        fprintf("\n--- Solving %dx%d System ---\n", n, n);
        solve_random_system(n);
    end
end
% Function to create and solve random matrix systems
function matrices = solve_random_system(n) 
    non_zero_percent = 0.99; % this is the amount of non-zeros%
    factor = 101;            % this is to scale it to [1, 100]
    
    % first input it into matrices
    A = full(round(sprand(n, n, non_zero_percent) * factor)); 
    B = round(rand(n, 1) * factor);
    
    disp("Original Matrix: ")
    printMatrix([A, B])
    disp("Solve Matrix: ")
    printMatrix(rref([A, B]))
end

%% Bonus: compare my solver to MATLAB's backslash and LU factorization
% create a 500x500 matrix

A500 = full(round(sprand(500, 500, 1.0) * 501));  % 100% sparsity (all non-zero)
B500 = round(rand(500, 1) * 501);

disp("Check how fast each method are for a 500x500 Matrix")
% My Solver (RREF)
disp("My Solver: ")
tic;
sol_rref = rref([A500, B500]);  
toc;

% MATLAB Backslash Solver
disp("Backslash Solver: ")
tic;
x = A500 \ B500;
toc;

% LU Decomposition
disp("LU Decomposition: ")
tic;
[L, U, P] = lu(A500);
toc;

% Solving using LU
disp("LU Factorization:")
tic;
y = L \ (P * B500);
x_lu = U \ y;
toc;

% compare solutions
% RREF get last column
x_rref = sol_rref(:, end);

% Compare the results
disp("Norm of difference between RREF and backslash: ")
disp(norm(x_rref - x))

disp("Norm of difference between RREF and LU solution: ")
disp(norm(x_rref - x_lu))

disp("Norm of difference between backslash and LU solution: ")
disp(norm(x - x_lu))
