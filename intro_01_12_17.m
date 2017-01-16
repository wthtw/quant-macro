clear
x = 2;
x = x + 1;

v1 = [1 4 3 2];

v2 = [ 2 1 4 1]';   % transpose

m1 = [1 3 4 1; 1 3 1 1]; % 2x4 matrix

m2 = [ [1 3 4 1]' ; [1 2 1 1]' ];    % 8x1 column vector

v3 = 0:0.2:1;   % gives a row vector of #'s from 0 to 1 by 0.2 increments

v4 = [0:0.2:1]; % same as above

v5 = [0:0.2:1]';    % the above's transpose

% important matrices
v6 = zeros(1,6);    % 1x6 zero row vector
v7 = zeros(6,1);    % 6x1 zero column vector

m4 = zeros(2,6);     % 2x6 zero matrix
m5 = ones(2,3);      % 2x3 matrix of 1's

m1;
m6 = m1 + 2;     % element by element addition: add 2 to each entry of m1

m1;
m6;
m7 = m1 + m6;       % matrix addition of m1 and m6

m8 = m1 - m6;       % matrix subtraction of m6 from m1

v2;
m6;
v8 = m6 * v2;       % product of 2x4 matrix and 4x1 vector

m1;
m6;
m9 = m1 .* m6;      % entry-wise product of 2x4 matrix and 2x4 matrix

m1;
m6;
m10 = m1 ./ m6;     % entry-wise quotient of 2x4 matrix and 2x4 matrix

m1;
m10 = m1 .^2;       % entry-wise exponentiation of 2x4 matrix elements

m11 = [1 2; 2 3];
m12 = m11^2;         % squaring a square matrix (via matrix multiplication)

m11;
m13 = m11.^2;        % entry-wise squaring of m11's elements

% simple variables

x2 = 2 * 3 + 1 - 2 / 3;      % x2 = 6.3333

x2 = 2 * (3 + 1) - 2 / 3;        % x2 = 7.3333

size(m6);                % gives ans = 2 4
[x y] = size(m6);        % [x y] is assigned to [2 4]
                         % this is useful to avoid mistakes involving
                         % mismatched dimensions of vectors and matrices.

x;
y;
m14 = zeros(x,y);        % since x = 2, y = 4, zeros(x,y) creates a 2x4
                         % matrix of zeros.
m14 = zeros(size(m6));   % does the same as above.

v1;
x4 = sum(v1);            % sums all the entries of vector v1. x4 = 10

m6;
v9 = sum(m6);            % sums all the columns of m6 to a 1x4 row vector
v10 = sum(m6,2);         % sums all the rows of m6 to a 2x1 column vector
                         % note that sum(m6) == sum(m6,1)

m6;
m10 = exp(m6);           % entry-wise exp() of entries in m6, 2x4 matrix.
m11 = log(m6);           % entry-wise log() of entries in m6, 2x4 matrix.


% commonly used technique below (maximizing functions):

v11 = [ 1 4 5 10 2 55 2 4 5 ];
max(v11);            % finds the largest value in the vector v11

[x5 ix5] = max(v11); % gives x5 = 55, ix5 = 6 (index of 55 in v11)
                     % used frequently... [max of obj fn, index of grid]

% maximization example: very commonly used in macro

x = [0:0.1:1];
y = x - x.^2;
[ymax iymax] = max(y);   % ymax = 0.2500, iymax = 6. gives maximum of y and 
                         % which point in the x g rid where it is maximized.
xymax = x(iymax);        % gives the value of x which maximizes y using
                         % the position of y, iymax.

% maximization example: where maximizing point is not in the grid space

x2 = [0:0.01:1];             % we change from 0.1 to 0.01.
y2 = x2 - x2.^1.3;           % change from 2 to 1.3
[y2max iy2max] = max(y2);    % [0.0962 43]
xy2max = x2(iy2max);         % xy2max = 42

% how fine should our grids be? if the statistic we're computing doesn't
% change, we've found a good enough grid.

% extra notes:
% 'save myvariables v5 v6 x4' saves the chosen variables in a file
% myvariables.mat
% you can clear the memory, then load the file myvariables.mat to get the
% three variables, using 'load myvariables'
% you can also load from a text file. the numbers should fit either a
% matrix format or a vector format.
% command 'clear' clears the memory. usually we start a matlab file with
% 'clear' on the first line.