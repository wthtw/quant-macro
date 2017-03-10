% basic graphing

x = [1:10];
y = x.^2;       % entry-wise squaring of elements in grid x plots out y=x^2
y2 = x;         % y = x for elements x in x
clf;            % clear current figures
plot(x, y, 'd-', x, y2, 'r-');       % plot y and y2 as functions of x
hold on;                            % retains plots in current axes
plot(x, y - 10, 'linewidth', 4);    % plots x, y - 10 using x, y above

clf;
plot(x, y - 10, 'linewidth', 4);    % plot x and y - 10 using x, y above
xlabel('time', 'fontsize', 20);     % label x axis, 'time' with custom font
ylabel('GDP');                      % label y axis
title('GDP over time');             % title the graph

% for loops

for i=1:10
    yx(i) = i - i^2;
    [i yx(i)];
    pause
end
yx;
% pause

% while loops
% example: convergence to a fixed point

dist = 1;
x(1) = 1;       % start with a guess for the solution
i = 1;
while dist > 0.000001   % while dist > some convergence criterion
    x(i + 1) = 1 + x(i)^0.5;    % update the guess
    [x(i) x(i + 1)]             % print initial and revised guess
    dist = abs(x(i) - x(i + 1)) % set new (shorter) distance
    i = i + 1
    pause
end
