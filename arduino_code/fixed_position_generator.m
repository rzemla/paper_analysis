
laps = 50;

every = 40;

x = 5:every:185;

x_rep = repmat(x,1,laps);

dlmwrite('every40.txt',x_rep)