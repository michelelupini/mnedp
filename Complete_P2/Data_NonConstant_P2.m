U = @(x,y) x.*y - x.*y.^2 + 1;
grad_U = @(x,y) [y-y.^2; x-2*x.*y];

nu = @(x,y) x;
beta = @(x,y) [x;1];
gamma = @(x,y) y;
f = @(x,y) 2*x.^2 - x.*y + x + y.^2 - x.*y.^3 ;
gD = @(x,y) 1;
gN = @(x,y) y - y.^2;