clear all

global zita csi eta omega;
[zita,csi,eta,omega] = int_nodes_weights(5);

global Hat_Phi Hat_Phi_Gradient Hat_Phi_Matrix Hat_Phi_Tensor;
Hat_Phi = @(x,y) [
    2*x.*(x-0.5); 
    2*y.*(y-0.5); 
    2*(1-x-y).*(0.5-x-y);
    4*x.*(1-x-y);
    4*y.*x
    4*y.*(1-x-y)
    ]; 
Hat_Phi_Gradient = @(x,y) [
    4*x-1,  0,      4*x+4*y-3,  4-8*x-4*y,    4*y,  -4*y;
    0,      4*y-1,  4*x+4*y-3,  -4*x          4*x,  4-4*x-8*y
    ];
Hat_Phi_Matrix = zeros(6,length(omega));
for q = 1:length(omega)
  Hat_Phi_Matrix(:,q) = Hat_Phi(csi(q),eta(q));
end
Hat_Phi_Tensor = zeros(2,6,length(omega));
for q = 1:length(omega)
    Hat_Phi_Tensor(:,:,q) = Hat_Phi_Gradient(csi(q),eta(q));
end

% Refining_Area = [0.02 0.005 0.002 0.001];
% Refining_Area = [0.005 0.001 0.0005 0.0001];
%Refining_Area = 0.0001; 
% Refining_Area = 0.0002;
Refining_Area = [0.005 0.0025 0.0012 0.001];
% Refining_Area = 0.02;

N_iter = length(Refining_Area);
E_inf = zeros(N_iter,1); % norm(u-u_h,inf)
E0 = zeros(N_iter,1); % norm(u-u_h,2)
E1 = zeros(N_iter,1); % norm(grad(u)-grad(u_h),2) 

condition_number_2 = zeros(N_iter,1);
condition_number_inf = zeros(N_iter,1);


%----------------------------------------------------------------------------

for iter = 1:N_iter

Sample_Square_Dirichlet_Neumann_P2
P2
Data_NonConstant_P2

% global ele X Y nT;
elements = geom.elements.triangles(:,:);
XY = geom.elements.coordinates(:,:);
nT = geom.nelements.nTriangles;
% Nnode = geom.nelements.nVertexes;
Nnode = length(geom.pivot.pivot);

u = zeros(Nnode,1);
U_ex = U(XY(:,1),XY(:,2));

[A,b,A_D,u_D] = linear_system_nonconstant_P2(geom,nu,beta,gamma,f,gD);
b_N = neumann_nonconstant_P2(geom,gN);
u0 = A\(b - A_D*u_D + b_N);

for j = 1:Nnode
  jj = geom.pivot.pivot(j);
  if jj > 0
    u(j) = u0(jj);
  else
    u(j) = u_D(-jj);
  end
end

condition_number_2(iter) = cond(A,2);
condition_number_inf(iter) = cond(A,inf);

elements_trisurf = [
    elements(:,1), elements(:,5), elements(:,4);
    elements(:,5), elements(:,2), elements(:,6)
    elements(:,6), elements(:,3), elements(:,4)
    elements(:,4), elements(:,5), elements(:,6)
    ];
figure
trisurf(elements_trisurf,XY(:,1),XY(:,2),U_ex)
figure
trisurf(elements_trisurf,XY(:,1),XY(:,2),u)
figure
trisurf(elements_trisurf,XY(:,1),XY(:,2),U_ex-u)

for e = 1:nT
  ele = elements(e,:);
  X = XY(ele,1);
  Y = XY(ele,2);
  AreaT = geom.support.TInfo(e).Area;
  
  deltaX(1) = X(3) - X(2);
  deltaX(2) = X(1) - X(3);
  deltaX(3) = X(2) - X(1);
  deltaY(1) = Y(2) - Y(3);
  deltaY(2) = Y(3) - Y(1);
  deltaY(3) = Y(1) - Y(2);
  
  B = [deltaX(2), -deltaX(1); -deltaY(2), deltaY(1)];
  B_inv = (1/(2*AreaT)) * [deltaY(1), deltaX(1); deltaY(2), deltaX(2)];
  B_invT = (B_inv)';
  
  for q = 1:length(omega)
    Fe_vector = [X(3), Y(3)]' + B*[csi(q); eta(q)];
    uh_E = 0;
    grad_uh_E = 0;
    for k = 1:6
      uh_E = uh_E + u(ele(k))*Hat_Phi_Matrix(k,q);
      grad_uh_E = grad_uh_E + u(ele(k))*Hat_Phi_Tensor(:,k,q);
    end
    E0(iter) = E0(iter) + 2*AreaT*omega(q)*( U(Fe_vector(1),Fe_vector(2)) - uh_E )^2;
    E1(iter) = E1(iter) + 2*AreaT*omega(q)*( grad_U(Fe_vector(1),Fe_vector(2)) - B_invT*grad_uh_E )' ...
        * ( grad_U(Fe_vector(1),Fe_vector(2)) - B_invT*grad_uh_E );
  end
end

E0(iter) = sqrt(E0(iter));
E1(iter) = sqrt(E1(iter));
E_inf(iter) = norm(U_ex-u,inf);

clear BC Domain geom RefiningOpions;

end

clear AreaT B B_inv B_invT deltaX deltaY e ele Fe_vector grad_uh_E jj k uh_E X Y

%----------------------------------------------------------------------------


p0 = polyfit(log(Refining_Area),log(E0),1);
p1 = polyfit(log(Refining_Area),log(E1),1);
p_inf = polyfit(log(Refining_Area),log(E_inf),1);