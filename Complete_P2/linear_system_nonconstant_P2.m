function [A,b,A_D,u_D] = linear_system_nonconstant_P2(geom,nu,beta,gamma,f,gD)

global csi eta omega;
global Hat_Phi_Matrix Hat_Phi_Tensor;

deltaX = zeros(3,1);
deltaY = zeros(3,1);

pivot = geom.pivot.pivot;
nT = geom.nelements.nTriangles;

Ndof = max(pivot);
N_D = size(geom.pivot.Di,1);
A = zeros(Ndof,Ndof);
b = zeros(Ndof,1);
A_D = zeros(Ndof,N_D);
u_D = zeros(N_D,1);

for e = 1:nT
  ele = geom.elements.triangles(e,:);
  X = geom.elements.coordinates(ele,1);
  Y = geom.elements.coordinates(ele,2);
  
  deltaX(1) = X(3) - X(2);
  deltaX(2) = X(1) - X(3);
  deltaX(3) = X(2) - X(1);
  deltaY(1) = Y(2) - Y(3);
  deltaY(2) = Y(3) - Y(1);
  deltaY(3) = Y(1) - Y(2);
  
  % Xb = geom.support.TInfo(e).CG(1);
  % Yb = geom.support.TInfo(e).CG(2);
  
  AreaT = geom.support.TInfo(e).Area;
  
  B = [deltaX(2), -deltaX(1); -deltaY(2), deltaY(1)];
  B_inv = (1/(2*AreaT)) * [deltaY(1), deltaX(1); deltaY(2), deltaX(2)];
  B_invT = (B_inv)';
  
  for j = 1:6
    jj = pivot(ele(j));
    if jj > 0
      for k = 1:6
        kk = pivot(ele(k));
        if kk > 0
          for q = 1:length(omega)
            Fe_vector = [X(3), Y(3)]' + B*[csi(q); eta(q)];
            A(jj,kk) = A(jj,kk) ...
                + 2*AreaT*omega(q)*( nu(Fe_vector(1),Fe_vector(2))*Hat_Phi_Tensor(:,k,q)'*B_inv*B_invT*Hat_Phi_Tensor(:,j,q) ...
                + beta(Fe_vector(1),Fe_vector(2))'*B_invT*Hat_Phi_Tensor(:,k,q)*Hat_Phi_Matrix(j,q) ...
                + gamma(Fe_vector(1),Fe_vector(2))*Hat_Phi_Matrix(k,q)*Hat_Phi_Matrix(j,q) );
          end
          % A(jj,kk) = A(jj,kk) ...
          % + 0.25*nu*(deltaY(k)*deltaY(j)+deltaX(k)*deltaX(j))/AreaT ...
          % + (beta(1)*deltaY(k)+beta(2)*deltaX(k))/6 ...
          % + gamma*AreaT*(1+(j==k))/12;
        else
          for q = 1:length(omega)
            Fe_vector = [X(3), Y(3)]' + B*[csi(q); eta(q)];
            A_D(jj,-kk) = A_D(jj,-kk) ...
                + 2*AreaT*omega(q)*( nu(Fe_vector(1),Fe_vector(2))*Hat_Phi_Tensor(:,k,q)'*B_inv*B_invT*Hat_Phi_Tensor(:,j,q) ...
                + beta(Fe_vector(1),Fe_vector(2))'*B_invT*Hat_Phi_Tensor(:,k,q)*Hat_Phi_Matrix(j,q) ...
                + gamma(Fe_vector(1),Fe_vector(2))*Hat_Phi_Matrix(k,q)*Hat_Phi_Matrix(j,q) );
          end
          % A_D(jj,-kk) = A_D(jj,-kk) ...
          % + 0.25*nu*(deltaY(k)*deltaY(j)+deltaX(k)*deltaX(j))/AreaT ...
          % + (beta(1)*deltaY(k)+beta(2)*deltaX(k))/6 ...
          % + gamma*AreaT*(1+(j==k))/12;
        end
      end
      for q = 1:length(omega) 
        Fe_vector = [X(3), Y(3)]' + B*[csi(q); eta(q)];
        b(jj) = b(jj) + 2*AreaT*omega(q)*f(Fe_vector(1),Fe_vector(2))*Hat_Phi_Matrix(j,q);
      end
      %b(jj) = b(jj) + f(Xb,Yb)*AreaT/3;
    end
  end
end

for i = 1:N_D
  u_D(i) = gD(geom.pivot.Di(i,1),geom.pivot.Di(i,2));
end

end