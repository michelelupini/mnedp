function [b_N] = neumann_nonconstant_P2(geom,gN)

Ndof = max(geom.pivot.pivot);
b_N = zeros(Ndof,1);
for e = 1:length(geom.pivot.Ne(:,1))
  l = geom.pivot.Ne(e,1);
  i_b = geom.elements.borders(l,1);
  i_e = geom.elements.borders(l,2);
  i_m = geom.elements.borders(l,5);
  iib = geom.pivot.pivot(i_b);
  iie = geom.pivot.pivot(i_e);
  iim = geom.pivot.pivot(i_m);
  Vb = geom.elements.coordinates(i_b,:);
  Ve = geom.elements.coordinates(i_e,:);
  normE = norm(Ve - Vb);
  
  x_t = @(t) Vb(1) + (Ve(1)-Vb(1)).*t;
  y_t = @(t) Vb(2) + (Ve(2)-Vb(2)).*t;
  syms x
  
  if iib > 0
    % b_N(iib) = b_N(iib) + normE*(2*gN(Vb(1),Vb(2))+gN(Ve(1),Ve(2)))/6;
    % b_N(iib) = b_N(iib) + normE*integral(@(x) gN(x_t(x),y_t(x))*2.*(1-x).*(0.5-x),0,1);
    b_N(iib) = b_N(iib) + normE*vpaintegral(gN(x_t(x),y_t(x))*2.*(1-x).*(0.5-x),0,1);
  end
  if iie > 0
    % b_N(iie) = b_N(iie) + normE*(gN(Vb(1),Vb(2))+2*gN(Ve(1),Ve(2)))/6;
    % b_N(iie) = b_N(iie) + normE*integral(@(x) gN(x_t(x),y_t(x))*2.*x.*(x-0.5),0,1);
    b_N(iie) = b_N(iie) + normE*vpaintegral(gN(x_t(x),y_t(x))*2.*x.*(x-0.5),0,1);
  end
  if iim > 0 
    % b_N(iim) = b_N(iim) + normE*integral(@(x) gN(x_t(x),y_t(x))*4.*x.*(1-x),0,1);
    b_N(iim) = b_N(iim) + normE*vpaintegral(gN(x_t(x),y_t(x))*4.*x.*(1-x),0,1);
  end
end

end