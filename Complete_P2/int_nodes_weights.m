% zita da buttare
% omega pesi
function [zita,csi,eta,omega] = int_nodes_weights(order,type)

% N1=csi cioè x,  il nodo uno e` (1,0)
% N2=eta cioè y,  il nodo uno e` (0,1)
% N3=zita, il nodo uno e` (0,0)
  
  persistent orders2D xw2D;
  
  try
    orders2D(order);
  catch me
    orders2D(order)=0;
  end
    
  if(orders2D(order)==1)
    csi   = xw2D(order).csi;
    eta   = xw2D(order).eta;
    zita  = xw2D(order).zita;
    omega = xw2D(order).omega;
  else
    if(order==5)
      xw2D(order).csi=[1/3,...
                       (6+sqrt(15))/21,...
                       (9-2*sqrt(15))/21,...
                       (6+sqrt(15))/21,...
                       (6-sqrt(15))/21,...
                       (9+2*sqrt(15))/21,...
                       (6-sqrt(15))/21];
      xw2D(order).eta=[1/3,...
                       (6+sqrt(15))/21,...
                       (6+sqrt(15))/21,...
                       (9-2*sqrt(15))/21,...
                       (6-sqrt(15))/21,...
                       (6-sqrt(15))/21,...
                       (9+2*sqrt(15))/21];
      xw2D(order).zita=1-xw2D(order).csi-xw2D(order).eta;

      xw2D(order).omega=[ 9/80,...
                          (155+sqrt(15))/2400,...
                          (155+sqrt(15))/2400,...
                          (155+sqrt(15))/2400,...
                          (155-sqrt(15))/2400,...
                          (155-sqrt(15))/2400,...
                          (155-sqrt(15))/2400];
    
      csi   = xw2D(order).csi;
      eta   = xw2D(order).eta;
      zita  = xw2D(order).zita;
      omega = xw2D(order).omega;
      orders2D(order)=1;
    end
  end

return