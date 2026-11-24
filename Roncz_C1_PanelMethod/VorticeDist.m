function [u_s,w_s,u_e,w_e] = VorticeDist(gamma_start,gamma_end,x,z,x_start,z_start,x_end,z_end)
% ControlPanel = 1 si el panel se influye a si mismo
   alpha=atan2((z_end-z_start),(x_start-x_end));

   MatrizGiro = [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];
   
       % Matriz de rotación
    MatrizGiro(abs(MatrizGiro) < 1e-14) = 0;  % limpiar casi-ceros
   
   Xp = MatrizGiro*[x-x_start;z-z_start];
   Xpj = MatrizGiro*[x_start-x_start;z_start-z_start];
   Xpj1 = MatrizGiro*[x_end-x_start;z_end-z_start];
   
   x_j_p = Xpj(1);
   % z_j_p = Xpj(2);

   x_p   = Xp(1);
   z_p   = Xp(2)+10^-5;

   x_j1_p = Xpj1(1);
   % z_j1_p = Xpj1(2);

   % gamma_0 = ((gamma_end-gamma_start)/x_j1_p);
   % 
   % u_p = gamma_start/(2*pi)*(atan2(z_p,(x_p-x_j1_p))-atan2(z_p,(x_p-x_j_p)))-...
   %     gamma_0/(4*pi)*(z_p*log(((x_p-x_j_p)^2+z_p^2)/...
   %     ((x_p-x_j1_p)^2+z_p^2))-2*x_p*(atan2(z_p,(x_p-x_j1_p))-atan2(z_p,(x_p-x_j_p))));
   % 
   % w_p = -gamma_start/(4*pi)*log(((x_p-x_j_p)^2+z_p^2)/((x_p-x_j1_p)^2+z_p^2))-...
   %     gamma_0/(2*pi)*(x_p/2*log(((x_p-x_j_p)^2+z_p^2)/((x_p-x_j1_p)^2+z_p^2))+...
   %     (x_j_p-x_j1_p)+z_p*(atan2(z_p,(x_p-x_j1_p))-atan2(z_p,(x_p-x_j_p))));

   gamma_0 = ((0-gamma_start)/x_j1_p);

   u_p_a = gamma_start/(2*pi)*(atan2(z_p,(x_p-x_j1_p))-atan2(z_p,(x_p-x_j_p)))-...
       gamma_0/(4*pi)*(z_p*log(((x_p-x_j_p)^2+z_p^2)/...
       ((x_p-x_j1_p)^2+z_p^2))-2*x_p*(atan2(z_p,(x_p-x_j1_p))-atan2(z_p,(x_p-x_j_p))));

   w_p_a = -gamma_start/(4*pi)*log(((x_p-x_j_p)^2+z_p^2)/((x_p-x_j1_p)^2+z_p^2))-...
       gamma_0/(2*pi)*(x_p/2*log(((x_p-x_j_p)^2+z_p^2)/((x_p-x_j1_p)^2+z_p^2))+...
       (x_j_p-x_j1_p)+z_p*(atan2(z_p,(x_p-x_j1_p))-atan2(z_p,(x_p-x_j_p))));
    
   gamma_0 = ((gamma_end-0)/x_j1_p);

   u_p_b = 0/(2*pi)*(atan2(z_p,(x_p-x_j1_p))-atan2(z_p,(x_p-x_j_p)))-...
       gamma_0/(4*pi)*(z_p*log(((x_p-x_j_p)^2+z_p^2)/...
       ((x_p-x_j1_p)^2+z_p^2))-2*x_p*(atan2(z_p,(x_p-x_j1_p))-atan2(z_p,(x_p-x_j_p))));

   w_p_b = -0/(4*pi)*log(((x_p-x_j_p)^2+z_p^2)/((x_p-x_j1_p)^2+z_p^2))-...
       gamma_0/(2*pi)*(x_p/2*log(((x_p-x_j_p)^2+z_p^2)/((x_p-x_j1_p)^2+z_p^2))+...
       (x_j_p-x_j1_p)+z_p*(atan2(z_p,(x_p-x_j1_p))-atan2(z_p,(x_p-x_j_p))));

   % V_p = [u_p;w_p];
   % V   = MatrizGiro'*V_p;

   % u = V(1);
   % w = V(2);

   V_p_a = -[u_p_a;w_p_a];
   V_a   = MatrizGiro'*V_p_a;

   u_s = V_a(1);
   w_s = V_a(2);

   V_p_b = -[u_p_b;w_p_b];
   V_b   = MatrizGiro'*V_p_b;

   u_e = V_b(1);
   w_e = V_b(2);
end