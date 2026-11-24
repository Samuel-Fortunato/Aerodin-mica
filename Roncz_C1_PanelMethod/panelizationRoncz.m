function [x_start,z_start,x_end,z_end,x_control,z_control,n_x,n_z,t_x,t_z,DeltaS] = panelizationRoncz(filename,N_deseado,mostrar_figura)
% Geometria Lednicer; N_deseado = numero total de paineis (par).
% Com apenas filename conserva os 70 paineis originais, sem interpolacao.
if nargin < 2 || isempty(N_deseado), N_deseado = 70; end
if nargin < 3 || isempty(mostrar_figura), mostrar_figura = true; end
if ~isscalar(N_deseado) || ~isfinite(N_deseado) || N_deseado < 4 || mod(N_deseado,2) ~= 0
    error('N_deseado deve ser um inteiro par >= 4.');
end
fid = fopen(filename,'r');
if fid == -1, error('Nao foi possivel abrir %s',filename); end
fgetl(fid);
counts = sscanf(fgetl(fid),'%f');
if numel(counts) < 2, fclose(fid); error('Cabecalho Lednicer invalido.'); end
upper = fscanf(fid,'%f %f',[2,counts(1)])';
lower = fscanf(fid,'%f %f',[2,counts(2)])';
fclose(fid);
if size(upper,1) ~= counts(1) || size(lower,1) ~= counts(2)
    error('Faltam coordenadas no ficheiro Lednicer.');
end
% Contorno anti-horario: bordo de fuga superior -> bordo de ataque -> fuga inferior
xu = flipud(upper(:,1)); zu = flipud(upper(:,2));
xl = lower(:,1);       zl = lower(:,2);
% 36 nos por superficie = 35 paineis por superficie = 70 no total.
N_original = size(upper,1) + size(lower,1) - 2;
if N_deseado ~= N_original
    su = [0; cumsum(hypot(diff(xu),diff(zu)))]; su = su / su(end);
    sl = [0; cumsum(hypot(diff(xl),diff(zl)))]; sl = sl / sl(end);
    qu = linspace(0,1,N_deseado/2+1)';
    ql = linspace(0,1,N_deseado/2+1)';
    xu0=xu; zu0=zu; xl0=xl; zl0=zl;
    xu=interp1(su,xu0,qu,'pchip'); zu=interp1(su,zu0,qu,'pchip');
    xl=interp1(sl,xl0,ql,'pchip'); zl=interp1(sl,zl0,ql,'pchip');
end
xraw = [xu;xl(2:end)]; zraw = [zu;zl(2:end)];
c = max(xraw)-min(xraw);
x = (xraw-min(xraw))/c; z = zraw/c;
x_start = x(1:end-1); z_start = z(1:end-1);
x_end = x(2:end); z_end = z(2:end);
x_control=(x_start+x_end)/2; z_control=(z_start+z_end)/2;
dx=x_end-x_start; dz=z_end-z_start;
DeltaS=hypot(dx,dz);
if any(DeltaS<=0), error('Painel de comprimento zero.'); end
t_x=dx./DeltaS; t_z=dz./DeltaS;
n_x=t_z; n_z=-t_x;
if mostrar_figura
    figure('Name','Panelizacao Roncz','Color','w');
    plot(x,z,'k-o','LineWidth',1.5,'MarkerFaceColor','r','MarkerSize',4); hold on;
    plot(x_control,z_control,'ks','MarkerSize',4,'MarkerFaceColor','k');
    axis equal; grid on; box on;
    xlabel('$x/c\ [-]$','Interpreter','latex');
    ylabel('$z/c\ [-]$','Interpreter','latex');
    title(sprintf('Panelizacao Roncz (%d paineis)',numel(DeltaS)),'Interpreter','none');
    legend('Nos','Pontos de controlo','Location','best','Interpreter','none');
    hold off;
end
end
