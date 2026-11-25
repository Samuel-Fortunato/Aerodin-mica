clearvars -except panelModeOverride;
close all; clc;
%% C1: escolha do modo
% normal: analise original com 70 paineis
% conver: convergencia de CL a 5 graus
% angles: CL(alpha) de -20 a 20 graus
% casos: Cp minimo, estagnacao aproximada e linhas de corrente a 0, 5 e 20 graus

filename = '../airfoil-points.txt';
V_inf_modul = 1;
N_final = 260; % atualizar depois do estudo de convergencia

fprintf('\n===== METODO DOS PAINEIS - RONCZ =====\n');
fprintf('1 - Analise normal\n');
fprintf('2 - Estudo de convergencia\n');
fprintf('3 - CL vs angulo de ataque\n');
fprintf('4 - Casos especificos\n');

if exist('panelModeOverride', 'var')
    opcao = panelModeOverride;
    clear panelModeOverride;
else
    opcao = input('Selecione uma opcao (1-4): ');
end

switch opcao
    case 1
        modo = 'normal';
    case 2
        modo = 'conver';
    case 3
        modo = 'angles';
    case 4
        modo = 'casos';
    otherwise
        error('Opcao invalida. Introduza um numero de 1 a 4.');
end

switch lower(modo)
    case 'conver'
        N_testes = 20:30:800;
        CL_testes = nan(size(N_testes)); erro_testes = nan(size(N_testes));
        tempos = nan(size(N_testes));
        for k = 1:numel(N_testes)
            tic;
            r = resolver_paineis_C1(filename,N_testes(k),5,V_inf_modul);
            CL_testes(k)=r.CL; erro_testes(k)=r.erro_normal;
            tempos(k)=toc;
        end

        CL_ref = CL_testes(end);

        Erro_ref_pct = abs((CL_testes - CL_ref)./CL_ref)*100;

        Tabela = table(N_testes(:), CL_testes(:), Erro_ref_pct(:), ...
            'VariableNames', {'N_paineis','CL','Erro_ref_pct'});

        disp(Tabela);
        
        variacao = [NaN,100*abs(diff(CL_testes))./max(abs(CL_testes(2:end)),eps)];
        disp(table(N_testes',CL_testes',variacao',erro_testes',tempos', ...
             'VariableNames',{'N','CL','Variacao_pct','ErroNormalRel','Tempo_s'}));


        %--------------------------------------

        figure('Name','Erro relativo de convergencia','Color','w');

        semilogy(N_testes, Erro_ref_pct, 'ko-', ...
            'LineWidth',1.5, ...
            'MarkerFaceColor','k', ...
            'MarkerSize',5);

        hold on;

        % Criterio de convergencia: 0.5 %
        yline(0.5,'r--', ...
            'Tolerancia = 0.5%', ...
            'LineWidth',1.5);

        % Numero selecionado
        xline(260,'b--', ...
            'N = 260', ...
            'LineWidth',1.5);

        xlabel('Numero de paineis, N');
        ylabel('Erro relativo de C_L [%]');

        title('Estudo de convergencia - Perfil Roncz', ...
            'Interpreter','none');

        grid on;
        box on;

        set(gca,'FontName','Times New Roman','FontSize',14);
        
        %---------------
        figure('Name','Convergencia CL','Color','w');
        plot(N_testes,CL_testes,'ko-','MarkerFaceColor','k','LineWidth',1.5);
        grid on; box on; xlabel('Numero de paineis, N');
        ylabel('$C_L$','Interpreter','latex');
        title('Convergencia de C_L para alpha = 5 deg','Interpreter','none');
        return
    case 'angles'
        angulos = -20:5:20;
        CL_a = zeros(size(angulos)); tempo_a=zeros(size(angulos));
        for k=1:numel(angulos)
            tic; r = resolver_paineis_C1(filename,N_final,angulos(k),V_inf_modul);
            CL_a(k)=r.CL; tempo_a(k)=toc;
        end
        %........................

        %% Comparacao CL: Metodo dos Paineis vs CFD

        AoA_CFD = [0 5 10 15 20];

        CL_CFD = [0.0594616944
            0.5976379630
            1.0645298200
            1.3789810000
            0.9559518660];

        figure('Name','Comparacao CL - Paineis e CFD','Color','w');

        hold on;

        h1 = plot(angulos,CL_a,'b-o', ...
            'LineWidth',1.8, ...
            'MarkerSize',6);

        h2 = plot(AoA_CFD,CL_CFD,'r-s', ...
            'LineWidth',1.8, ...
            'MarkerSize',7);

        xlabel('$\alpha\ [^\circ]$','Interpreter','latex');
        ylabel('$C_L\ [-]$','Interpreter','latex');

        legend([h1,h2], ...
            {'Metodo dos Paineis','CFD'}, ...
            'Location','northwest');

        grid on;
        box on;

        set(gca,'FontName','Times New Roman','FontSize',15);

        hold off;
        %...............


        disp(table(angulos',CL_a',tempo_a','VariableNames',{'AoA_graus','CL','Tempo_s'}));
        figure('Name','CL versus AoA','Color','w');
        plot(angulos,CL_a,'bo-','LineWidth',1.6,'MarkerFaceColor','b');
        xlabel('$\alpha\ [^\circ]$','Interpreter','latex');
        ylabel('$C_L$','Interpreter','latex');grid on;box on;
        idx = angulos >= -5 & angulos <= 5;
        p = polyfit(angulos(idx),CL_a(idx),1);
        fprintf('Declive local = %.6g /grau = %.6g /rad\n',p(1),p(1)*180/pi);
        fprintf('Alpha de sustentacao nula (ajuste) = %.3f graus\n',-p(2)/p(1));
        return
    case 'casos'
        for a0=[0 5 20]
            r=resolver_paineis_C1(filename,N_final,a0,V_inf_modul);
            %---------

            fprintf(['AoA=%g deg: CL=%.6f, Cp_min=%.6f, ' ...
                'x_Cp_min/c=%.5f, x_est/c=%.5f, z_est/c=%.5f\n'], ...
                    a0,r.CL,r.Cp_min,r.x_Cp_min,r.x_est,r.z_est);
            %--------------
            
            desenhar_linhas_corrente_sustentador(r.gamma,r.x_start,r.z_start, ...
                r.x_end,r.z_end,V_inf_modul,deg2rad(a0));
            %--------------
            figure('Name',sprintf('Estagnacao - AoA %.1f',a0), ...
                'Color','w');

            plot([r.x_start; r.x_end(end)], ...
                [r.z_start; r.z_end(end)], ...
                'k-','LineWidth',1.5);

            hold on;

            if isfinite(r.x_est)
                plot(r.x_est,r.z_est,'ro', ...
                    'MarkerSize',8, ...
                    'MarkerFaceColor','r');
            end

            axis equal;
            grid on;

            xlabel('$x/c$','Interpreter','latex');
            ylabel('$z/c$','Interpreter','latex');

            title(sprintf('Ponto de estagnacao - AoA %.1f deg',a0), ...
                'Interpreter','none');

            legend('Perfil Roncz','Estagnacao estimada', ...
                'Location','best');

            hold off;
            %......


            


        end
 

        return
    case 'normal'
        % Executar abaixo o codigo original
    otherwise
        error('Modo invalido. Escolha normal, conver, angles ou casos.');
end

%% Analise normal (codigo anterior conservado)
x = input('Angulo de ataque (graus): ');
AoA_rad = deg2rad(x);
[x_start,z_start,x_end,z_end,x_control,z_control,n_x,n_z,t_x,t_z, DeltaS] = ...
    panelizationRoncz(filename,N_final,true);
graficar_panelizationVectors_Roncz(x_start, z_start, x_end, z_end, ...
    x_control, z_control, ...
    t_x, t_z, n_x, n_z)

%% Passo 2 e 3: Aplicação VorticeDist e matriz a 
N = length(DeltaS); % # PANEIS

a=zeros(N+1,N+1);

for i=1:N 
for j=1:N
    [u_start,w_start, ~, ~] = VorticeDist(1,0,x_control(i),z_control(i), x_start(j), z_start(j),x_end(j),z_end(j));
    [~, ~,u_end,w_end] = VorticeDist(0,1,x_control(i), z_control(i), x_start(j), z_start(j),x_end(j),z_end(j));
    
    a(i,j) = a(i,j) + u_start* n_x(i) + w_start * n_z(i);

    a(i,j+1) = a(i,j+1) + u_end * n_x(i) + w_end*n_z(i);
end
end
a(N+1,1) =1;
a(N+1,N+1)=1;

FNSC = zeros(N+1,1);

[u_inf, w_inf] = EsUniforme(V_inf_modul, AoA_rad);
for i=1:N
    FNSC(i,1)= -[u_inf w_inf] * [n_x(i)  n_z(i)]';
end

FNSC(N+1,1) = 0;

gamma = a \ FNSC;

fig = desenhar_linhas_corrente_sustentador(gamma, x_start, z_start, x_end, z_end, V_inf_modul, AoA_rad);

%% Passo 4: cálculo dos coeficientes de pressão
V_velocidades = zeros(N,1);
V_normal = zeros(N,1);
V_tang = zeros(N,1);
Cp = zeros(N,1);

for i=1:N
   u_total = u_inf;
   w_total = w_inf;
   for  j=1:N
    [u_start, w_start, ~, ~] = VorticeDist(1,0,x_control(i), z_control(i), x_start(j),z_start(j), x_end(j), z_end(j));
    [~, ~, u_end, w_end] = VorticeDist(0,1,x_control(i),z_control(i), ...
        x_start(j),z_start(j),x_end(j),z_end(j));
    u_total = u_total + gamma(j)*u_start + gamma(j+1)*u_end;
    w_total = w_total + gamma(j)*w_start + gamma(j+1)*w_end;
   end
V_velocidades(i) = sqrt(u_total^2 + w_total^2);
V_tang(i) = u_total*t_x(i) + w_total*t_z(i);
V_normal(i) = u_total*n_x(i) + w_total*n_z(i);
end

V_inf =sqrt(u_inf^2 + w_inf^2);

Cp = 1 - (V_tang./V_inf).^2;

% Cp com modulo de velocidade
Cp_total = 1 - (V_velocidades./V_inf).^2;

% Cp com velocidade tangencial
Cp_tang = 1 - (V_tang./V_inf).^2;

% Diferencia entre os resultados
erro_Cp = max(abs(Cp_total - Cp_tang));

erro_normal = max(abs(V_normal));

fprintf('Error maximo de impermeabilidad: %.6e\n', erro_normal);
fprintf('Diferencia maxima entre Cp: %.6e\n', erro_Cp);

%% Grafica do Cp - Perfil Roncz

% Separar extradós e intradós
idx_upper = z_control >= 0;
idx_lower = z_control < 0;

x_upper = x_control(idx_upper);
Cp_upper = Cp(idx_upper);

x_lower = x_control(idx_lower);
Cp_lower = Cp(idx_lower);

% Ordenar de borde de ataque a borde de fuga
[x_upper, iu] = sort(x_upper);
Cp_upper = Cp_upper(iu);

[x_lower, il] = sort(x_lower);
Cp_lower = Cp_lower(il);

figure('Name','Distribuição do coeficiente de pressão C_p ao longo do perfil','Color','w');
hold on;
grid on;
box on;

% Extradós (azul)
h1 = plot(x_upper, -Cp_upper, 'o-', ...
    'Color', 'b', ...
    'LineWidth', 1.8, ...
    'MarkerSize', 5, ...
    'MarkerFaceColor', 'none', ...
    'DisplayName', 'Extradós');

% Intradós (rojo)
h2 = plot(x_lower, -Cp_lower, 'o-', ...
    'Color', 'r', ...
    'LineWidth', 1.8, ...
    'MarkerSize', 5, ...
    'MarkerFaceColor', 'none', ...
    'DisplayName', 'Intradós');

% Etiquetas
xlabel('$x/c\ [-]$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$-C_p\ [-]$', 'Interpreter', 'latex', 'FontSize', 18);

% Título
title('Distribuição de C_p em x/c - Perfil Roncz', ...
    'Interpreter', 'none', 'FontSize', 18, 'FontWeight', 'bold');

% Texto del angulo de ataque arriba a la derecha
text(0.97, 0.95, sprintf('$\\alpha = %.1f^{\\circ}$', x), ...
    'Units', 'normalized', ...
    'Interpreter', 'latex', ...
    'FontSize', 18, ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top');

% Leyenda
legend([h1,h2], {'Extradós','Intradós'}, ...
    'Location', 'northwest', ...
    'Interpreter', 'none', ...
    'FontSize', 13);

% Limites
xlim([0 1]);

% Ajuste automatico razonable en y
ymax = max([-Cp_upper; -Cp_lower]);
ymin = min([-Cp_upper; -Cp_lower]);
ylim([floor(ymin)-0.2, ceil(ymax)+0.2]);

% Estilo de ejes
set(gca, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 15, ...
    'LineWidth', 1);

hold off;


%% Passo 5: Calculo do coeficiente de sustentacao CL

c = 1;

% Coeficientes de força nos eixos x e z
Cx = -sum(Cp .* n_x .* DeltaS) / c;

Cz = -sum(Cp .* n_z .* DeltaS) / c;

% Sustentacion perpendicular al escoamento livre
CL = -Cx*sin(AoA_rad) + Cz*cos(AoA_rad);

% Coeficiente de resistencia de presion
CD = Cx*cos(AoA_rad) + Cz*sin(AoA_rad);

% Mostrar resultados
fprintf('\n--- Coeficientes Aerodinâmicos ---\n');
fprintf('Ángulo de ataque: %.2f deg\n', rad2deg(AoA_rad));
fprintf('Cx = %.6f\n', Cx);
fprintf('Cz = %.6f\n', Cz);
fprintf('CL = %.6f\n', CL);
fprintf('CD (pressão) = %.6f\n', CD);
%% Mapa de contornos de Cp alrededor del perfil Roncz

% Flujo libre
[u_inf, w_inf] = EsUniforme(V_inf_modul, AoA_rad);
V_inf = hypot(u_inf, w_inf);

% ---------- 1. Crear malla alrededor del perfil ----------
x_min_dom = -0.2;
x_max_dom = 1.2;
z_min_dom = -0.35;
z_max_dom = 0.35;

nxg = 250;   % puntos en x
nzg = 180;   % puntos en z

xg = linspace(x_min_dom, x_max_dom, nxg);
zg = linspace(z_min_dom, z_max_dom, nzg);

[X, Z] = meshgrid(xg, zg);

% Inicializar campos
U = zeros(size(X));
W = zeros(size(X));
Cp_field = zeros(size(X));

% ---------- 2. Poligono del perfil para enmascarar interior ----------
x_profile = [x_start; x_end(end)];
z_profile = [z_start; z_end(end)];

% Si el contorno no está cerrado, cerrarlo:
if x_profile(1) ~= x_profile(end) || z_profile(1) ~= z_profile(end)
    x_profile = [x_profile; x_profile(1)];
    z_profile = [z_profile; z_profile(1)];
end

% ---------- 3. Calcular velocidad en cada punto de la malla ----------
for ii = 1:size(X,1)
    for jj = 1:size(X,2)

        xp = X(ii,jj);
        zp = Z(ii,jj);

        % Empezar con velocidad de flujo libre
        u_total = u_inf;
        w_total = w_inf;

        % Sumar contribucion de todos los paneles
        for k = 1:N

            [u_start, w_start, ~, ~] = VorticeDist(1,0, ...
                xp, zp, ...
                x_start(k), z_start(k), ...
                x_end(k), z_end(k));

            [~, ~, u_end, w_end] = VorticeDist(0,1, ...
                xp, zp, ...
                x_start(k), z_start(k), ...
                x_end(k), z_end(k));

            u_total = u_total + gamma(k)*u_start + gamma(k+1)*u_end;
            w_total = w_total + gamma(k)*w_start + gamma(k+1)*w_end;
        end

        U(ii,jj) = u_total;
        W(ii,jj) = w_total;

        Vmag = hypot(u_total, w_total);
        Cp_field(ii,jj) = 1 - (Vmag / V_inf)^2;
    end
end

% ---------- 4. Enmascarar el interior del perfil ----------
in = inpolygon(X, Z, x_profile, z_profile);
Cp_field(in) = NaN;

% ---------- 5. Graficar contornos tipo CFD ----------
figure('Color','w');
hold on;

% Numero de niveles de contorno
niveles = linspace(min(Cp_field(:),[],'omitnan'), ...
                   max(Cp_field(:),[],'omitnan'), 60);

contourf(X, Z, Cp_field, niveles, 'LineColor', 'none');

% Perfil en blanco con borde negro
patch(x_profile, z_profile, [1 1 1], ...
    'EdgeColor', 'k', 'LineWidth', 1.5);

axis equal;
xlim([x_min_dom x_max_dom]);
ylim([z_min_dom z_max_dom]);

xlabel('$x/c$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$z/c$', 'Interpreter', 'latex', 'FontSize', 14);
title(sprintf('Distribuicao de Cp ao redor do perfil Roncz (AoA = %.1f deg)', x), ...
    'Interpreter', 'none');

colormap(jet);
cb = colorbar;
cb.Label.String = 'Color Map: C_p';
cb.Label.FontSize = 12;

box on;
grid off;
hold off;

%% Funcao local: sem graficos para testes repetidos
function r = resolver_paineis_C1(filename,N_deseado,AoA_graus,V_inf_modul)
[xs,zs,xe,ze,xc,zc,nx,nz,tx,tz,dS] = panelizationRoncz(filename,N_deseado,false);
N=numel(dS);
A=zeros(N+1,N+1);
% Guardar influencias unitarias para evitar recalcular VorticeDist no Cp
Us=zeros(N,N); Ws=zeros(N,N); Ue=zeros(N,N); We=zeros(N,N);
for i=1:N
    for j=1:N
        [us,ws,~,~]=VorticeDist(1,0,xc(i),zc(i),xs(j),zs(j),xe(j),ze(j));
        [~,~,ue,we]=VorticeDist(0,1,xc(i),zc(i),xs(j),zs(j),xe(j),ze(j));
        Us(i,j)=us; Ws(i,j)=ws; Ue(i,j)=ue; We(i,j)=we;
        A(i,j)=A(i,j)+us*nx(i)+ws*nz(i);
        A(i,j+1)=A(i,j+1)+ue*nx(i)+we*nz(i);
    end
end
A(N+1,1)=1; A(N+1,N+1)=1;
alpha=deg2rad(AoA_graus);
[u_inf,w_inf]=EsUniforme(V_inf_modul,alpha);
b=[-(u_inf*nx+w_inf*nz);0];
gamma=A\b;
Ut=u_inf+Us*gamma(1:N)+Ue*gamma(2:N+1);
Wt=w_inf+Ws*gamma(1:N)+We*gamma(2:N+1);
Vt=Ut.*tx+Wt.*tz;
Vn=Ut.*nx+Wt.*nz;
Cp=1-(Vt/V_inf_modul).^2;
Cx=-sum(Cp.*nx.*dS); Cz=-sum(Cp.*nz.*dS);
r.CL=-Cx*sin(alpha)+Cz*cos(alpha);
r.CD_pressure=Cx*cos(alpha)+Cz*sin(alpha);
r.erro_normal=max(abs(Vn))/V_inf_modul;
r.residuo=max(abs(A*gamma-b));
r.gamma=gamma; r.Cp=Cp; r.Vt=Vt; r.Cp_min=min(Cp);
[~,im]=min(Cp); r.x_Cp_min=xc(im);
%---------------------

%% Localizacao do ponto de estagnacao

% Identificar mudancas de sinal da velocidade tangencial
idx = find(Vt(1:end-1).*Vt(2:end) < 0);

% Inicializar coordenadas
r.x_est = NaN;
r.z_est = NaN;

% Procurar a mudanca de sinal mais proxima do bordo de ataque
if ~isempty(idx)

    % Coordenadas interpoladas para cada mudanca de sinal
    x_candidatos = zeros(length(idx),1);
    z_candidatos = zeros(length(idx),1);

    for k = 1:length(idx)

        i = idx(k);

        % Fator de interpolacao linear
        lambda = -Vt(i)/(Vt(i+1)-Vt(i));

        % Coordenadas do ponto onde Vt = 0
        x_candidatos(k) = xc(i) + ...
            lambda*(xc(i+1)-xc(i));

        z_candidatos(k) = zc(i) + ...
            lambda*(zc(i+1)-zc(i));

    end

    % Selecionar candidato mais proximo do bordo de ataque
    [~,imin] = min(x_candidatos);

    r.x_est = x_candidatos(imin);
    r.z_est = z_candidatos(imin);
    fprintf('AoA = %.1f deg: %d mudancas de sinal de Vt\n', ...
        AoA_graus, numel(idx));
end
%------------------
r.x_start=xs; r.z_start=zs; r.x_end=xe; r.z_end=ze;
end
