function fig = desenhar_linhas_corrente_sustentador(gamma, x_start, z_start, x_end, z_end, V_inf_modul, AoA_rad, ax, bounds)
% DESENHAR_LINHAS_CORRENTE_SUSTENTADOR Integra o campo de velocidades e representa
% graficamente as linhas de corrente em redor de um perfil (problema sustentador).
%
% SINTAXE:
%   fig = desenhar_linhas_corrente_sustentador(gamma, x_start, z_start, x_end, z_end)
%   fig = desenhar_linhas_corrente_sustentador(gamma, x_start, z_start, x_end, z_end, V_inf_modul, AoA_rad)
%
% PARÂMETROS DE ENTRADA:
%   gamma       - Vetor de circulação (intensidade dos vórtices) nos nós [N+1 x 1]
%   x_start     - Coordenadas x do ponto inicial de cada painel [N x 1]
%   z_start     - Coordenadas z do ponto inicial de cada painel [N x 1]
%   x_end       - Coordenadas x do ponto final de cada painel [N x 1]
%   z_end       - Coordenadas z do ponto final de cada painel [N x 1]
%   V_inf_modul - (Opcional) Módulo da velocidade do escoamento incidente (por defeito: 1 m/s)
%   AoA_rad     - (Opcional) Ângulo de ataque em radianos (por defeito: 0 rad)
%   ax          - (Opcional) Eixos existentes para desenhar as linhas de corrente
%   bounds      - (Opcional) Dominio de integracao [xmin xmax zmin zmax]
%
% PARÂMETROS DE SAÍDA:
%   fig         - Identificador (handle) da figura gerada

    % ---------------------------------------------------------------------
    % VERIFICAÇÃO E DEFINIÇÃO DE PARÂMETROS OPCIONAIS POR DEFEITO
    % ---------------------------------------------------------------------
    if nargin < 6 || isempty(V_inf_modul)
        V_inf_modul = 1; % Velocidade de referência por defeito (1 m/s)
    end
    if nargin < 7 || isempty(AoA_rad)
        AoA_rad = 0;     % Ângulo de ataque por defeito (0 rad)
    end

    % Avaliar a velocidade diretamente nos pontos da integracao, sem
    % interpolar uma grelha com NaN junto ao contorno do perfil.
    x_start = x_start(:); z_start = z_start(:);
    x_end = x_end(:); z_end = z_end(:);
    gamma = gamma(:);
    x_contorno = [x_start; x_end(end)];
    z_contorno = [z_start; z_end(end)];
    dx = x_end-x_start; dz = z_end-z_start;
    L = hypot(dx,dz);
    C = -dx./L; S = dz./L;
    C(abs(C)<1e-14) = 0;
    S(abs(S)<1e-14) = 0;
    G = gamma(1:end-1);
    dG = -diff(gamma)./L;
    [u_inf,w_inf] = EsUniforme(V_inf_modul,AoA_rad);
    velocity = @(q) panelVelocity(q,x_start,z_start,L,C,S,G,dG,u_inf,w_inf);

    % Confirmar a equivalencia com a implementacao original dos vortices.
    for q = [-0.3 0.3 1.2; 0.1 0.2 -0.15]
        v = [u_inf; w_inf];
        for j = 1:numel(L)
            [us,ws,ue,we] = VorticeDist(gamma(j),gamma(j+1), ...
                q(1),q(2),x_start(j),z_start(j),x_end(j),z_end(j));
            v = v + [us+ue;ws+we];
        end
        assert(norm(velocity(q)-v) < 1e-9*max(1,norm(v)), ...
            'Campo de velocidades inconsistente com VorticeDist.');
    end

    distance = @(q) surfaceDistance(q,x_start,z_start,dx,dz,L, ...
        x_contorno,z_contorno);
    if nargin < 9 || isempty(bounds)
        bounds = [-0.5 1.5 -0.6 0.6];
    end
    opts = odeset('RelTol',1e-7,'AbsTol',1e-9,'MaxStep',0.01, ...
        'Events',@(t,q) stopTrace(t,q,velocity,distance,bounds,V_inf_modul));

    disp('A integrar as linhas de corrente diretamente no campo dos paineis...');

    % =====================================================================
    % PASSO 4: GERAÇÃO E CONFIGURAÇÃO DA FIGURA GRÁFICA
    % =====================================================================
    if nargin < 8 || isempty(ax)
        fig = figure('Name', 'Linhas de Corrente (Sustentador)', ...
            'Position', [200, 200, 800, 450]);
        ax = axes(fig);
    else
        fig = ancestor(ax, 'figure');
    end
    hold(ax, 'on');

    % Semear em torno do perfil e integrar em ambos os sentidos.
    % Assim ha linhas junto a ambas as superficies mesmo a alpha elevado.
    spacing = 0.055;
    xMid = 0.4;
    cross = (x_start <= xMid & x_end > xMid) | ...
            (x_end <= xMid & x_start > xMid);
    zz = z_start(cross) + (xMid-x_start(cross))./dx(cross).*dz(cross);
    zUpper = max(zz);
    zLower = min(zz);
    zSeeds = [zLower-0.018:-spacing:bounds(3)+0.002, ...
              zUpper+0.018:spacing:bounds(4)-0.002];

    paths = cell(0,1);
    for z = zSeeds
        seed = [xMid; z];
        upstream = traceFlow(seed,-1,velocity,opts);
        downstream = traceFlow(seed,1,velocity,opts);
        paths{end+1,1} = [flipud(upstream); downstream(2:end,:)];
    end

    % Completar as linhas que entram e saem sem cruzar a secao central.
    zIn = bounds(3)+spacing/2:spacing:bounds(4)-0.002;
    xIn = repmat(bounds(1)+0.002,size(zIn));
    if abs(w_inf) > 0.02*V_inf_modul
        stepX = spacing*V_inf_modul/abs(w_inf);
        xExtra = bounds(1)+stepX/2:stepX:bounds(2)-0.002;
        zExtra = repmat(bounds(3)+0.002,size(xExtra));
        if w_inf < 0
            zExtra(:) = bounds(4)-0.002;
        end
        xIn = [xIn xExtra];
        zIn = [zIn zExtra];
    end
    for k = 1:numel(xIn)
        seed = [xIn(k);zIn(k)];
        xy = traceFlow(seed,1,velocity,opts);
        if ~isempty(xy) && (max(xy(:,1)) < xMid || min(xy(:,1)) > xMid)
            paths{end+1,1} = xy;
        end
    end

    % Agrupar as setas num unico objeto grafico para exportacao vetorial.
    arrowX = [];
    arrowZ = [];

    for k = 1:numel(paths)
        xy = paths{k};

        % stream2 pode devolver NaN ao chegar a uma celula mascarada.
        % Terminar no primeiro ponto invalido, sem ligar trocos separados.
        firstBad = find(~all(isfinite(xy), 2), 1);
        if ~isempty(firstBad)
            xy = xy(1:firstBad-1,:);
        end
        if size(xy,1) < 2
            continue
        end

        plot(ax, xy(:,1), xy(:,2), 'k-', 'LineWidth', 0.8, ...
            'HandleVisibility', 'off');

        % Setas tangentes, sem apagar segmentos da linha.
        arc = [0; cumsum(hypot(diff(xy(:,1)), diff(xy(:,2))))];
        [arc, idx] = unique(arc, 'stable');
        xy = xy(idx,:);
        if numel(arc) < 2 || any(~isfinite(arc))
            continue
        end

        for a = 0.18:0.35:arc(end)-0.04
            q = interp1(arc, xy, [a-0.015; a+0.015]);
            d = q(2,:) - q(1,:);
            len = hypot(d(1),d(2));
            if isfinite(len) && len > 0 && distance(q(2,:)') > 0.03
                d = d/len;
                n = [-d(2), d(1)];
                tip = q(2,:);
                left = tip - 0.025*d + 0.012*n;
                right = tip - 0.025*d - 0.012*n;
                arrowX = [arrowX, left(1), tip(1), right(1), NaN];
                arrowZ = [arrowZ, left(2), tip(2), right(2), NaN];
            end
        end
    end

    plot(ax, arrowX, arrowZ, 'k-', 'LineWidth', 0.8, ...
        'HandleVisibility', 'off');

    fill(ax, x_contorno, z_contorno, [0.8 0.8 0.8], ...
        'EdgeColor', 'k', 'LineWidth', 1.5);

    axis(ax, 'equal');
    xlim(ax, [-0.4, 1.4]);
    ylim(ax, [-0.5, 0.5]);
    grid(ax, 'on');
    box(ax, 'on');

    xlabel(ax, 'x');
    ylabel(ax, 'z');

    hold(ax, 'off');
    disp('Linhas de corrente geradas com sucesso.');
end

function v = panelVelocity(q,xs,zs,L,C,S,G,dG,uInf,wInf)
    % Mesma expressao de VorticeDist, somando as duas contribuicoes
    % lineares de cada painel e avaliando todos os paineis em simultaneo.
    x = C.*(q(1)-xs)-S.*(q(2)-zs);
    z = S.*(q(1)-xs)+C.*(q(2)-zs)+1e-5;
    beta = atan2(z,x+L)-atan2(z,x);
    ell = log((x.^2+z.^2)./((x+L).^2+z.^2));
    up = G.*beta/(2*pi)-dG.*(z.*ell-2*x.*beta)/(4*pi);
    wp = -G.*ell/(4*pi)-dG.*(x.*ell/2+L+z.*beta)/(2*pi);
    v = [uInf-sum(C.*up+S.*wp); wInf+sum(S.*up-C.*wp)];
end

function d = surfaceDistance(q,xs,zs,dx,dz,L,xc,zc)
    lambda = ((q(1)-xs).*dx+(q(2)-zs).*dz)./L.^2;
    lambda = max(0,min(1,lambda));
    d = min(hypot(q(1)-xs-lambda.*dx,q(2)-zs-lambda.*dz));
    if inpolygon(q(1),q(2),xc,zc)
        d = -d;
    end
end

function [value,isterminal,direction] = stopTrace(~,q,velocity,distance,bounds,Vinf)
    value = [q(1)-bounds(1); bounds(2)-q(1); ...
             q(2)-bounds(3); bounds(4)-q(2); ...
             distance(q)-2e-4; norm(velocity(q))-1e-5*Vinf];
    isterminal = ones(6,1);
    direction = -ones(6,1);
end

function xy = traceFlow(seed,signFlow,velocity,opts)
    [~,xy] = ode45(@(~,q) flowDirection(q,signFlow,velocity),[0 5],seed,opts);
end

function dq = flowDirection(q,signFlow,velocity)
    v = velocity(q);
    dq = signFlow*v/max(norm(v),eps);
end
