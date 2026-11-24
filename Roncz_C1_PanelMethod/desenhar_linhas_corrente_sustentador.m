function fig = desenhar_linhas_corrente_sustentador(gamma, x_start, z_start, x_end, z_end, V_inf_modul, AoA_rad)
% DESENHAR_LINHAS_CORRENTE_SUSTENTADOR Calcula o campo de velocidades e representa
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

    % Número total de painéis
    N = length(x_start);

    disp('A calcular o campo de velocidades para as linhas de corrente (sustentador)...');

    % =====================================================================
    % PASSO 1: CRIAÇÃO DA MALHA ESPACIAL (GRELHA) EM REDOR DO PERFIL
    % =====================================================================
    x_grid = linspace(-0.5, 1.5, 80); 
    z_grid = linspace(-0.6, 0.6, 60); 
    [X, Z] = meshgrid(x_grid, z_grid);

    % Pré-alocação das matrizes de componentes de velocidade (U, W)
    U = zeros(length(x_grid), length(z_grid));
    W = zeros(length(x_grid), length(z_grid));

    % =====================================================================
    % PASSO 2: CÁLCULO DO CAMPO DE VELOCIDADES TOTAIS
    % =====================================================================
    u = 0; 
    w = 0;

    for i = 1:length(x_grid)
        for k = 1:length(z_grid)
            
            % Contribuição induzida por cada painel de vórtices
            for j = 1:N
                [u_s, w_s, u_e, w_e] = VorticeDist(gamma(j), gamma(j+1), x_grid(i), z_grid(k), ...
                    x_start(j), z_start(j), x_end(j), z_end(j));
                u = u + u_s + u_e;
                w = w + w_s + w_e;
            end
            
            % Contribuição do escoamento uniforme incidente
            [u_uniforme, w_uniforme] = EsUniforme(V_inf_modul, AoA_rad);
            
            % Armazenar a velocidade total no ponto atual da malha
            U(i, k) = u + u_uniforme;
            W(i, k) = w + w_uniforme;
            
            % Reiniciar acumuladores de velocidade
            u = 0;
            w = 0;
        end  
    end

    % =====================================================================
    % PASSO 3: MASCARAMENTO DO INTERIOR DO PERFIL
    % =====================================================================
    x_contorno = [x_start(:); x_end(end)];
    z_contorno = [z_start(:); z_end(end)];

    in = inpolygon(X', Z', x_contorno, z_contorno);

    U(in) = NaN; 
    W(in) = NaN;

    % =====================================================================
    % PASSO 4: GERAÇÃO E CONFIGURAÇÃO DA FIGURA GRÁFICA
    % =====================================================================
    fig = figure('Name', 'Linhas de Corrente (Sustentador)', 'Position', [200, 200, 800, 450]);
    hold on;

    h = streamslice(x_grid, z_grid, U', W'); 
    set(h, 'Color', 'k');

    fill(x_contorno, z_contorno, [0.8 0.8 0.8], 'EdgeColor', 'k', 'LineWidth', 1.5);

    axis equal;
    xlim([-0.4, 1.4]);
    ylim([-0.5, 0.5]);
    grid on; 
    box on;

    % title(sprintf('Linhas de Corrente - Sustentador (\\alpha = %.1f^\\circ)', AoA_rad * 180 / pi));
    xlabel('x');
    ylabel('z');

    hold off;
    disp('Linhas de corrente geradas com sucesso.');
end