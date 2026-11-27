clear; close all; clc;

%% 0. Configuracao

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
figDir = fullfile(scriptDir, 'figures');

if ~isfolder(figDir)
    mkdir(figDir);
end

% Cores dos metodos
colorThin  = [0 0.4470 0.7410];
colorPanel = [0.8500 0.3250 0.0980];
colorCFD   = [0.4660 0.6740 0.1880];

%% 0.1. Atualizar automaticamente os resultados dos metodos

thinDir = fullfile(rootDir, 'thin-airfoil');
panelDir = fullfile(rootDir, 'panel-method');

runThinAirfoil(thinDir);
addpath(panelDir);
exportResults("all");

% Os resultados CFD sao importados das simulacoes externas.

%% 1. Perfis delgados

T = readtable(fullfile(rootDir, ...
    'thin-airfoil', 'results.csv'), ...
    'TextType', 'string');

idx = T.kind == "series" & T.key == "cl_curve";

alphaThin = T.x(idx);
clThin = T.value(idx);

timeThin = T.value(T.kind == "scalar" & ...
    T.key == "execution_time_s");

%% 2. Metodo dos paineis

T = readtable(fullfile(rootDir, ...
    'panel-method', 'results.csv'));

alphaPanel = T.alpha_deg;
clPanel = T.cl;

timePanel = sum(T.time_s);

%% 3. CFD

CFD = readtable(fullfile(rootDir, 'CFD', 'results.csv'));

alphaCFD = CFD.aoa;
clCFD = CFD.cl;

%% 4. Comparacao da sustentacao

fig = figure('Name', 'Comparacao de metodos', ...
    'Color', 'w', 'WindowStyle', 'docked');

ax = axes(fig);
hold(ax, 'on');

plot(ax, alphaThin, clThin, '-o', ...
    'Color', colorThin, ...
    'LineWidth', 1.5, ...
    'MarkerSize', 5, ...
    'MarkerFaceColor', 'w', ...
    'MarkerIndices', find(ismember(alphaThin, -20:5:20)), ...
    'DisplayName', 'Perfis delgados');

plot(ax, alphaPanel, clPanel, '-s', ...
    'Color', colorPanel, ...
    'LineWidth', 1.5, ...
    'MarkerSize', 5, ...
    'MarkerFaceColor', 'w', ...
    'DisplayName', 'Método dos painéis');

plot(ax, alphaCFD, clCFD, '-^', ...
    'Color', colorCFD, ...
    'LineWidth', 1.5, ...
    'MarkerSize', 5, ...
    'MarkerFaceColor', 'w', ...
    'DisplayName', 'CFD');

xline(ax, 0, 'k:', 'HandleVisibility', 'off');
yline(ax, 0, 'k:', 'HandleVisibility', 'off');

xlabel(ax, '$\alpha$ [$^\circ$]', ...
    'Interpreter', 'latex');

ylabel(ax, '$c_l$', 'Interpreter', 'latex');

xlim(ax, [-20 20]);
xticks(ax, -20:5:20);

grid(ax, 'on');
box(ax, 'on');

ax.FontSize = 12;
ax.GridAlpha = 0.12;

legend(ax, 'Location', 'northwest', 'Box', 'off');

hold(ax, 'off');

%% 5. Exportacao

exportgraphics(ax, ...
    fullfile(figDir, 'lift_comparison.pdf'), ...
    'ContentType', 'vector');



%% 6. Analise quantitativa

methods = ["Perfis delgados"; "Paineis"; "CFD"];

alphaData = {alphaThin; alphaPanel; alphaCFD};
clData = {clThin; clPanel; clCFD};

% Angulos comuns aos tres metodos
alphaCommon = (0:5:20)';

% Intervalo utilizado no ajuste linear
fitRange = [0 10];

% Tolerancia relativa para identificar a zona linear
tolSlope = 0.10;

nMethods = numel(methods);

slope_deg = zeros(nMethods, 1);
slope_rad = zeros(nMethods, 1);
alpha_L0_deg = zeros(nMethods, 1);
R2 = zeros(nMethods, 1);
linear_max_deg = zeros(nMethods, 1);

for k = 1:nMethods

    % Selecionar os angulos comuns
    [found, idx] = ismember(alphaCommon, alphaData{k});
    assert(all(found), 'Faltam angulos de ataque.');

    a = alphaCommon;
    cl = clData{k}(idx);
    cl = cl(:);

    % Ajuste linear
    idxFit = a >= fitRange(1) & a <= fitRange(2);

    p = polyfit(a(idxFit), cl(idxFit), 1);
    clFit = polyval(p, a(idxFit));

    slope_deg(k) = p(1);
    slope_rad(k) = p(1)*180/pi;
    alpha_L0_deg(k) = -p(2)/p(1);

    R2(k) = 1 - sum((cl(idxFit)-clFit).^2) ...
        / sum((cl(idxFit)-mean(cl(idxFit))).^2);

    % Identificacao da regiao aproximadamente linear
    localSlope = diff(cl)./diff(a);
    refSlope = localSlope(1);

    relativeError = abs(localSlope-refSlope) ...
        / abs(refSlope);

    firstBad = find(relativeError > tolSlope, 1);

    if isempty(firstBad)
        linear_max_deg(k) = a(end);
    else
        linear_max_deg(k) = a(firstBad);
    end

end

%% 7. Comparacao dos tempos de calculo

time_s = [timeThin; timePanel; NaN];

n_angles = [
    numel(alphaThin);
    numel(alphaPanel);
    numel(alphaCFD)
];

%% 8. Exportacao dos resultados quantitativos

results = table(methods, ...
    slope_deg, slope_rad, alpha_L0_deg, ...
    R2, linear_max_deg, time_s, n_angles);

disp(results);

writetable(results, ...
    fullfile(scriptDir, 'results.csv'));

% Atualizar os valores usados na tabela de tempos do relatorio.
% As macros guardam os tempos medidos, sem arredondamento de apresentacao.
assert(all(isfinite(time_s(1:2)) & time_s(1:2) > 0), ...
    'Tempos de calculo invalidos.');
configDir = fullfile(rootDir, 'relatorio', 'config');
timingFile = fullfile(configDir, 'tempos_calculo.tex');
timingText = sprintf([ ...
    '%% Gerado automaticamente por Comparison/main.m a partir dos CSV.\n' ...
    '%% Unidades: perfis delgados em ms; paineis em s.\n' ...
    '\\newcommand{\\TempoPerfisDelgadosMs}{%.15g}\n' ...
    '\\newcommand{\\TempoPaineisS}{%.15g}\n' ...
    '\\newcommand{\\NumAngulosPerfisDelgados}{%d}\n' ...
    '\\newcommand{\\NumAngulosPaineis}{%d}\n'], ...
    1000*time_s(1), time_s(2), n_angles(1), n_angles(2));

% Substituicao atomica para evitar leituras incompletas pelo latexmk.
tmpTimingFile = [tempname(configDir) '.tex'];
writelines(timingText, tmpTimingFile);
[ok, msg] = movefile(tmpTimingFile, timingFile, 'f');
assert(ok, 'Erro ao exportar os tempos para LaTeX: %s', msg);

% Atualizar os restantes valores numericos do relatorio.
exportReportResults();



%% 9. Estudo de convergencia dos paineis

conv = readtable(fullfile(panelDir, 'convergence.csv'));

N_final = 260;
tolPanel = 0.5;

figConv = figure('Name', 'Convergencia dos paineis', ...
    'Color', 'w', 'WindowStyle', 'docked');

layout = tiledlayout(figConv, 2, 1, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

% Coeficiente de sustentacao

axCl = nexttile(layout);

plot(axCl, conv.N, conv.cl, '-o', ...
    'Color', colorPanel, ...
    'LineWidth', 1.5, ...
    'MarkerSize', 5, ...
    'MarkerFaceColor', 'w');

xline(axCl, N_final, 'k--', ...
    sprintf('N = %d', N_final));

yline(axCl, conv.cl(end), 'k:', ...
    'Referencia');

xlabel(axCl, 'Número de painéis, N');
ylabel(axCl, '$c_l$', 'Interpreter', 'latex');

grid(axCl, 'on');
box(axCl, 'on');

axCl.FontSize = 12;
axCl.GridAlpha = 0.12;

% Erro relativo

axError = nexttile(layout);

plot(axError, conv.N, conv.relative_error_pct, '-o', ...
    'Color', colorPanel, ...
    'LineWidth', 1.5, ...
    'MarkerSize', 5, ...
    'MarkerFaceColor', 'w');

yline(axError, tolPanel, 'k--', ...
    sprintf('Tolerância = %.1f%%', tolPanel));

xline(axError, N_final, 'k--', ...
    sprintf('N = %d', N_final));

xlabel(axError, 'Número de painéis, N');
ylabel(axError, 'Erro relativo [%]');

grid(axError, 'on');
box(axError, 'on');

axError.FontSize = 12;
axError.GridAlpha = 0.12;

linkaxes([axCl, axError], 'x');

xlim(axCl, [0 820]);
xticks(axCl, 0:100:800);
xticks(axError, 0:100:800);

% Exportacao individual

exportgraphics(axCl, ...
    fullfile(figDir, 'panel_convergence_cl.pdf'), ...
    'ContentType', 'vector');

exportgraphics(axError, ...
    fullfile(figDir, 'panel_convergence_error.pdf'), ...
    'ContentType', 'vector');


%% 10. Painelizacao

filename = fullfile(rootDir, 'airfoil-points.txt');

[xs, zs, xe, ze, xc, zc] = ...
    panelizationRoncz(filename, N_final, false);

% Tamanho fixo (largo e baixo): com axis equal, uma janela alta deixa
% muito espaco em branco por cima do perfil
figPanel = figure('Name', 'Painelização', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2 2 18 6]);
if exist('theme', 'file')   % R2025a+: evita o tema escuro
    theme(figPanel, 'light');
end

axPanel = axes(figPanel);

hold(axPanel, 'on');

% Paineis e nos

plot(axPanel, [xs; xe(end)], [zs; ze(end)], ...
    'k-o', ...
    'LineWidth', 1, ...
    'MarkerSize', 2.5, ...
    'MarkerFaceColor', 'w', ...
    'DisplayName', 'Painéis e nós');

% Pontos de controlo

plot(axPanel, xc, zc, '.', ...
    'Color', colorPanel, ...
    'MarkerSize', 6, ...
    'DisplayName', 'Pontos de controlo');

axis(axPanel, 'equal');

xlim(axPanel, [-0.05 1.05]);
ylim(axPanel, [min([zs; ze])-0.03, ...
               max([zs; ze])+0.03]);

xlabel(axPanel, '$x/c$', 'Interpreter', 'latex');
ylabel(axPanel, '$z/c$', 'Interpreter', 'latex');

grid(axPanel, 'on');
box(axPanel, 'on');

axPanel.FontSize = 12;
axPanel.GridAlpha = 0.12;

legend(axPanel, 'Location', 'northoutside', ...
    'Orientation', 'horizontal', 'Box', 'off');

hold(axPanel, 'off');

exportgraphics(axPanel, ...
    fullfile(figDir, 'panelization.pdf'), ...
    'ContentType', 'vector');



%% 11. Linhas de corrente - Metodo dos paineis

data = load(fullfile(panelDir, 'cases.mat'));
cases = data.cases;

for k = 1:numel(cases)

    r = cases(k);

    figStream = figure('Name', sprintf( ...
        'Linhas de corrente - %g graus', r.alpha_deg), ...
        'Color', 'w', 'WindowStyle', 'docked');

    layoutStream = tiledlayout(figStream, 1, 1, ...
        'Padding', 'loose', 'TileSpacing', 'loose');

    axStream = nexttile(layoutStream);

    % Referencial do escoamento livre: rotacao de -alpha em torno
    % do meio da corda, apenas para representacao grafica.
    alpha = deg2rad(r.alpha_deg);
    R = [cos(alpha), sin(alpha); -sin(alpha), cos(alpha)];
    pivot = [0.5; 0];
    xView = [-0.4 1.4];
    zView = [-0.5 0.5];

    % Integrar num dominio suficientemente amplo para cobrir a janela
    % depois da rotacao, evitando cortes nas linhas de corrente.
    corners = [xView(1) xView(1) xView(2) xView(2); ...
               zView(1) zView(2) zView(1) zView(2)];
    originalCorners = pivot + R.'*(corners-pivot);
    margin = 0.12;
    bounds = [min(originalCorners(1,:))-margin, ...
              max(originalCorners(1,:))+margin, ...
              min(originalCorners(2,:))-margin, ...
              max(originalCorners(2,:))+margin];

    desenhar_linhas_corrente_sustentador( ...
        r.gamma, ...
        r.x_start, r.z_start, ...
        r.x_end, r.z_end, ...
        r.V_inf, alpha, axStream, bounds);

    % Rodar linhas, setas e contorno sem alterar os resultados do solver.
    % findall inclui linhas e setas com HandleVisibility='off'.
    hLines = findall(axStream, 'Type', 'line');
    assert(~isempty(hLines), 'Nenhuma linha de corrente encontrada.');
    for h = hLines(:).'
        x = h.XData; z = h.YData;
        q = pivot + R*([x(:).'; z(:).']-pivot);
        h.XData = reshape(q(1,:),size(x));
        h.YData = reshape(q(2,:),size(z));
    end
    for h = findall(axStream, 'Type', 'patch').'
        x = h.XData; z = h.YData;
        q = pivot + R*([x(:).'; z(:).']-pivot);
        h.XData = reshape(q(1,:),size(x));
        h.YData = reshape(q(2,:),size(z));
    end

    hold(axStream, 'on');

    % Ponto de estagnacao
    if isfinite(r.x_est) && isfinite(r.z_est)

        qStag = pivot + R*([r.x_est; r.z_est]-pivot);
        hStag = plot(axStream, qStag(1), qStag(2), 'o', ...
            'Color', colorPanel, ...
            'MarkerFaceColor', colorPanel, ...
            'MarkerSize', 6, ...
            'DisplayName', 'Ponto de estagnação');

        lgd = legend(axStream, hStag, 'Box', 'off');
        lgd.Layout.Tile = 'south';

    end

    % Configuracao dos eixos

    axis(axStream, 'equal');

    xlim(axStream, xView);
    ylim(axStream, zView);

    xlabel(axStream, '$x''/c$', 'Interpreter', 'latex');
    ylabel(axStream, '$z''/c$', 'Interpreter', 'latex');

    grid(axStream, 'on');
    box(axStream, 'on');

    axStream.FontSize = 12;
    axStream.GridAlpha = 0.12;

    hold(axStream, 'off');

    % Exportacao

    drawnow;
    exportgraphics(layoutStream, ...
        fullfile(figDir, sprintf( ...
        'panel_streamlines_%02d.pdf', r.alpha_deg)), ...
        'ContentType', 'vector');

end

%% 12. Resultados dos casos especificos

panelCases = readtable(fullfile(panelDir, 'cases.csv'));

disp(panelCases);

%% 13. Distribuicao do coeficiente de pressao - Metodo dos paineis

cpPanel = readtable(fullfile(panelDir, 'cp.csv'), 'TextType', 'string');

for a = [0 5 20]
    T = cpPanel(cpPanel.alpha_deg == a,:);
    upper = sortrows(T(T.surface == "upper",:), 'x_c');
    lower = sortrows(T(T.surface == "lower",:), 'x_c');
    assert(~isempty(upper) && ~isempty(lower), ...
        'Distribuicao de Cp em falta para alpha = %g graus.', a);

    figCp = figure('Name', sprintf('Cp - Paineis - %g graus', a), ...
        'Color', 'w', 'WindowStyle', 'docked');
    axCp = axes(figCp);
    hold(axCp, 'on');

    plot(axCp, upper.x_c, upper.cp, '-', ...
        'Color', colorThin, 'LineWidth', 1.5, ...
        'DisplayName', 'Extradorso');
    plot(axCp, lower.x_c, lower.cp, '-', ...
        'Color', colorPanel, 'LineWidth', 1.5, ...
        'DisplayName', 'Intradorso');

    % Convencao aerodinamica: Cp negativo no topo, sem alterar os valores.
    axCp.YDir = 'reverse';
    xlim(axCp, [0 1]);
    xlabel(axCp, '$x/c$', 'Interpreter', 'latex');
    ylabel(axCp, '$C_p$', 'Interpreter', 'latex');
    title(axCp, sprintf('$\\alpha = %g^{\\circ}$', a), ...
        'Interpreter', 'latex');
    grid(axCp, 'on');
    box(axCp, 'on');
    axCp.FontSize = 12;
    axCp.GridAlpha = 0.12;
    legend(axCp, 'Location', 'best', 'Box', 'off');
    hold(axCp, 'off');

    exportgraphics(axCp, ...
        fullfile(figDir, sprintf('panel_cp_%02d.pdf', a)), ...
        'ContentType', 'vector');
end

%% Funcoes auxiliares

function runThinAirfoil(thinDir)
    % O clear do script original fica isolado deste script.
    run(fullfile(thinDir, 'main.m'));
end
