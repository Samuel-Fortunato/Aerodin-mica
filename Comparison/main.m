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
