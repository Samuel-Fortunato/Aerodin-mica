clear;
close all;
clc;

%% 0. Configuracao

scriptDir = fileparts(mfilename('fullpath'));
filename = fullfile(scriptDir, '..', 'airfoil-points.txt');
outputFile = fullfile(scriptDir, 'results.csv');
figDir = fullfile(scriptDir, 'figures');
showLiftCurve = false; % Grafico individual opcional, para diagnostico.

N = 1000;
N_test = [50 100 200 500 1000 2000 5000 10000];
alpha_deg = -20:0.5:20;
alpha = deg2rad(alpha_deg);

% Tolerancias propostas para o estudo de convergencia (erros absolutos).
tol_cl = 0.001;
tol_alpha_deg = 0.01;

fprintf('\n=== METODO DOS PERFIS DELGADOS ===\n');

%% 1. Importar geometria

fid = fopen(filename, 'r');
assert(fid ~= -1, 'Nao foi possivel abrir: %s', filename);
name = fgetl(fid);
counts = sscanf(fgetl(fid), '%f');
if numel(counts) ~= 2 || any(counts < 2) || any(counts ~= fix(counts))
	fclose(fid);
	error('Cabecalho Lednicer invalido.');
end

[upper, nu] = fscanf(fid, '%f %f', [2, counts(1)]);
[lower, nl] = fscanf(fid, '%f %f', [2, counts(2)]);
fclose(fid);
assert(nu == 2*counts(1) && nl == 2*counts(2), ...
	'Numero de coordenadas incorreto.');
upper = upper';
lower = lower';
assert(all(isfinite(upper(:))) && all(isfinite(lower(:))), ...
	'Coordenadas nao finitas.');

xu = upper(:, 1); zu = upper(:, 2);
xl = lower(:, 1); zl = lower(:, 2);

% Corrigir apenas o arredondamento dos extremos.
tol = 1e-5;
assert(all(abs([xu(1), xl(1)]) < tol) && ...
	all(abs([xu(end), xl(end)] - 1) < tol), ...
	'Verificar normalizacao das coordenadas da corda.');
xu([1 end]) = [0; 1];
xl([1 end]) = [0; 1];
assert(all(diff(xu) > 0) && all(diff(xl) > 0), ...
	'As coordenadas x devem ser estritamente crescentes.');

%% 2. Resolver o problema (unica implementacao matematica)

[I0, A1, alpha_L0, cl, data] = ...
	solveThinAirfoil(xu, zu, xl, zl, N, alpha);

% Declive teorico da curva de sustentacao.
cl_alpha_rad = 2*pi;
cl_alpha_deg = cl_alpha_rad*pi/180;

%% 3. Estudo de convergencia

I0_test = zeros(size(N_test));
A1_test = zeros(size(N_test));
alpha_L0_test = zeros(size(N_test));
cl_5_test = zeros(size(N_test));

for k = 1:numel(N_test)
	[I0_test(k), A1_test(k), alpha_L0_test(k), cl_test] = ...
		solveThinAirfoil(xu, zu, xl, zl, N_test(k), alpha);
	cl_5_test(k) = cl_test(alpha_deg == 5);
end

alpha_L0_test_deg = rad2deg(alpha_L0_test);
% Referencia numerica: maior N testado, nao a solucao exata.
error_cl = abs(cl_5_test - cl_5_test(end));
error_alpha_L0 = abs(alpha_L0_test_deg - alpha_L0_test_deg(end));
error_I0 = abs(I0_test - I0_test(end));
error_A1 = abs(A1_test - A1_test(end));

% Verificar o N usado na solucao principal.
error_cl_N = abs(cl(alpha_deg == 5) - cl_5_test(end));
error_alpha_N = abs(rad2deg(alpha_L0) - alpha_L0_test_deg(end));
converged_N = error_cl_N <= tol_cl && error_alpha_N <= tol_alpha_deg;

%% 4. Tempo de calculo

% Mesma funcao e mesmas cinco saidas da solucao principal.
% Exclui importacao, graficos, convergencia e exportacao.
solver = @() solveThinAirfoil(xu, zu, xl, zl, N, alpha);
execution_time = timeit(solver, 5);

%% 5. Apresentar e exportar resultados

% So abrir o CSV depois de os calculos terminarem com sucesso.
fid_out = fopen(outputFile, 'w');
assert(fid_out ~= -1, 'Nao foi possivel abrir o ficheiro de resultados.');
try
	fprintf(fid_out, 'kind,key,x,value,unit\n');

	fprintf('\n[1] Geometria do perfil\n');
	fprintf('%-27s : %s\n', 'Perfil', strtrim(name));
	printResult(fid_out, 'upper_points', numel(xu), '1', 'Pontos do extradorso', '%d');
	printResult(fid_out, 'lower_points', numel(xl), '1', 'Pontos do intradorso', '%d');

	fprintf('\n[2] Linha media\n');
	printResult(fid_out, 'N', N, '1', 'Pontos de interpolacao', '%d');
	printResult(fid_out, 'camber_max', max(data.zc), '1', 'Linha media maxima zc/c', '%+.4f');
	printResult(fid_out, 'camber_min', min(data.zc), '1', 'Linha media minima zc/c', '%+.4f');
	printResult(fid_out, 'thickness_max', max(data.thickness), '1', 'Espessura maxima e/c', '%.4f');

	fprintf('\n[3] Declive da linha media\n');
	printResult([], 'slope_max', max(data.dzdx), '1', 'Declive maximo', '%+.4f');
	printResult([], 'slope_min', min(data.dzdx), '1', 'Declive minimo', '%+.4f');
	printResult([], 'slope_LE', data.dzdx(1), '1', 'Declive no BA (aprox.)', '%+.4f');
	printResult([], 'slope_TE', data.dzdx(end), '1', 'Declive no BF (aprox.)', '%+.4f');

	fprintf('\n[4] Mudanca de coordenadas\n');
	printResult([], 'N_theta', numel(data.theta), '1', 'Pontos de integracao', '%d');

	fprintf('\n[5] Coeficientes de Fourier\n');
	printResult(fid_out, 'I0', I0, '1', 'I0', '%+.6f');
	printResult(fid_out, 'A1', A1, '1', 'A1', '%+.6f');
	printResult([], 'A0_alpha0', -I0, '1', 'A0 para alpha = 0', '%+.6f');

	fprintf('\n[6] Resultados aerodinamicos\n');
	printResult(fid_out, 'alpha_L0_deg', rad2deg(alpha_L0), 'deg', 'Angulo sustentacao nula', '%+.4f');
	printResult(fid_out, 'cl_alpha_rad', cl_alpha_rad, 'rad^-1', 'Declive da curva', '%.4f');
	printResult(fid_out, 'cl_alpha_deg', cl_alpha_deg, 'deg^-1', 'Declive da curva', '%.4f');
	printResult([], 'cl_0', cl(alpha_deg == 0), '1', 'cl a 0 graus', '%+.4f');
	printResult([], 'cl_5', cl(alpha_deg == 5), '1', 'cl a 5 graus', '%+.4f');
	printResult(fid_out, 'cl_curve', cl, '1', 'Curva de sustentacao', '', alpha_deg);

	fprintf('\n[7] Estudo de convergencia\n');
	T = table(N_test(:), I0_test(:), A1_test(:), ...
		alpha_L0_test_deg(:), cl_5_test(:), error_cl(:), error_alpha_L0(:), ...
		'VariableNames', {'N', 'I0', 'A1', 'AlphaL0_deg', 'cl_5', ...
		'Error_cl', 'Error_AlphaL0_deg'});
	disp(T);
	printResult(fid_out, 'convergence_I0', I0_test, '1', '', '', N_test);
	printResult(fid_out, 'convergence_A1', A1_test, '1', '', '', N_test);
	printResult(fid_out, 'convergence_alpha_L0', alpha_L0_test_deg, 'deg', '', '', N_test);
	printResult(fid_out, 'convergence_cl_5', cl_5_test, '1', '', '', N_test);
	printResult(fid_out, 'convergence_N_ref', N_test(end), '1', 'N de referencia', '%d');
	printResult(fid_out, 'convergence_tol_cl', tol_cl, '1', 'Tolerancia cl', '%.4g');
	printResult(fid_out, 'convergence_tol_alpha_deg', tol_alpha_deg, 'deg', 'Tolerancia alpha L0', '%.4g');
	printResult(fid_out, 'convergence_error_cl_N', error_cl_N, '1', 'Erro cl a 5 graus (N)', '%.6g');
	printResult(fid_out, 'convergence_error_alpha_N', error_alpha_N, 'deg', 'Erro alpha L0 (N)', '%.6g');
	if converged_N
		fprintf('N = %d satisfaz ambas as tolerancias.\n', N);
	else
		fprintf('N = %d NAO satisfaz ambas as tolerancias.\n', N);
	end

	fprintf('\n[8] Desempenho computacional\n');
	printResult(fid_out, 'execution_time_s', execution_time, 's', 'Tempo de calculo', '%.8f');
	printResult([], 'execution_time_ms', 1000*execution_time, 'ms', 'Tempo de calculo', '%.4f');
	printResult(fid_out, 'alpha_count', numel(alpha), '1', 'Angulos de ataque calculados', '%d');
catch ME
	fclose(fid_out);
	rethrow(ME);
end
fclose(fid_out);
fprintf('\nResultados guardados em:\n%s\n', outputFile);

%% 6. Visualizacao e exportacao vetorial (fora da medicao de tempo)

if ~isfolder(figDir)
	mkdir(figDir);
end

% Geometria e declive com coordenada horizontal partilhada.
figGeometry = figure('Name', 'Geometria e declive', ...
	'Color', 'w', 'WindowStyle', 'docked');
layout = tiledlayout(figGeometry, 2, 1, ...
	'TileSpacing', 'compact', 'Padding', 'compact');

axGeometry = nexttile(layout);
hold(axGeometry, 'on');
plot(axGeometry, data.x, data.zu_interp, 'b-', 'LineWidth', 1.2);
plot(axGeometry, data.x, data.zl_interp, 'r-', 'LineWidth', 1.2);
plot(axGeometry, xu, zu, 'b.', 'MarkerSize', 10);
plot(axGeometry, xl, zl, 'r.', 'MarkerSize', 10);
plot(axGeometry, data.x, data.zc, 'k-', 'LineWidth', 1.6);
yline(axGeometry, 0, 'k--');
hold(axGeometry, 'off');
axis(axGeometry, 'equal');
pbaspect(axGeometry, [5 1 1]);
grid(axGeometry, 'on');
xlim(axGeometry, [0 1]);
ylim(axGeometry, 'auto');
xlabel(axGeometry, 'x/c');
ylabel(axGeometry, 'z/c');
title(axGeometry, 'Geometria e linha media');
legend(axGeometry, 'Extradorso interpolado', 'Intradorso interpolado', ...
	'Pontos do extradorso', 'Pontos do intradorso', ...
	'Linha media', 'Linha de corda', ...
	'Location', 'northoutside', 'NumColumns', 3);

axSlope = nexttile(layout);
plot(axSlope, data.x, data.dzdx, 'k-', 'LineWidth', 1.5);
yline(axSlope, 0, 'k--');
grid(axSlope, 'on');
xlim(axSlope, [0 1]);
xlabel(axSlope, 'x/c');
ylabel(axSlope, 'dz_c/dx');
title(axSlope, 'Declive da linha media');

% Sincronizar apenas x: as escalas verticais sao independentes.
linkaxes([axGeometry, axSlope], 'x');
xlim(axGeometry, [0 1]);
drawnow;
exportgraphics(axGeometry, ...
	fullfile(figDir, 'geometry.pdf'), ...
	'ContentType', 'vector');
exportgraphics(axSlope, ...
	fullfile(figDir, 'camber_slope.pdf'), ...
	'ContentType', 'vector');

% A curva individual fica opcional. Os 81 pontos continuam no CSV.
if showLiftCurve
	figure('Name', 'Curva de sustentacao', 'Color', 'w');
	plot(alpha_deg, cl, 'b-', 'LineWidth', 1.5);
	xline(0, 'k--'); yline(0, 'k--'); grid on;
	xlabel('\alpha [deg]'); ylabel('c_l');
	title('Roncz - Teoria dos perfis delgados');
end

% Convergencia numerica.
figConvergence = figure('Name', 'Convergencia', ...
	'Color', 'w', 'WindowStyle', 'docked');
layoutConv = tiledlayout(figConvergence, 2, 1, ...
	'TileSpacing', 'compact', 'Padding', 'compact');
axCl = nexttile(layoutConv);
semilogx(axCl, N_test, cl_5_test, 'bo-', 'LineWidth', 1.5);
grid(axCl, 'on');
ylabel(axCl, 'c_l a 5 graus');
title(axCl, 'Convergencia da sustentacao');
axAlpha = nexttile(layoutConv);
semilogx(axAlpha, N_test, alpha_L0_test_deg, 'ro-', 'LineWidth', 1.5);
grid(axAlpha, 'on');
xlabel(axCl, 'N');
ylabel(axAlpha, '\alpha_{L=0} [deg]');
title(axAlpha, 'Convergencia do angulo de sustentacao nula');
linkaxes([axCl, axAlpha], 'x');
drawnow;
exportgraphics(axCl, ...
    fullfile(figDir, 'convergence_cl.pdf'), ...
    'ContentType', 'vector');
exportgraphics(axAlpha, ...
    fullfile(figDir, 'convergence_alpha_L0.pdf'), ...
    'ContentType', 'vector');

fprintf('\nGraficos guardados em:\n%s\n', figDir);
