function exportResults(mode)

if nargin < 1
	mode = "angles";
end

mode = string(mode);
outputDir = fileparts(mfilename('fullpath'));

switch mode
	case "angles"
		T = exportAngles();

		writetable(T, ...
			fullfile(outputDir, 'results.csv'));

	case "convergence"
		T = exportConvergence();

		writetable(T, ...
			fullfile(outputDir, 'convergence.csv'));

	case "all"
		exportResults("angles");
		exportResults("convergence");
		exportResults("cases");

	case "cases"
		cases = exportCases();

		save(fullfile(outputDir, 'cases.mat'), 'cases');

		T = table( ...
			[cases.alpha_deg]', ...
			[cases.CL]', ...
			[cases.Cp_min]', ...
			[cases.x_Cp_min]', ...
			[cases.x_est]', ...
			[cases.z_est]', ...
			'VariableNames', ...
			{'alpha_deg', 'cl', 'cp_min', ...
			'x_cp_min', 'x_stag', 'z_stag'});

		writetable(T, ...
			fullfile(outputDir, 'cases.csv'));

		% Distribuicao de Cp nos pontos de controlo de cada superficie.
		writetable(pressureTable(cases), ...
			fullfile(outputDir, 'cp.csv'));

	otherwise
		error('Modo de exportacao invalido.');
end

end

%% Exportacao da curva de sustentacao

function T = exportAngles()

panelModeOverride = 3;
panelExportOnly = true;

run(fullfile(fileparts(mfilename('fullpath')), 'AirfoilRoncz_main.m'));

T = table(angulos(:), CL_a(:), tempo_a(:), ...
	'VariableNames', ...
	{'alpha_deg', 'cl', 'time_s'});

end

%% Exportacao da convergencia

function T = exportConvergence()

panelModeOverride = 2;
panelExportOnly = true;

run(fullfile(fileparts(mfilename('fullpath')), 'AirfoilRoncz_main.m'));

T = table( ...
	N_testes(:), ...
	CL_testes(:), ...
	Erro_ref_pct(:), ...
	erro_testes(:), ...
	tempos(:), ...
	'VariableNames', ...
	{'N', 'cl', 'relative_error_pct', ...
	'normal_error', 'time_s'});

end



%% Exportacao dos casos especificos

function cases = exportCases()

panelModeOverride = 4;
panelExportOnly = true;

run(fullfile(fileparts(mfilename('fullpath')), 'AirfoilRoncz_main.m'));

cases = caseResults;

end

%% Exportacao da distribuicao de pressao (sem recalcular o solver)

function T = pressureTable(cases)

blocks = cell(numel(cases),1);

for k = 1:numel(cases)
	r = cases(k);
	N = numel(r.Cp);
	assert(mod(N,2) == 0, 'Numero de paineis deve ser par.');

	% panelizationRoncz ordena o contorno: extradorso TE -> LE,
	% seguido do intradorso LE -> TE, com N/2 paineis por superficie.
	% A separacao por sinal de z seria incorreta num perfil reflex.
	surface = [repmat("upper",N/2,1); repmat("lower",N/2,1)];
	xc = (r.x_start(:) + r.x_end(:))/2;
	zc = (r.z_start(:) + r.z_end(:))/2;

	blocks{k} = table(repmat(r.alpha_deg,N,1), (1:N)', ...
		surface, xc, zc, r.Cp(:), ...
		'VariableNames', {'alpha_deg','panel_index','surface', ...
		'x_c','z_c','cp'});
end

T = vertcat(blocks{:});

end
