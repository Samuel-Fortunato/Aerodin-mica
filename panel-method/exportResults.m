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
