function exportResults()
    panelModeOverride = 3;
    AirfoilRoncz_main;

    T = table(angulos(:), CL_a(:), tempo_a(:), ...
        'VariableNames', {'alpha_deg', 'cl', 'time_s'});

    outputDir = fileparts(mfilename('fullpath'));
    writetable(T, fullfile(outputDir, 'results.csv'));

end
