function printResult(fid, key, value, unit, label, fmt, x)
% Imprime um resultado na consola e exporta-o para CSV.
% unit = '1' identifica grandezas adimensionais.
% Sem x: escalar. Com x: serie. label vazio: exportacao silenciosa.

    if nargin < 5, label = key; end
    if nargin < 6, fmt = '%.6f'; end

    if nargin < 7
        assert(isscalar(value), 'O resultado deve ser escalar.');
        if ~isempty(label)
            suffix = '';
            if ~strcmp(unit, '1') && ~isempty(unit)
                suffix = [' ' unit];
            end
            fprintf('%-27s : %s%s\n', label, sprintf(fmt, value), suffix);
        end
        fprintf(fid, 'scalar,%s,,%.15g,%s\n', key, value, unit);
    else
        assert(numel(x) == numel(value), ...
            'Os vetores devem ter o mesmo tamanho.');
        for k = 1:numel(x)
            fprintf(fid, 'series,%s,%.15g,%.15g,%s\n', ...
                key, x(k), value(k), unit);
        end
        if ~isempty(label)
            fprintf('%-27s : %d pontos exportados\n', label, numel(x));
        end
    end
end
