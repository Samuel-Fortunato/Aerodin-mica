function printResult(fid, key, value, unit, label, fmt, x)
% Imprime um resultado na consola e/ou exporta-o para CSV.
% fid vazio: apenas imprimir.
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
        if ~isempty(fid)
            fprintf(fid, 'scalar,%s,,%.15g,%s\n', key, value, unit);
        end
    else
        assert(numel(x) == numel(value), ...
            'Os vetores devem ter o mesmo tamanho.');
        if ~isempty(fid)
            for k = 1:numel(x)
                fprintf(fid, 'series,%s,%.15g,%.15g,%s\n', ...
                    key, x(k), value(k), unit);
            end
        end
        if ~isempty(label)
            if isempty(fid)
                fprintf('%-27s : %d pontos\n', label, numel(x));
            else
                fprintf('%-27s : %d pontos exportados\n', label, numel(x));
            end
        end
    end
end
