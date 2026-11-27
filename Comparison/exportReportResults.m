function exportReportResults()
% Exportar valores numericos do relatorio sem executar os solvers.
root = fileparts(fileparts(mfilename('fullpath')));
thin = readtable(fullfile(root,'thin-airfoil','results.csv'),'TextType','string');
panel = readtable(fullfile(root,'panel-method','results.csv'));
conv = readtable(fullfile(root,'panel-method','convergence.csv'));
cases = readtable(fullfile(root,'panel-method','cases.csv'));
cp = readtable(fullfile(root,'panel-method','cp.csv'));
cfd = readtable(fullfile(root,'CFD','results.csv'));
comp = readtable(fullfile(root,'Comparison','results.csv'),'TextType','string');
scalar = @(key) thin.value(thin.kind=="scalar" & thin.key==key);
curve = thin(thin.key=="cl_curve",:);

% Dados geometricos e convergencia dos perfis delgados.
v = struct();
v.PontosExtradorso = scalar("upper_points");
v.PontosIntradorso = scalar("lower_points");
v.NPerfisDelgados = scalar("N");
v.NRefPerfisDelgados = scalar("convergence_N_ref");
v.EspessuraPct = 100*scalar("thickness_max");
v.CurvaturaMaxPct = 100*scalar("camber_max");
v.CurvaturaMinPct = 100*scalar("camber_min");
v.ErroClPerfisDelgados = scalar("convergence_error_cl_N");
v.ErroAlphaPerfisDelgados = scalar("convergence_error_alpha_N");
v.TolClPerfisDelgados = scalar("convergence_tol_cl");
v.TolAlphaPerfisDelgados = scalar("convergence_tol_alpha_deg");
tc = sortrows(thin(thin.key=="convergence_cl_5",:),'x');
ta = sortrows(thin(thin.key=="convergence_alpha_L0",:),'x');
assert(isequal(tc.x,ta.x),'Malhas de convergencia diferentes.');
v.NMinTestePerfisDelgados = min(tc.x);
v.NMinConvergenciaPerfisDelgados = min(tc.x( ...
    abs(tc.value-tc.value(end))<=v.TolClPerfisDelgados & ...
    abs(ta.value-ta.value(end))<=v.TolAlphaPerfisDelgados));

% Discretizacao realmente usada nos casos de pressao.
nPanels = arrayfun(@(a) sum(cp.alpha_deg==a),cases.alpha_deg);
assert(all(nPanels==nPanels(1)) && nPanels(1)>0);
v.NPaineis = nPanels(1);
v.NIncognitasPaineis = v.NPaineis+1;
v.NPaineisOriginais = v.PontosExtradorso+v.PontosIntradorso-2;
conv = sortrows(conv,'N');
idx = find(conv.N==v.NPaineis);
assert(isscalar(idx) && idx>1);
v.NMinTestePaineis = conv.N(1);
v.PassoTestePaineis = conv.N(2)-conv.N(1);
v.NRefPaineis = conv.N(end);
v.NAnteriorPaineis = conv.N(idx-1);
v.ClRefPaineis = conv.cl(end);
v.ClConvergenciaPaineis = conv.cl(idx);
v.ErroAnteriorPaineis = conv.relative_error_pct(idx-1);
v.ErroConvergenciaPaineis = conv.relative_error_pct(idx);
v.TempoRefPaineis = conv.time_s(end);
v.TempoConvergenciaPaineis = conv.time_s(idx);
v.ReducaoTempoPaineis = 100*(1-v.TempoConvergenciaPaineis/v.TempoRefPaineis);

% Coeficientes de sustentacao nos angulos comuns.
angles = [0 5 10 15 20];
suffix = ["Zero","Cinco","Dez","Quinze","Vinte"];
for k = 1:numel(angles)
    a = angles(k);
    s = suffix(k);
    v.("ClDelgados"+s) = curve.value(curve.x==a);
    v.("ClPaineis"+s) = panel.cl(panel.alpha_deg==a);
    v.("ClCFD"+s) = cfd.cl(cfd.aoa==a);
    v.("DifClPaineisDelgados"+s) = v.("ClPaineis"+s)-v.("ClDelgados"+s);
    v.("DifClPaineisCFD"+s) = v.("ClPaineis"+s)-v.("ClCFD"+s);
    if k>1
        v.("IncrementoClCFD"+s) = v.("ClCFD"+s)-v.("ClCFD"+suffix(k-1));
    end
end
v.CdCFDVinte = cfd.cd(cfd.aoa==20);
v.ClMaxCFD = max(cfd.cl);
v.NumAngulosCFD = height(cfd);
v.ReducaoDecliveCFD = 100*(1-v.IncrementoClCFDDez/v.IncrementoClCFDCinco);

% Regressao linear comum aos tres metodos.
names = ["Perfis delgados","Paineis","CFD"];
prefix = ["Delgados","Paineis","CFD"];
for k=1:numel(names)
    r = comp(comp.methods==names(k),:);
    v.("AlphaZero"+prefix(k)) = r.alpha_L0_deg;
    v.("Declive"+prefix(k)) = r.slope_rad;
end
v.DifDeclivePct = 100*(v.DeclivePaineis/v.DecliveDelgados-1);
v.LimiteLinearCFD = comp.linear_max_deg(comp.methods=="CFD");

% Casos especificos: somente os valores disponiveis no CSV.
for k=[1 2 5]
    r = cases(cases.alpha_deg==angles(k),:);
    s = suffix(k);
    v.("EstagnacaoPaineis"+s) = r.x_stag;
    v.("CpMinPaineis"+s) = r.cp_min;
    v.("XCpMinPaineis"+s) = r.x_cp_min;
end

% Preservar a precisao; o arredondamento e feito pelo LaTeX.
lines = ["% Gerado automaticamente por Comparison/exportReportResults.m."; ...
         "% Nao editar: executar exportReportResults para atualizar a partir dos CSV."];
keys = fieldnames(v);
for k=1:numel(keys)
    key = keys{k};
    value = v.(key);
    assert(isscalar(value) && isfinite(value),'Resultado em falta: %s',key);
    lines(end+1,1) = string(sprintf('\\newcommand{\\%s}{%.15g}',key,value));
end
out = fullfile(root,'relatorio','config','resultados_calculo.tex');
tmp = [tempname(fileparts(out)) '.tex'];
writelines(lines,tmp);
[ok,msg] = movefile(tmp,out,'f');
assert(ok,'Erro ao exportar resultados LaTeX: %s',msg);
end
