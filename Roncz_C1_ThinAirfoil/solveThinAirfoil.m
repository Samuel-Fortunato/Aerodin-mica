function [I0, A1, alpha_L0, cl, data] = ...
    solveThinAirfoil(xu, zu, xl, zl, N, alpha)
% SOLVETHINAIRFOIL Calcula a sustentacao pela teoria dos perfis delgados.
% Coordenadas normalizadas pela corda; alpha e alpha_L0 em radianos.
% data contem os resultados intermedios para graficos e exportacao.

    assert(isscalar(N) && N >= 3 && N == fix(N), ...
        'N deve ser um inteiro >= 3.');

    %% Linha media
    x = linspace(0, 1, N);
    s = sqrt(x);
    zu_interp = interp1(sqrt(xu), zu, s, 'pchip');
    zl_interp = interp1(sqrt(xl), zl, s, 'pchip');
    zc = (zu_interp + zl_interp)/2;
    thickness = zu_interp - zl_interp;

    assert(all(isfinite(zc)) && all(thickness >= -1e-10), ...
        'Geometria interpolada invalida.');

    %% Declive em relacao a x/c
    dzdx = gradient(zc, x);

    %% Mudanca de coordenadas para integracao
    theta = linspace(0, pi, N);
    x_theta = (1 - cos(theta))/2;
    dzdx_theta = interp1(x, dzdx, x_theta, 'pchip');
    assert(all(isfinite(dzdx_theta)), 'Declive invalido em theta.');

    %% Coeficientes de Fourier
    I0 = trapz(theta, dzdx_theta)/pi;
    A1 = (2/pi)*trapz(theta, dzdx_theta.*cos(theta));
    A0 = alpha - I0;

    %% Resultados aerodinamicos
    alpha_L0 = I0 - A1/2;
    cl = 2*pi*(A0 + A1/2);

    %% Resultados intermedios (opcional)
    if nargout >= 5
        data.x = x;
        data.zu_interp = zu_interp;
        data.zl_interp = zl_interp;
        data.zc = zc;
        data.thickness = thickness;
        data.dzdx = dzdx;
        data.theta = theta;
        data.x_theta = x_theta;
        data.dzdx_theta = dzdx_theta;
        data.A0 = A0;
    end
end
