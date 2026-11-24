function [u_uniforme,w_uniforme] = EsUniforme(V_inf_modul,AoA_rad)
% =========================================================================
% FUNÇÃO: EsUniforme
% =========================================================================
% Calcula as componentes do campo de velocidades (u, w) de um escoamento
% uniforme constante em todo o domínio 2D.
%
% Referência Teórica: Katz, J., & Plotkin, A. (2001). Low-Speed Aerodynamics.
%
% PARÂMETROS DE ENTRADA:
%   V_inf_modul : Módulo da velocidade do escoamento livre [m/s]
%   AoA_rad     : Ângulo de ataque (Angle of Attack, \alpha) em radianos [rad]
%
% PARÂMETROS DE SAÍDA:
%   u_uniforme  : Componente horizontal da velocidade (eixo X) [m/s]
%   w_uniforme  : Componente vertical da velocidade (eixo Z) [m/s]
%
% FUNDAMENTAÇÃO MATEMÁTICA:
%   Um escoamento uniforme é caracterizado por ter velocidade vetorial constante
%   em todos os pontos da malha. Projetando o módulo V_inf nos eixos cartesianos:
%       u = V_inf * cos(AoA)
%       w = V_inf * sin(AoA)
% =========================================================================

% Construção do vetor velocidade 2D [Vx; Vz] através da decomposição trigonométrica
V_inf = V_inf_modul * [cos(AoA_rad); sin(AoA_rad)];

% Componente horizontal da velocidade (eixo X)
u_uniforme = V_inf(1);

% Componente vertical da velocidade (eixo Z)
w_uniforme = V_inf(2);
end