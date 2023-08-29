% Pure longitudinal force FX0
% this function remap the sclar function to its vectorial form
function [Calfa_vec] = MF96_CorneringStiffness2(kappa_vec, alpha_vec, phi_vec, Fz_vec, tyre_data)

  
  Calfa_vec = zeros(size(kappa_vec));
  for i = 1:length(kappa_vec)
   % precode
   [alpha__y, Bx, Cx, Dx, Ex, SVy] = MF96_FY0_coeffs(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);
   % main code
    Calfa_vec(i) = magic_formula_stiffness(alpha__y, Bx, Cx, Dx, Ex, SVy);
  end
  
 end
