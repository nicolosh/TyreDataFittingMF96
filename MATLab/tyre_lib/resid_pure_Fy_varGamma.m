function res = resid_pure_Fy_varGamma(P,FY,ALPHA,GAMMA,FZ,tyre_data)

    % ----------------------------------------------------------------------
    %% Compute the residuals - least squares approach - to fit the Fx curve 
    %  with Fz=Fz_nom, IA=0. Pacejka 1996 Magic Formula
    % ----------------------------------------------------------------------

    % Define MF coefficients

    %Fz0 = 200*4.44822; % Nominal load 200 lbf
    
    tmp_tyre_data = tyre_data;
    % assigne computed parameter
    tmp_tyre_data.pDy3 = P(1);
    tmp_tyre_data.pEy3 = P(2);
    tmp_tyre_data.pEy4 = P(3);
    tmp_tyre_data.pHy3 = P(4);
    tmp_tyre_data.pKy3 = P(5);
    tmp_tyre_data.pVy3 = P(6);
   
        
    % Longitudinal Force (Pure Longitudinal Slip) Equations
    res = 0;
    for i=1:length(ALPHA)
       fy0  = MF96_FY0(0, ALPHA(i), GAMMA(i), FZ, tmp_tyre_data);
       res = res+(fy0-FY(i))^2;
    end
    
    % Compute the residuals
    res = res/sum(FY.^2);

end

