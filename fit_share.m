function [omega_cg_t,omega_cl_t,omega_ch_t] = fit_share(year,omega_cg_t_raw,omega_cl_t_raw,omega_ch_t_raw)

    % Perform linear regression (fit a line)
    for i=1:9
        coefficients = polyfit(year, omega_cg_t_raw(i,:), 1);
        omega_cg_t(i,:) = polyval(coefficients, year);
    end

    for i=1:9
        coefficients = polyfit(year, omega_cl_t_raw(i,:), 1);
        omega_cl_t(i,:) = polyval(coefficients, year);
    end

    for i=1:9
        coefficients = polyfit(year, omega_ch_t_raw(i,:), 1);
        omega_ch_t(i,:) = polyval(coefficients, year);
    end

end

