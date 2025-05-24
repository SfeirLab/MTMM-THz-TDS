function [t_cs] = MTMM(d, lambda0, theta0, nr, ns, flag, dlimit, nk)
    % MTMM - Modified Transfer Matrix Method for s-polarized transmission
    N = numel(d);                          % Total number of layers
    t_cs = zeros(numel(lambda0), numel(theta0)); % Preallocate output
    idx = find(isnan(nk));                % Locate the sample layer (NaN)

    for a = 1:numel(lambda0)
        k0 = 2 * pi / lambda0(a);         % Free-space wave number

        % Assign refractive index of sample or reference
        if flag == 0
            n_s = nr;
        else
            n_s = ns(a);
        end

        % Construct full refractive index profile for this lambda
        n = nk;
        n(idx) = n_s;

        M_sPol = eye(2);                  % Initialize transfer matrix

        for c = 1:N-1
            k_x = n(c) * k0;              % Wavevector in layer c
            phi = k_x * d(c);             % Phase delay

            % Default Transfer matrix (with reflection and transmission)
            D = [(n(c) + n(c+1)) / (2 * n(c)), (n(c) - n(c+1)) / (2 * n(c));
                 (n(c) - n(c+1)) / (2 * n(c)), (n(c) + n(c+1)) / (2 * n(c))];

            % Suppress reflection if the layer is too thick (beyond dlimit)
            if d(c+1) > dlimit(c+1)
                D = [(n(c) + n(c+1)) / (2 * n(c)), 0;
                     (n(c) - n(c+1)) / (2 * n(c)), 0];
            end

            % Propagation matrix for layer c
            if c == 1
                P = eye(2);               % No propagation in semi-infinite first medium
            else
                P = [exp(1i * phi), 0; 0, exp(-1i * phi)];
            end

            % Update total transfer matrix
            M_sPol = M_sPol * P * D;
        end

        % Transmission coefficient
        t_cs(a, :) = 1 / M_sPol(1, 1);
    end
end
