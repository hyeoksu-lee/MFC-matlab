function [k, E] = compute_energy_spectrum(u, v, w, Lx, Lz)
    % COMPUTE_ENERGY_SPECTRUM - Computes 1D kinetic energy spectrum from 3D velocity fields
    % Inputs:
    %   u  - x-velocity field (3D array, size Nx x Ny)
    %   v  - x-velocity field (3D array, size Nx x Ny)
    %   w  - z-velocity field (3D array, size Nx x Ny)
    %   Lx - Domain length in x-direction
    %   Lz - Domain length in z-direction
    % Outputs:
    %   k  - Wavenumber array
    %   E  - Energy spectrum E(k)

    % Get grid dimensions
    [Nx, Nz] = size(u);

    % Mean energy
    K = mean(u.^2 + v.^2 + w.^2)/2;

    % Wavenumber increments
    dkx = 2 * pi / Lx;  % Wavenumber step in x
    dkz = 2 * pi / Lz;  % Wavenumber step in z

    % Create wavenumber grids
    kx = dkx * [0:Nx/2-1, -Nx/2+1:-1];  % Wavenumber array in x (shifted)
    kz = dkz * [0:Nz/2-1, -Nz/2+1:-1];  % Wavenumber array in y (shifted)
    [KX, KZ] = meshgrid(kx, kz);        % 2D wavenumber grid
    K = sqrt(KX.^2 + KZ.^2);            % Magnitude of wavenumber

    % Compute FFT of velocity fields
    u_hat = fft2(u);  % 2D FFT of u
    v_hat = fft2(v);  % 2D FFT of v
    w_hat = fft2(w);  % 2D FFT of w

    % Compute energy in wavenumber space (kinetic energy = 1/2 * |u|^2 + |v|^2 + |w|^2)
    E_hat = (u_hat.*conj(u_hat) + v_hat.*conj(v_hat) + w_hat.*conj(w_hat))/(Nx^2 * Nz^2);

    % Define 1D wavenumber bins
    kmax = min(max(kx), max(kz));  % Maximum wavenumber (Nyquist limit)
    k = 0:dkx:kmax;                % 1D wavenumber array
    E = zeros(size(k));            % Initialize energy spectrum

    % Bin energy into 1D spectrum by averaging over shells
    for i = 1:length(k)-1
        % Find indices where wavenumber falls in the shell [k(i), k(i+1)]
        mask = (K >= k(i)) & (K < k(i+1));
        if any(mask(:))
            E(i) = sum(E_hat(mask)) / sum(mask(:));  % Average energy in the shell
        end
    end
end