%% Computer Project 2
% BIOENG 1005
% Paul Kullmann


delta_z = 1.5*10e-2;
delta_t = 2.5*10e-11;
k = 1:1000;

function bone(delta_z, delta_t)
    % Compute FDTD for bone marrow slab:

    % Parameters:
    mu = 4*pi*10e-7;
    epsilon = 7.2103*(8.854*10e-12);
    sigma = 0;

    mu_v, epsilon_v, sigma_v = createGeometry(mu, epsilon, sigma);

end

function fat(delta_z, delta_t)
    % Compute FDTD for fat slab:

    % Parameters:
    mu = 4*pi*10e-7;
    epsilon = 5.9215*(8.854*10e-12);
    sigma = 0.03687;

    mu_v, epsilon_v, sigma_v = createGeometry(mu, epsilon, sigma);


end

function bladder(delta_z, delta_t)
    % Compute FDTD for bladder:

    % Parameters:
    mu = 4*pi*10e-7;
    epsilon = 20.093*(8.854*10e-12);
    sigma = 0.31684;

    mu_v, epsilon_v, sigma_v = createGeometry(mu, epsilon, sigma);
end

function mu_v, epsilon_v, sigma_v = createGeometry(mu, epsilon, sigma);
    % Create the 1D FDTD geometry vectors from the given values
    mu_0 = 4*pi*10e-7;
    epsilon_0 = 8.854*10e-12;
    
    mu_v(1:499) = mu_0;
    mu_v(500:600) = mu;
    mu_v(601:1000) = mu_0;

    epsilon_v(1:499) = epsilon_0;
    epsilon_v(500:600) = epsilon;
    epsilon_v(601:1000) = epsilon_0;

    sigma_v(1:499) = 0;
    sigma_v(500:600) = sigma;
    sigma_v(601:1000) = 0;
   
end