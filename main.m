%% Computer Project 2
% BIOENG 1005
% Paul Kullmann


delta_z = 1.5*10e-2;
k = 1:1000;

main_fcn;

function main_fcn()

    E = zeros(1,1000);
    H = zeros(1,999);

    t_step = 2.5*10e-11;
    z_step = 1.5*10e-2;

    num_time_steps = numel(0:t_step:22.5*10e-9);
    num_z_steps = 1000;


    source = zeros(1,num_time_steps);
    T = 10e-9;
    t_vec = (0:num_time_steps-1) * t_step;
    mask = (t_vec >= 0) & (t_vec <= 6*T);
    source(mask) = exp(-((t_vec(mask) - 3*T) / (T^2)));



    for t = 1:num_time_steps

        %% Update H field:
        for k = 1:num_z_steps-1 
            H(k) = -(t_step/mu(k)/z_step)*(E(k+1)-E(k-1)) + H(k+1);
        end

        %% Update E field:
        for k = 2:num_z_steps-1
            E(k) = -(1/z_step)*(H(k+1)-H(k-1))/(epsilon(k)/t_step + sigma(k)/2) + E(k-1)*(epsilon(k)/delta_t - sigma(k)/2)/(epsilon(k)/delta_t + sigma(k)/2);
        end

        %% Source:
        E(1) = E(1) + source(t);


        %% Boundary Conditions:
        E(num_z_steps) = 0;

    end
    

end

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

function E = pulse(t)
    % Return E at location 0 as the gaussian pulse for 0-6T

    T = 10e-9;

    if t > 6*T
        E = 0;
        return
    elseif t >=0 && t <= 6*T
        E = exp(-(t-3*T)^2/T^2);
        return
    end
end