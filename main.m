%% Computer Project 2
% BIOENG 1005
% Paul Kullmann

figure();
xline(501,'k--'); hold on;
t = text(501, 0.5, 'Interface 1');
t.Rotation = 90;
main_fcn();



function [E, H] = main_fcn()
    

    E = zeros(1,1000);
    H = zeros(1,999);

    t_step = 2.5e-11;
    z_step = 1.5e-2;

    num_time_steps = numel(0:t_step:100e-9);
    num_z_steps = 1000;


    source = zeros(1,num_time_steps);
    T = 10e-9;
    t_vec = (0:num_time_steps-1) * t_step;
    mask = (t_vec >= 0) & (t_vec <= 6*T);
    source(mask) = exp(-((t_vec(mask) - 3*T).^2) / (T^2));
    disp(max(source));

    muX = 4*pi*1e-7;
    epsilonX = 5.9215*(8.854e-12);
    sigmaX = 0.03687;

    [mu, epsilon, sigma] = createGeometry(muX, epsilonX, sigmaX);

    for t = 1:num_time_steps

        %% Update fields:
        for k = 1:num_z_steps-1 
            % when k = 1 that means k+1/2
            H(k) = -t_step/(mu(k)*z_step) * (E(k+1) - E(k)) + H(k);
            
        end

        for k = 2:num_z_steps-1
            E(k) = (-1/z_step)*(H(k)-H(k-1))/(epsilon(k)/t_step + sigma(k)/2) + E(k) * ...
            (epsilon(k)/t_step - sigma(k)/2)/(epsilon(k)/t_step + sigma(k)/2);
        end

        %% Source:
        E(5) = source(t);


        %% Boundary Conditions:
        E(end) = 0;

        cla;
        plot(E); hold on; 
        xline(501,'k--');
        xline(601, 'k--');
        ylim([-1 1])
        yyaxis right;
        cla;
        plot(H); 
        ylim([-1 1])
        yyaxis left;
        pause(0.005)


    end
end


function [mu_v, epsilon_v, sigma_v] = createGeometry(mu, epsilon, sigma)
    % Create the 1D FDTD geometry vectors from the given values
    mu_0 = 4*pi*1e-7;
    epsilon_0 = 8.854e-12;
    
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