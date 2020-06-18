function out = cost_function(x, Mass, M, n ,d, pre, J, Pt)
    %assign mu and epsilon
    mu = x(1);
    epsilon =x(2);

    %run APPoisson Juniors to get predicted velocity
    v_calc = APPoisson_juniors(Mass, n, d, pre, mu, epsilon);
    
    %do error calculations
    p_hat = (M*(J*(v_calc - pre')))'; %compute error in contanct point momentum
    out = norm(Pt - p_hat'); %norm of error
end