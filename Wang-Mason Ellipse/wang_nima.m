% Wang/Mason 2D impact model
% mu is the coefficient of friction
% epsilon is the coefficient of restitution
function [v_plus, z] = wang_nima(M, n, s, v, ha, mu, epsilon)

  % update v
  v = v + ha;

  % setup the impulse
  z = [0 0]';

  % setup dt
  dt = 1e-4;

  % setup near zero
  NEAR_ZERO = 1e-8;

  % setup v_plus
  v_plus = v;

  % keep integrating until the normal velocity is zero 
  while (n'*v_plus < 0)

    % get tangent velocity 
    if (abs(s'*v_plus) < NEAR_ZERO)
      % compute the impulse
      dV = [dt 0]';
      dz = [n s] \ M*[n s]*dV;  
%dV = [n s]*M \ [n s]' * dz

      % tangent velocity is effectively zero; see whether impulse is within
      % friction cone 
      if (mu*dz(1) < abs(dz(2)))
        % not within the friction cone; get the direction of impending
        % slip
        dz(2) = mu*dz(1);
        dV = M\[n s]*dz;

        % reverse direction if necessary
        if (dV(2)*dz(2) > 0)
          dz(2) = -dz(2);
        end
      end

      % update z
      z = z + dz;
    else
      % sliding update to z is easy
      z(1) = z(1) + dt;
      z(2) = z(2) - sign(s'*v_plus)*dt*mu;
    end

    % update v_plus
    v_plus = v + M \ [n s] * z; 
  end

  % continue integrating
  zstar = z;
  while (zstar(1) < z(1)*(1+epsilon))
    % get tangent velocity 
    if (abs(s'*v_plus) < NEAR_ZERO)
      % compute the impulse
      dV = [dt 0]';
      dz = [n s]\M*[n s]*dV;  

      % tangent velocity is effectively zero; see whether impulse is within
      % friction cone 
      if (mu*dz(1) < abs(dz(2)))
        % not within the friction cone; get the direction of impending
        % slip
        dz(2) = mu*dz(1);
        dV = M\[n s]*dz;

        % reverse direction if necessary
        if (dV(2)*dz(2) > 0)
          dz(2) = -dz(2);
        end
      end

      % update zstar
      zstar = zstar + dz;
    else
      zstar(1) = zstar(1) + dt;
      zstar(2) = zstar(2) - sign(s'*v_plus)*dt*mu;
    end

    % update v_plus
    v_plus = v + M \ [n s] * zstar; 
  end

  % setup z
  z = zstar;
  z(3) = 0;
  if (z(2) < 0)
    z(3) = -z(2);
    z(2) = 0;
  end
  z(4) = abs(s'*v_plus);

