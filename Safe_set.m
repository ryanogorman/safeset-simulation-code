function poly_safe_set = Safe_set(ego, lead, road, i)
% Inputs:
% ego: structure that contains all data for the ego vehicle, including
% dmin, Mode, Option, MPC params, reference velocity, and all useful 
% physical properties.
% lead: structure that contains all data for the lead vehicle
% road: structure that contains all data for the road, in particular
% position and grade arrays
% Outputs:
% poly_safe_set: fitted polynomial coefficients for determining guaranteed 
% safe following distance

N = ego.Params.N;
% presetting variables b/c of some weird simulink thing
% worth checking if I can remove this later
sLead_predict = zeros(1,N+1);
S1_emg = zeros(1,300);
S2_emg = zeros(1,300);
V1_emg = zeros(1,300);
V2_emg = zeros(1,300);
dsafe = zeros(1,300);
%
poly_safe_set = zeros(N+1,3); % Initializing variable. This is output when time is not a multiple of 0.2
delta_t = 1 / ego.Params.freq_MPC; % Time step used to calculate safe-set
dmin = ego.l/2 + lead.l/2 + ego.dmin;
g = 9.81;

time = (i - 1)/ego.Params.freq_sim;
mult = ego.Params.freq_sim/ego.Params.freq_MPC;

if abs(floor(time*ego.Params.freq_MPC) - time*ego.Params.freq_MPC) < 0.001 
    
    % The OPTION will determine how much we know about the lead car's
    % velocity profile.
    if ego.OPTION == 1
        % The exact future velocity is known
        V01 = lead.velocity_profile(i:mult:i+mult*10); % preview of v1 for i to i+Np during the prediction horizon
    elseif ego.OPTION == 2
        % Current velocity known, future velocity assumed constant.
        V01 = lead.velocity(i)*ones(1,N+1);
    % possibility for other OPTIONs aswell
    else
        error('Invalid Option')
    end
    
    sLead_predict(1) = lead.states.pos; % come back to this later
    for ti = 1:N
        sLead_predict(ti+1) = sLead_predict(ti) + V01(ti) * delta_t; % recall from MPC code that s(k+1) = s(k) + v(k)*delta_t
    end
    S01 = sLead_predict; % preview of s1 for i to i+Np during the prediction horizon
    % emergency
    for j = 1:N+1
        k = 1;
        S1_emg(1) = S01(j);
        V1_emg(1)= V01(j); 
        dt2 = delta_t / 10; % discretization time for computing the safe set
        while V1_emg(k) >=0 % while the car has forward motion
            if abs(ego.MODE - 1) <= 0.001
                a01 = interp1(road.position, road.Grade, S1_emg(k), 'linear', 0); % road angle (rad)
                a01 = atan(a01/100);
            else
                a01 = 0;
            end
            V1_emg(k+1)= V1_emg(k)+dt2/lead.m*(lead.Umin-lead.kr*lead.m*g*cos(a01)-lead.ka*V1_emg(k)^2-lead.m*g*sin(a01));
            S1_emg(k+1) = S1_emg(k)+dt2*V1_emg(k); 
            k = k+1;
        end
        % K IS WHEN THE FRONT CAR STOPS
       h = 1;
       S2_emg(h) = S1_emg(k)- dmin; % S2_emg(1) is where the ego car must be when it comes to a complete stop
       V2_emg(h) = 0;
       dsafe(h) = (S1_emg(1)-S2_emg(h)); % safe distance
       % backward integration of the follower vehicle EOM
       while V2_emg(h) < ego.Vmax 
        if abs(ego.MODE - 1) <= 0.001
            a02 = interp1(road.position, road.Grade, S2_emg(h),'linear', 0); % road angle (rad) for ego car
            a02 = atan(a02);
        else
            a02 = 0;
        end
        V2_emg(h+1)= V2_emg(h)-dt2/ego.m*(ego.Umin-ego.kr*ego.m*g*cos(a02)-ego.ka*V2_emg(h)^2-ego.m*g*sin(a02));
        S2_emg(h+1) = S2_emg(h)-dt2* V2_emg(h);
        h = h+1; 
        dsafe(h) = (S1_emg(1)-S2_emg(h));
       end
%        figure(11);hold all
%        plot(dsafe,V2_emg) % this is what we are fitting a polynomial to
       % "hold" is just a way of getting rid of all the trivial entries at
       % the end of V2_emg and dsafe. I did this to work around some
       % simulink issue, so don't worry about it too much.
       hold = max(find(V2_emg));
        poly_safe_set(j,:) = polyfit(V2_emg(1:hold),dsafe(1:hold),2); 
    end
end
