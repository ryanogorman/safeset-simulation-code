function controller_ACC = MPC_setup(Params, ego, lead_l)
%%% Author: Ryan O'Gorman <ryanogorman@berkeley.edu> %%%
% Params structure will contain the fields:
% N: MPC horizon length
% freq_MPC: frequency of MPC (Hz)
% freq_sim: frequency of dynamics (Hz)
% Q: Value punishing velocity deviance from reference velocity
% P: Terminal cost for velocity deviance
% R_i: Value punishing input force
% R_j: Value punishing jerk
% Ego structure will contain all data for the ego vehicle

if mod(Params.freq_sim/Params.freq_MPC,1)
   error('freq_sim must be divisible by freq_MPC')
end

delta_t = 1 / Params.freq_MPC;
g = 9.81;

Q = Params.Q;
R = Params.R;
RR = Params.RR;
P = Params.P;
N = Params.N;

nx = 1; % Number of states
nu = 1; % Number of inputs
% Control input u is based on (kN)
ny = 1; % Number of outputs

% MPC data
num_poly_coeff = 3; % 3 is the number of polynomial coefficients

u = sdpvar(repmat(nu,1,N),repmat(1,1,N));
x = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
r = sdpvar(repmat(ny,1,N+1),repmat(1,1,N+1));
s = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));

theta_var = sdpvar(repmat(ny,1,N+1),repmat(1,1,N+1));
sLead_var = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
vLead_var = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
poly_safe_set_var = sdpvar(N+1,num_poly_coeff);
poly_safe_set_var1 = sdpvar(N+1,num_poly_coeff);
U_i = sdpvar(1);

%epsi = sdpvar(1);
constraints = [];
objective = 0;

% The following are based on basic discrete non-linearized kinematics
% -- i.e. s(k+1) = s(k) + x(k)*delta_t
% & x(k+1) = x(k) + sum(Forces(k))/m*delta_t
% However, position 's' is defined as the distance to the collision
% point where the paths of the two vehicles meet. Thus, if either
% vehicle is 20 meters away from the collision point its position is 
% s = 20, if it is 20 meters past the collision point its position is
% s = -20. As a result, the kinematic equations are slightly modified
% so that s(k+1) = s(k) - x(k)*delta_t <-- Note the change in sign
for k = 1:N
    if k == 1
        objective = objective + (x{k}-vLead_var{k})'*Q*(x{k}-vLead_var{k}) + u{k}'*R*u{k} + (u{k}-U_i)'*RR*(u{k}-U_i);
    else
        objective = objective + (x{k}-vLead_var{k})'*Q*(x{k}-vLead_var{k}) + u{k}'*R*u{k} + (u{k}-u{k-1})'*RR*(u{k}-u{k-1});
    end
    constraints = [constraints, s{k+1} == s{k} + delta_t*x{k}];
    constraints = [constraints, x{k+1} == x{k}+ delta_t./ego.m*(u{k}*1000 - ego.ka*x{k}^2 - ego.m*g*ego.kr*cos(theta_var{k})-ego.m*g*sin(theta_var{k}))];
        
    constraints = [constraints, ego.Vmin<=x{k+1}<=ego.Vmax]; 
    constraints = [constraints,  u{k}<=ego.Umax/1000, ego.Umin/1000<=u{k}]; 
    %constraints = [constraints, epsi >= 0]; 
       
%         if ((poly_safe_set_var(k+1,1)*(x{k+1})^2+poly_safe_set_var(k+1,2)*(x{k+1})+poly_safe_set_var(k+1,3))) >= 4.425
%             constraints = [constraints, poly_safe_set_var1(k,:) == poly_safe_set_var(k,:)];
%         else
%             constraints = [constraints, poly_safe_set_var1(k,1:2) == 0, poly_safe_set_var1(k,3) == 4.425];
%         end
        
    %constraints = [constraints, ((poly_safe_set_var1(k+1,1)*(x{k+1})^2+poly_safe_set_var1(k+1,2)*(x{k+1})+poly_safe_set_var1(k+1,3))) >= 4.425];
    constraints = [constraints, sLead_var{k+1} == sLead_var{k} + delta_t*vLead_var{k}];
    constraints = [constraints , (sLead_var{k+1} - s{k+1}) >= (poly_safe_set_var1(k+1,1)*(x{k+1})^2+poly_safe_set_var1(k+1,2)*(x{k+1})+poly_safe_set_var1(k+1,3))];
        
    constraints = [constraints, (sLead_var{k+1} - s{k+1}) <= ego.l/2 + lead_l/2 + 10];
    constraints = [constraints, (sLead_var{k+1} - s{k+1}) >= ego.l/2 + lead_l/2 + 0.1]; % Hard collision constraint
end
  
objective = objective + (x{N+1}-vLead_var{N+1})'*P*(x{N+1}-vLead_var{N+1}); % Terminal cost
    
parameters_in = {x{1}, s{1}, [r{:}], [vLead_var{:}], sLead_var{1}, [theta_var{:}], poly_safe_set_var, U_i};

solutions_out = {[u{:}], [x{:}], [s{:}], [sLead_var{:}]};
% ipopt is currently being used as the solver
controller_ACC = optimizer(constraints, objective, sdpsettings('solver','fmincon','verbos',0),parameters_in,solutions_out);

