% This function plots the data that results from the AHRS dmu sim
% in Chapter 10 of the book 
% Aided Navigation: GPS and high rate sensors
% Jay A. Farrell, 2008, Mc Graw-Hill
% 
% This software is distibuted without a written or implied warranty. 
% The software is for educational purposes and is not intended for
% use in applications. Adaptation for applications is at the
% users/developers risk.
function [sys,x0,str,ts] = sfun_ahrs(t,x,u,flag,p)
% u is angular rate in rad/s and -f in m/s/s
persistent accel_ind P_r1 P_r2
T = 0.0058072;
d2r = pi/180;
if isempty(P_r1),
    accel_ind = 0;
    P_r1 = 0;
    P_r2 = 0;
end

% The state vector is x =[rho,xg,xa]

switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    %approx_hist_init;
    [sys,x0,str,ts]=InitSizes(u,T);

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1,
    [sys]=[];
    
  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2,
    [sys,accel_ind,P_r1,P_r2]=mdlUpdate(t,x,u,T);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    [sys]=mdlOutputs(t,x,u,accel_ind,P_r1,P_r2);

  %%%%%%%%%%%%%%%%%%%%%%%
  % GetTimeOfNextVarHit %
  %%%%%%%%%%%%%%%%%%%%%%%
  case 4,
    sys=mdlGetTimeOfNextVarHit(t,x,u);

  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
  case 9,
    sys=mdlTerminate(t,x,u);

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    error(['Unhandled flag = ',num2str(flag)]);

end

% end sfuntmpl

%
%=============================================================================
% mdlInitSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts]=InitSizes(u,T)
sizes = simsizes;

sizes.NumContStates  = 0;   
sizes.NumDiscStates  = 10;  % b, x_g, x_a
sizes.NumOutputs     = 18;  % b, [ind,P_res_1,P_res_2], x_a, Euler
sizes.NumInputs      = 6;  
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed
sys = simsizes(sizes);
x0  = zeros(10,1);  
x0(1) = 1;
x0(2) = 0;
x0(3) = 0;
x0(4) = 0;
x0(1:4) = x0(1:4)/norm(x0(1:4));
str = [];
ts  = [T 0];





%
%=============================================================================
% mdlUpdate
% Return the discrete states.
%=============================================================================
%
function [sys,stationary,P_r1,P_r2]=mdlUpdate(t,x,u,Ts)
persistent Phi Qd Fg Fa Pm Hm ym Rm a_filt T nrm_acc_avg int_T init int_w int_g
persistent sigma_xg sigma_nu_g sigma_xa w_ie Re Ra ya Ha T_a ge sigma_nu_a
persistent nrm_f_m2 nrm_f_m1 fcnt accl_flt_nrm att_init



T   = Ts;        % sample period
T_start = 3*Ts;  % waiting time prior to starting
b   = x(1:4);      % quat for rot from n to b
x_g = x(5:7);      % rad/s, 
x_a = x(8:10);      % m/s/s

w_ip_p = u(1:3);                % rad/s,   body frame angular rate
f_ip_p = -u(4:6);                % m/s/s,   body frame specific force

if t < T_start, %  isempty(init),   % initialize constants
    T_a        = T_start +1;    % last acceleration update time
    init       = 1
    att_init   = 0;
    fcnt       = 1;
    int_T      = 0;
    int_w      = zeros(3,1);
    int_g      = [0;0;0];
    a_filt     = exp(-0.1*T);   % digial filter with 10 s time constant
    sigma_nu_g = 2.2e-3;     % rad/s/rt_Hz, angle drift rate
    sigma_nu_a = 2.2e-2;     % m/s/s/rt_Hz, velocity drift rate
    lambda_g = 1/1000;       % 1/sec,   correlation time
    lambda_a = 1/1000;       % 1/sec,   correlation time
    Pxg       = 2e-6;        % rad^2/s^2, ss bias cov
    sigma_xg  = sqrt(2*lambda_g*Pxg)      % rad/s/s/rt_Hz, bias drift rate
    Pxa       = 2e-4;        % m^2/s^4, ss bias cov
    sigma_xa  = sqrt(2*lambda_a*Pxa) % m/s/s/s/rt_Hz, bias drift rate
    Phi = eye(9,9);
    Qd  = zeros(9,9);
    Fg  = -lambda_g*eye(3,3);
    Fa  = -lambda_a*eye(3,3);
    ge = 9.78;              % gravity magnitude
    Hm = zeros(1,9);
    Hm(1,3) = 1;
    ym = [1,0,0]';
    Rm = (1*pi/180)^2;      % rad^2, equiv to 1 deg, magnetometer noise
    Ha = zeros(3,9);
    Ha(1,2) = ge;
    Ha(2,1) =  -ge;
    ya = [0,0,ge]';
    Ra = (sigma_nu_a)^2;      % rad^2, equiv to 1 deg, accel noise
    
    Pp = zeros(9,9);
    Pm = Pp;
    accl_flt_nrm = zeros(2,1);
    
    nrm_acc_avg = norm(f_ip_p)-9.78;
    w_ie = 7.3e-5;          % rate of ECEF rel. inertial, rps
    Re   = 6e6;             % Earth radius, m
end
% initialize residual variance
P_r1 = 0;
P_r2 = 0;

% compute average acceleration norm
nrm_f = norm(f_ip_p)-9.78;
gyro_nrm = norm(w_ip_p);
% implement a bandpass filter for accelerometer norm
if fcnt== 1,
    accl_flt_nrm(1) = 0.1535*nrm_f;
    nrm_f_m1 = nrm_f;
    fcnt = 2;
elseif fcnt== 2,
    accl_flt_nrm(2) = accl_flt_nrm(1);
    accl_flt_nrm(1) = 1.681*accl_flt_nrm(1)+ ...
        0.1535*nrm_f;
    nrm_f_m2 = nrm_f_m1;
    nrm_f_m1 = nrm_f;
    fcnt = 3;
else
    tmp = accl_flt_nrm(2);
    accl_flt_nrm(2) = accl_flt_nrm(1);
    accl_flt_nrm(1) = 1.681*accl_flt_nrm(1)-0.6839*tmp + ...
        0.1535*nrm_f - 0.1535*nrm_f_m2;
    nrm_f_m2 = nrm_f_m1;
    nrm_f_m1 = nrm_f;
end 
stationary = 0;
if and(abs(nrm_f)<0.1,and(gyro_nrm<1*pi/180,accl_flt_nrm(1)<0.02))
    stationary = 1;
elseif t>T_start,
    stationary = 0;
end

if att_init == 0,
    % wait until T_start
    % then average sensors until motion is detected
    if and(t>T_start,stationary),        % initial averaging while stationary
        int_T = int_T + T;
        int_w = int_w + w_ip_p*T;
        int_g = int_g + -f_ip_p*T;
        E(1) = atan2(int_g(2),int_g(3));
        E(2) = atan2(-int_g(1)/int_T,norm(int_g(2:3)/int_T));
        E(3) = 0;               % temporarily
        Rt2b = C_Rt2b(E);
        sys(1:4,1)    = Rot2Quat(Rt2b);
        sys(5:7,1)    = int_w/(int_T);         % x_g
        sys(8:10,1)   = [0;0;0];         % x_a
        Pp(1:3,1:3) = (diag([sigma_nu_a/ge,sigma_nu_a/ge,sqrt(Rm)])^2)/int_T;   %(1*pi/180)^2*eye(3,3); % rho uncertainty, rad
        Pp(4:6,4:6) = (sigma_nu_g)^2*eye(3,3)/int_T; % gyro bias uncertainty, rad/s, 4d/hr=0.001 d/s
        Pp(7:9,7:9) = 0*(2e-3)^2*eye(3,3); % accel bias uncertainty, m/s/s
        Pp(1,8) =  sqrt( Pp(1,1)*Pp(8,8) );
        Pp(2,9) =  sqrt( Pp(2,2)*Pp(9,9) );
        Pp(7:9,1:3) = Pp(1:3,7:9)';
        %diag(Pp)'
        %P_ant_symmetric = (Pp-Pp')
        Pm = Pp;
    elseif t>T_start,           % stop averaging as no longer stationary
        sys = x;         % x_a
        att_init = 1;
%        sys(5:7,1)    = 0*int_w/(int_T);         % x_g

        sprintf('Initialization ended due to motion detection at %0.5g',t)
    else
        sys = x;         % x_a
    end
    ha  = Ha(1:2,:);
    P_r1 = ha(1,:)*Pm*ha(1,:)'+Ra;
    P_r2 = ha(2,:)*Pm*ha(2,:)'+Ra;
else                            % state propagation while in motion
%     nrm_acc_avg = a_filt*nrm_acc_avg + (1-a_filt)*nrm_f;
%     % Because there is no pos and vel info available account for the 
%     % error in the nav frame rotation rate (at least approximately)
%     w_in_p =0*( (nrm_f-nrm_acc_avg)*T/Re + w_ie)^2; % wrong line
    
    w_ip_p = w_ip_p - x_g;      % correct gyro for bias, eqn. 11.25
    w_bn_b =  -w_ip_p;   % eqn. 11.24
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The following section intergrates the state vector
    W   = w_bn_b*T/2;
    w   = norm(W);
    if abs(w)>1
        w = w/abs(w);
        W,w,w_bn_b,T
        error('Integrated angle too large');
    end
    Wc = [0    -W(3) W(2)
        W(3)  0    -W(1)
        -W(2)  W(1)  0];
    W_mat = [0 -W'
        W Wc];
        if w == 0,
            sinwow = 1;
        else
            sinwow = sin(w)/w;
        end
        sys(1:4,1)    = (cos(w)*eye(4,4) + W_mat*sinwow)*b; % eqn. D.36 
        sys(5:7,1)    = x(5:7);         % x_g
        sys(8:10,1)   = x(8:10);        % x_a
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % The following section accumulates phi and Qd between measurements
        sigma_in = (norm(f_ip_p) - 9.8)*0;
        Z = zeros(3,3);
        I = eye(3,3);
        % compute rotation matrix half way through interval (predictor-corrector)
        Rn2b    =(quat2R(b)+quat2R(sys(1:4,1)))/2; 
        Rb2n    = Rn2b';
        F = [Z -Rb2n  Z
            Z  Fg    Z
            Z  Z     Fa];
            G = [I  Z    -Rb2n Z
                Z  I     Z     Z
                Z  Z     Z     I];
        Q = [I*sigma_in^2 Z          Z           Z
             Z        I*sigma_xg^2   Z           Z
             Z            Z     I*sigma_nu_g^2   Z
             Z            Z          Z       I*sigma_xa^2*0];
        [phi,q]=calc_Qd_phi(F,G*Q*G',T);
         Phi = phi*Phi;          % accumulate phi, see CH7
         Qd  = phi*(q+Qd)*phi';  % accumulate Qd,  see CH7
         % update cov at high rate only to make nice plots
         Pm  = Phi*Pm*Phi' + Qd;             % time propagate cov
         % prepare for next period of integration
         Phi = eye(9,9);       % reset Phi
         Qd  = zeros(9,9);     % reset Qd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following section incorporates KF measurements
         if t>T_start,
           nrm_accel = nrm_f;      
           if and(t>T_a+0.1,stationary),   % max 1 Hz update
               T_a = t;
               sprintf('KF Attitude correction at %0.5g',t)
                 Ha(:,7:9) = Rb2n;
                 ha  = Ha(1:2,:);
                 hatg_n =  Rb2n*(x_a-f_ip_p); % accel reading
                 Raa = (Ra + (10*nrm_accel)^2)*eye(2,2);    % account for acceleration as error
                 K   = Pm*ha'*inv(ha*Pm*ha'+Raa);     % KF gain
                 Pp  = Pm - K*(ha*Pm);               % meas update for cov
                 res = ya - hatg_n;
                 del_x = K*res(1:2);              % compute state correction
                 rho = del_x(1:3);
                 rho_cross = [0     -rho(3) rho(2)
                              rho(3) 0     -rho(1)
                             -rho(2) rho(1) 0];
                 Rn2b = Rn2b*(I-rho_cross);
                 sys(1:4)  = Rot2Quat(Rn2b);            % correct quaternion
                 sys(5:10) = sys(5:10) + del_x(4:9);    % correct state
                 del_x = 0*del_x;           % reset state correction
                 Pm = Pp;                   % get ready for time updates
           end
         end
          Ha(:,7:9) = Rb2n;
     ha  = Ha(1:2,:);
     P_r1 = ha(1,:)*Pm*ha(1,:)'+Ra;
     P_r2 = ha(2,:)*Pm*ha(2,:)'+Ra;
end  % else
     


%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function [sys]=mdlOutputs(t,x,u,ind_accel,Pr1,Pr2)
persistent d2r grav bias_a bias_g pos_o
persistent ym Rm sigma_nu_a ya Ra

if t<=0.005,
    d2r         = pi/180;
    grav        = [0;0;9.78];
    bias_a      = 0.01*randn(3,1);
    bias_g      = 0.0005*randn(3,1);
    
    ym = [1,0,0]';
    Rm = (1*pi/180)^2;      % rad^2, equiv to 1 deg, magnetometer noise
    
    sigma_nu_a = 2.2e-2;     % m/s/s/rt_Hz, velocity drift rate
    ge = 9.78;              % gravity magnitude
    ya = [0,0,ge]';
    Ra = (sigma_nu_a)^2*eye(2,2);      % rad^2, equiv to 1 deg, magnetometer noise
end

b   = x(1:4);      % rad
x_g = x(5:7);      % rad/s, 
x_a = x(8:10);     % m/s/s

w_ip_p = u(1:3);   % rad/s,   body frame angular rate
f_ip_p = -u(4:6);  % m/s/s,   body frame specific force

[E] = b2Euler(b);  % Euler angles

Rn2b    = quat2R(b); 
Rb2n    = Rn2b';

hat_gb =  (x_a-f_ip_p);
hat_gn =  Rb2n*hat_gb;
res_a = ya - hat_gn;
roll      = atan2(hat_gb(2),hat_gb(3));
ptch      = atan2(-hat_gb(1),norm(hat_gb(2:3)));

% format output
sys = [x(1:7);ind_accel;Pr1;Pr2;E(1);roll;E(2);ptch;res_a(1:2);sqrt(diag(Ra))];
if length(sys)~=18,
    sys
    x
    E
    diag(Ra)
    res_a
end

%
function sys=mdlGetTimeOfNextVarHit(t,x,u,p)
sys = [];




%
%=============================================================================
% mdlTerminate
% Perform any end of simulation tasks.
%=============================================================================
%
function sys=mdlTerminate(t,x,u)

%RGR
sys = [];
% end mdlTerminate



function [ans]=limit_pi(x)
two_pi=2*pi;
ans = x;
for i=1:length(x);
    while ans(i)>pi
        ans(i) = ans(i) - two_pi;
    end
    while ans(i)<-pi
        ans(i) = ans(i) + two_pi;
    end
end



% convert a quternion b to Euler angles
function [E] = b2Euler(b)
E(2,1) = asin( -2*(b(2)*b(4) + b(1)*b(3)) );
E(1,1) = atan2( 2*(b(3)*b(4)-b(1)*b(2)) , 1-2*(b(2)^2+b(3)^2) );
E(3,1) = atan2( 2*(b(3)*b(2)-b(1)*b(4)) , 1-2*(b(3)^2+b(4)^2) );



% convert a quaternion b to a rotation matrix
function [Rn2b]=quat2R(b)
if norm(b)~= 0,
    b = b/norm(b);
    B = b(1);
    Bv(:,1)= b(2:4);
    Bc = [0    -Bv(3) Bv(2)
         Bv(3)  0    -Bv(1)
        -Bv(2)  Bv(1) 0];
    Rn2b = (B*B-Bv'*Bv)*eye(3,3)+2*Bv*Bv'+2*B*Bc;
%     Rn2b = [(b(1)^2+b(2)^2-b(3)^2-b(4)^2) 2*(b(2)*b(3)-b(1)*b(4)) 2*(b(1)*b(3)+b(2)*b(4))
%             2*(b(2)*b(3)+b(1)*b(4)) (b(1)^2-b(2)^2+b(3)^2-b(4)^2) 2*(b(3)*b(4)-b(1)*b(2))
%             2*(b(2)*b(4)-b(1)*b(3)) 2*(b(1)*b(2)+b(3)*b(4))   b(1)^2-b(2)^2-b(3)^2+b(4)^2]
else
    Rn2b = eye(3,3)    % fault condition
    error('Norm b = 0');
end
% b_check1=Rot2Quat(Rn2b)
% n1 = norm(b_check1)


function [b]=Rot2Quat(R)
[U,S,V]=svd(R);
R = U*V';
if 1+R(1,1)+R(2,2)+R(3,3) > 0
    b(1,1)    = 0.5*sqrt(1+R(1,1)+R(2,2)+R(3,3));
    b(2,1)    = (R(3,2)-R(2,3))/4/b(1);
    b(3,1)    = (R(1,3)-R(3,1))/4/b(1);
    b(4,1)    = (R(2,1)-R(1,2))/4/b(1);
    b       = b/norm(b);    % renormalize
else
    R
    error('R diagonnal too negative.')
    b = zeros(4,1);
end


% convert Euler angles to a rotation matrix
function [Rt2b] = C_Rt2b(x)
c_r = cos(x(1));
s_r = sin(x(1));
c_p = cos(x(2));
s_p = sin(x(2));
c_y = cos(x(3));
s_y = sin(x(3));
Rt2b =[ c_y*c_p             s_y*c_p                 -s_p
    (-s_y*c_r+c_y*s_p*s_r) (c_y*c_r+s_y*s_p*s_r)   (c_p*s_r)
    ( s_y*s_r+c_y*s_p*c_r) (-c_y*s_r+s_y*s_p*c_r)  (c_p*c_r)];
% b_check2=Rot2Quat(Rt2b)
% n2 = norm(b_check2)



function [Re2t] = C_Re2t(lat,lng)
c_lat = cos(lat);
s_lat = sin(lat);
c_lng = cos(lng);
s_lng = sin(lng);
Re2t  =[ -s_lat*c_lng        -s_lat*s_lng       c_lat
         -s_lng             c_lng              0 
        -c_lat*c_lng        -c_lat*s_lng      -s_lat];



function [Ob2Euler] = C_Oteb(x)
c_r = cos(x(1));
s_r = sin(x(1));
c_p = cos(x(2));
t_p = tan(x(2));
Ob2Euler =[ 1             s_r*t_p     c_r*t_p
            0             c_r         -s_r
            0             s_r/c_p     c_r/c_p];


%
%	F the continuous time state transition matrix
%	Q the continuous time process noise covariance matrix
%	Ts the time step duration
%
%	phi the discrete time transition matrix
%	Qd the equivalent discrete time driving noise
%
function [phi,Qd]=calc_Qd_phi(F,Q,Ts)
[n,m] = size(F);
if n~=m,
	error('In calc_Qd_phi, the F matrix must be square');
end	%if

chi = [ -F      Q
		 0*eye(n) F']*Ts;
gamma=expm(chi);

phi = gamma((n+1):(2*n),(n+1):(2*n))';
Qd  = phi*gamma(1:n,(n+1):(2*n));