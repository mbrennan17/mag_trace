function [maglinen, maglines] = mag_trace(pos,int_field_function,ext_field_function,varargin)
%INPUTS:
%     pos: Positions (n, 3) (Rj) in IAU Frame (ie IAU_JUPITER) from which magnetic field line is traced, n position vectors.
%     int_field_function: internal field function handle, see example below
%     ext_field_function: external field function handle, see example below
% 
%     Optional Inputs:
%       Must provide a structure with the following input variables, otherwise default values will be used:
%           user_params.alt: Altitude [1] (km) of jupiter termination of field lines 
%           user_params.body_radii: Body radii [1,3] (km), ie Jupiter as [71492 71492 66854] for ellipsoid or [71492 71492 71492] for sphere;
%           user_params.xprop:  Distance [1] (Radii) of total propogation, (ie Rj for Jupiter)
%           user_params.rmax:  Max radius [1] (Radii) of field line from central body (suggested 200 Rj for Jupiter)
%           user_params.ext_rmin: Min radius [1] (Radii) for including external/current-sheet component in field line computation (ext_rmin = 0 for always use external/current sheet)
       
%OUTPUTS:
%       maglinen: {n}(m1, 3) (Radii) northern magnetic field line (of m1 pts) for each n position in IAU frame
%       maglines: {n}(m2, 3) (Radii) southern magnetic field line (of m2 pts) for each n position in IAU frame
%          Note m1 and m2 points are arbitrary/uncontrolled results from the field line tracing
%       
%       User can find magnetic footprints and/or terminal points of each magline trace, see following example:
%             for n=1:length(maglinen)
%                 magfpn0(n,:) = maglinen{n}(end,:);
%                 magfps0(n,:) = maglines{n}(end,:);
%             end


%EXAMPLE 1:
% pos = [0,5,0; 0 6 0; 0 7 0; 0 8 0; 0 199 0];
% int_field_function = @(x,y,z) (jovian_jrm33_order13_internal_xyz(x,y,z));
% ext_field_function = @(x,y,z) (con2020_model_xyz('analytic',x,y,z));
% [maglinen, maglines] = mag_trace_jup(pos, int_field_function, ext_field_function);

%EXAMPLE 2:
% pos = [0,5,0; 0 6 0; 0 7 0; 0 8 0; 0 199 0];
% ext_params = con2020_model_xyz('default_values'); % get defaults
% ext_params.xt__cs_tilt_degs = 0; % change tilt value
% int_field_function = @(x,y,z) (jovian_jrm33_order13_internal_xyz(x,y,z));
% ext_field_function = @(x,y,z) (con2020_model_xyz('hybrid',x,y,z,ext_params));
% input.alt       = 200; % (km)
% input.body_radii = [71492 71492 71492]; % (km) % size 1x3
% input.xprop     = 500;% (Rj)
% input.rmax      =  100; % (Rj)
% input.ext_rmin  =    0; % (Rj)
% input.error_check = 1;
% [maglinen, maglines] = mag_trace_jup(pos, int_field_function, ext_field_function);

% Created by Martin Brennan along with Chris Lawler and Rob Wilson,  July 2023


Defaults.alt       = 400; % (km)
Defaults.body_radii = [71492 71492 66854]; % (km) % size 1x3
Defaults.xprop     = 1000;% (Rj)
Defaults.rmax      =  200; % (Rj)
Defaults.ext_rmin  =    0; % (Rj)
Defaults.error_check = 1;

switch numel(varargin)
    case 1
        user_params = varargin{1};
        if ~isstruct(user_params)
            error('Must be a structure of terms to use in code')
        end
        alt         = double(user_params.alt);
        body_radii   = double(user_params.body_radii);
        xprop       = double(user_params.xprop);
        rmax        = double(user_params.rmax);
        ext_rmin    = double(user_params.ext_rmin);
        error_check = double(user_params.error_check);
    case 0

        alt = Defaults.alt ;
        body_radii = Defaults.body_radii;
        xprop =  Defaults.xprop;
        rmax =  Defaults.rmax;
        ext_rmin =  Defaults.ext_rmin;
        error_check = Defaults.error_check;
    otherwise
         error('ERROR: code accepts either 3 or 4 input arguments, not %d.',numel(varargin)+3)
end


if error_check

    % Check position input
    if (~isnumeric(pos)) || (size(pos,2) ~= 3), error('ERROR: Position variable (pos) must be an n x 3 (n rows, 3 columns) array of numbers'); end
    
    % Check if internal and external field models exist
    if (isempty(int_field_function)), error('ERROR: Internal field function not found'); end
    if (isempty(ext_field_function)), error('ERROR: External field function not found'); end

    % Check optional inputs
    if  (~isnumeric(alt      )) || (numel(alt      )~=1) || alt<-max(body_radii), error('ERROR: Altitude optional input (user_params.alt) must be a scalar double number, typically >0, but must be > -Rj');end
    if  (~isnumeric(body_radii)) || (numel(body_radii)~=3) || any(body_radii<=0), error('ERROR: Body radii optional input (user_params.body_radii) must be a 1x3 vector (ie [71492 71492 66854] or [71492 71492 71492])');end
    if  (~isnumeric(xprop)) || (numel(xprop)~=1) || xprop==0, error('ERROR: Total propogated distance optional input (user_params.xprop) must be a non-zero scalar double number in Rj');end
    if  (~isnumeric(rmax)) || (numel(rmax)~=1) || rmax<=0, error('ERROR: Max field line radius optional input (user_params.rmax) must be a non-zero positive scalar double number in Rj');end
    if  (~isnumeric(ext_rmin)) || (numel(ext_rmin)~=1) || ext_rmin<0, error('ERROR: Min radius for including external/current-sheet model optional input (user_params.ext_rmin) must be a positive scalar double number in Rj');end

end


rmax_sq     = rmax * rmax   ;
ext_rmin_sq = ext_rmin * ext_rmin;

% Manually Adjusted Solver options
alt_radii = (body_radii + alt) / body_radii(1);
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8,'Vectorized','on','Events', @(t, y) terminate_at_body(y, alt_radii)); % 1e-8 AbsTol = ~1m. 

% Initialize output field lines
maglinen = cell(size(pos,1),1);
maglines = maglinen;

% Loop through n positions and propogate along both field line directions.
for n = 1:size(pos,1)
    [~, yout1] = ode45(@(x,y) mag_field_evaluate(y, int_field_function, ext_field_function, rmax_sq, ext_rmin_sq), [0 -xprop], pos(n,:), opts);
    [~, yout2] = ode45(@(x,y) mag_field_evaluate(y, int_field_function, ext_field_function, rmax_sq, ext_rmin_sq), [0  xprop], pos(n,:), opts);

        maglinen{n} = yout1;
        maglines{n} = yout2;
end

end

function [rout] = mag_field_evaluate(rin, int_field_function, ext_field_function, rmax_sq, ext_rmin_sq)
    % Checks if integrator trying to propogate past defined max radii (rmax)
    R_sq = rin(1)*rin(1)+rin(2)*rin(2)+rin(3)*rin(3);

    if R_sq < rmax_sq
        
        % Internal Field
        rout = int_field_function(rin(1),rin(2),rin(3)); % single 1x3 vector

        % External Field included if > rmin
        if R_sq >= ext_rmin_sq 
                rout = rout + ext_field_function(rin(1),rin(2),rin(3)); % single 1x3 vector
        end
        % Convert to 3x1 and normalize
        rout = [rout(1);rout(2);rout(3)] / sqrt(rout(1)*rout(1)+rout(2)*rout(2)+rout(3)*rout(3));
    else
        rout = [0; 0; 0];
    end
end

function [value,isterminal, direction] = terminate_at_body(r, body_radii)
    % Checks if propogation has breached body surface/threshold
    % This will terminate if r = [0, 0, 0] too, as that's in the elipse
    % Ellipsoid equation: inequality checks inner or outer space
    value = ((r(1)/body_radii(1))^2 + (r(2)/body_radii(2))^2 + (r(3)/body_radii(3))^2) < 1;
    isterminal = 1; % Needed for odezero 
    direction  = 0; % Needed for odezero
end