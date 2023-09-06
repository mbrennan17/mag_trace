function [maglinen, maglines] = mag_trace(x_rj,y_rj,z_rj,int_field_function,ext_field_function,varargin)
%INPUTS:
%     x_rj: Position x-coordinates (n) (Rj) in cartesian frame required by field functions (ie IAU_JUPITER, J2000, etc) from which magnetic field line is traced, n position vectors.
%     y_rj: Position y-coordinates (n) (Rj) in cartesian frame required by field functions (ie IAU_JUPITER, J2000, etc) from which magnetic field line is traced, n position vectors.
%     z_rj: Position z-coordinates (n) (Rj) in cartesian frame required by field functions (ie IAU_JUPITER, J2000, etc) from which magnetic field line is traced, n position vectors.
%     int_field_function: internal field function handle, see example below
%     ext_field_function: external field function handle, see example below
%
%     Optional Inputs:
%       Must provide a structure with the following input variables, otherwise default values will be used:
%           user_params.alt: Altitude [1] (km) of jupiter termination of field lines
%           user_params.body_radii: Body radii [1,3] (km), ie Jupiter as [71492 71492 66854] for ellipsoid or [71492 71492 71492] for sphere;
%           user_params.xprop:  Distance [1] (Radii) of total propogation, (ie Rj for Jupiter)
%           user_params.rmax:  Max radius [1] (Radii) of field line from central body (suggested 200 Rj for Jupiter)
%           user_params.ext_rmin: Min radius [1] (Radii) for including external/current-sheet component in field line computation
%               (ext_rmin = 0 for always use external/current sheet, ext_rmin = rmax+1 for never use external/current sheet)
%
%OUTPUTS:
%       maglinen: {n}(m1, 3) (Radii) northern magnetic field line (of m1 pts) for each n position in cartesian frame
%       maglines: {n}(m2, 3) (Radii) southern magnetic field line (of m2 pts) for each n position in cartesian frame
%          Note m1 and m2 points are arbitrary/uncontrolled results from the field line tracing
%
%       User can find magnetic footprints and/or terminal points of each magline trace, see following example:
%             for n=1:length(maglinen)
%                 magfpn0(n,:) = maglinen{n}(end,:);
%                 magfps0(n,:) = maglines{n}(end,:);
%             end
%
%
%EXAMPLE 1:
% pos = [0,5,0; 0 6 0; 0 7 0; 0 8 0; 0 199 0];
% int_field_function = @(x,y,z) (jovian_jrm33_order13_internal_xyz(x,y,z));
% ext_field_function = @(x,y,z) (con2020_model_xyz('analytic',x,y,z));
% [maglinen, maglines] = mag_trace(pos(:,1),pos(:,2),pos(:,3), int_field_function, ext_field_function);
%
%EXAMPLE 2:
% pos = [0,5,0; 0 6 0; 0 7 0; 0 8 0; 0 199 0];
% ext_params = con2020_model_xyz('default_values'); % get defaults
% ext_params.xt__cs_tilt_degs = 0; % change tilt value
% int_field_function = @(x,y,z) (jovian_jrm33_order13_internal_xyz(x,y,z));
% ext_field_function = @(x,y,z) (con2020_model_xyz('hybrid',x,y,z,ext_params));
% input.alt         =  200; % (km)
% input.body_radii  = [71492 71492 71492]; % (km) % size 1x3
% input.xprop       =  500; % (Rj)
% input.rmax        =  100; % (Rj)
% input.ext_rmin    =    0; % (Rj)
% input.error_check =    1;
% [maglinen, maglines]  = mag_trace(pos(:,1),pos(:,2),pos(:,3), int_field_function, ext_field_function);
%
% Created by Martin Brennan along with Chris Lawler and Rob Wilson,  July 2023

switch numel(varargin)
    case 0
        Defaults.alt         =  400; % (km)
        Defaults.body_radii  = [71492 71492 66854]; % (km) % size 1x3
        Defaults.xprop       = 1000; % (Rj)
        Defaults.rmax        =  200; % (Rj)
        Defaults.ext_rmin    =    0; % (Rj)
        Defaults.error_check =    1;
        
        alt         = Defaults.alt ;
        body_radii  = Defaults.body_radii ;
        xprop       = Defaults.xprop      ;
        rmax        = Defaults.rmax       ;
        ext_rmin    = Defaults.ext_rmin   ;
        error_check = Defaults.error_check;
    case 1
        user_params = varargin{1};
        if ~isstruct(user_params)
            error('Must be a structure of terms to use in code')
        end
        alt         = double(user_params.alt        );
        body_radii  = double(user_params.body_radii );
        xprop       = double(user_params.xprop      );
        rmax        = double(user_params.rmax       );
        ext_rmin    = double(user_params.ext_rmin   );
        error_check = double(user_params.error_check);
    otherwise
        error('ERROR: code accepts either 3 or 4 input arguments, not %d.',numel(varargin)+3)
end


if error_check
    
    % Check position input
    if (~isnumeric(x_rj)) || (size(x_rj,2) ~= 1), error('ERROR: x-coordiante variable (x_rj) must be an n x 1 (n rows, 1 columns) array of numbers'); end
    if (~isnumeric(y_rj)) || (size(y_rj,2) ~= 1), error('ERROR: y-coordiante variable (y_rj) must be an n x 1 (n rows, 1 columns) array of numbers'); end
    if (~isnumeric(z_rj)) || (size(z_rj,2) ~= 1), error('ERROR: z-coordiante variable (z_rj) must be an n x 1 (n rows, 1 columns) array of numbers'); end
    
    if (length(x_rj)~=length(y_rj)) && (length(x_rj)~=length(z_rj)), error('ERROR: postion coordinates (x_rj,y_rj,z_rj) must be the same size (n elements)'); end
    
    % Check if internal and external field models exist
    if ~(exist('int_field_function','var'))      , error('ERROR: Internal field function not found'); end
    if ~(exist('ext_field_function','var'))      , error('ERROR: External field function not found'); end
    if ~isa(int_field_function,'function_handle'), error('ERROR: Internal field function is not a function handle'); end
    if ~isa(ext_field_function,'function_handle'), error('ERROR: External field function is not a function handle'); end
    
    
    % Check optional inputs
    if  (~isnumeric(alt       )) || (numel(alt       )~=1) || (alt<-max(body_radii)), error('ERROR: Altitude optional input (user_params.alt) must be a scalar double number, typically >0, but must be > -Rj'                                ); end
    if  (~isnumeric(body_radii)) || (numel(body_radii)~=3) || (any(body_radii<=0)  ), error('ERROR: Body radii optional input (user_params.body_radii) must be a 1x3 vector (ie [71492 71492 66854] or [71492 71492 71492])'                  ); end
    if  (~isnumeric(xprop     )) || (numel(xprop     )~=1) || (xprop         <=0   ), error('ERROR: Total propogated distance optional input (user_params.xprop) must be a non-zero positive scalar double number in Rj'                      ); end
    if  (~isnumeric(rmax      )) || (numel(rmax      )~=1) || (rmax          <=0   ), error('ERROR: Max field line radius optional input (user_params.rmax) must be a non-zero positive scalar double number in Rj'                           ); end
    if  (~isnumeric(ext_rmin  )) || (numel(ext_rmin  )~=1) || (ext_rmin      < 0   ), error('ERROR: Min radius for including external/current-sheet model optional input (user_params.ext_rmin) must be a positive scalar double number in Rj'); end
    
end

npos = length(x_rj); % number of positions
rmax_sq     = rmax     * rmax    ;
ext_rmin_sq = ext_rmin * ext_rmin;

% Manually Adjusted Solver options
alt_radii = (body_radii + alt) / body_radii(1);
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8,'Vectorized','on','Events', @(t, y) terminate_at_body(y, alt_radii)); % 1e-8 AbsTol = ~1m.

% Initialize output field lines
maglinen = cell(npos,1);
maglines = maglinen;

% Loop through n positions and propogate along both field line directions.
pos_in = [x_rj, y_rj, z_rj];
for n = 1:npos
    [~, pos_out1] = ode45(@(s,pos_out) mag_field_evaluate(pos_out, int_field_function, ext_field_function, rmax_sq, ext_rmin_sq), [0 -xprop], pos_in(n,:), opts);
    [~, pos_out2] = ode45(@(s,pos_out) mag_field_evaluate(pos_out, int_field_function, ext_field_function, rmax_sq, ext_rmin_sq), [0  xprop], pos_in(n,:), opts);
    
    maglinen{n} = pos_out1(n,:);
    maglines{n} = pos_out2(n,:);
end

end
%% Mag Field Evaluation
function [Bout] = mag_field_evaluate(rin, int_field_function, ext_field_function, rmax_sq, ext_rmin_sq)
% Checks if integrator trying to propogate past defined max radii (rmax)
R_sq = rin(1)*rin(1)+rin(2)*rin(2)+rin(3)*rin(3);

if R_sq <= rmax_sq
    % Internal Field: compute magnetic direction B-vector
    Bout = int_field_function(rin(1),rin(2),rin(3)); % single 1x3 vector
    
    % External Field included if > rmin: compute magnetic direction
    if R_sq >= ext_rmin_sq
        Bout = Bout + ext_field_function(rin(1),rin(2),rin(3)); % single 1x3 vector
    end
    % Convert B-vector to 3x1 and normalize
    Bout = [Bout(1);Bout(2);Bout(3)] / sqrt(Bout(1)*Bout(1)+Bout(2)*Bout(2)+Bout(3)*Bout(3));
else
    Bout = [0; 0; 0];
end
end

%% Termination at Body
function [value, isterminal, direction] = terminate_at_body(r, body_radii)
% Checks if propogation has breached body surface/threshold
% This will terminate if r = [0, 0, 0] too, as that's in the elipse
% Ellipsoid equation: inequality checks inner or outer space
value = ((r(1)/body_radii(1))^2 + (r(2)/body_radii(2))^2 + (r(3)/body_radii(3))^2) < 1;
isterminal = 1; % Needed for odezero
direction  = 0; % Needed for odezero
end