function M = galaxy_model()
%
% Author : nicolas.douillet9 (at) gmail.com, 2008-2024.


% Galaxy parameter meaning 
% nb_stars = 5e4;  % number of stars
% amplitude = 1;   % thickness of the galactic disk
% radius = 1;      % basis radius
% curvature = 9;   % the higher this value is, the greater the spiral is wrapped around itself
% nb_arms = 2;     % number of arms, integer
% way = -1;        % -1 / 1 : rotation way
% R0 = 3;          % radius translation parameter / zoom /unzoom in the spiral centre
% mean_val = 3;    % influence stars density at the centre
% sdn = 4;         % standard deviations number
% k = 1;           % phase at the origin
% a = 1;           % ellipsoidal X axis
% b = 1;           % ellipsoidal Y axis


% Parameter set value examples

default_param_set = [5e4, 1, 1, 9, 2, -1, 3, 3, 4, 1, 1, 1, 0.08];   % spiral galaxy

param_set2 =        [5e4, 3, 1, 360, 2, -1, 4, 3, 4, 3, 2, 1, 0.08]; % old metalic galaxy


[nb_stars, amplitude, radius, curvature, nb_arms, way, R0, mean_val, sdn, k, a, b] = struct('v', num2cell(default_param_set)).v;


% Noise parameters
noise = 0.08;
X_noise_percent = noise;
Y_noise_percent = noise;
Z_noise_percent = noise;

% Display parameters
profile_view = true;
starsize = 'adaptated'; % 'constant'
az = 0;  % -40; % azimut
el = 90; % 78; % elevation
zoom_lvl = 1.5;

% Data distribution
R = radius*(randn(1,nb_stars)) + mean_val;
% Pass band filtering on [mean_val-sdn*radius; mean_val+sdn*radius] where
% radius is the standard deviation of the distribution

e = find(abs(R)-mean_val < sdn*radius);
s = numel(e);
R = R(e);
Theta = 2*pi*randn(1,nb_stars);
Theta = Theta(e);
X = a*R.*cos(Theta);
Y = b*R.*sin(Theta);
Phi = angle(X+Y*1i);
% Spiral shape handle function
%
% f = @(A,s,c,r,m,angl,p) A*exp(-(r.^2)/sdn^3).*sin(s*c*log(r)+m*angl+p*pi); 
f = @(A,s,c,r,m,angl,p) A*exp(-(((r-mean_val).^2)/2/radius^2)/sdn^3).*sin(s*c*log(r)+m*angl+p*pi)/radius/sqrt(2*pi); 
Z = f(amplitude,way,curvature,R+R0,nb_arms,Phi,k);
% Normal noise
%
X_noise = X_noise_percent*amplitude*randn(1,s);
Y_noise = Y_noise_percent*amplitude*randn(1,s);
Z_noise = Z_noise_percent*amplitude*randn(1,s);
X = X + X_noise;
Y = Y + Y_noise;
Z = Z + Z_noise;
M = cat(2,X',Y',Z');

% Display
h = figure;
set(h,'Position',get(0,'ScreenSize'));
set(gcf,'Color',[0 0 0]);
% Color parameters
%
colormap('bone'); % 'bone', 'hot', 'autumn', 'jet', etc.
cmap = colormap;
% % RGB color filters
%
% cmap(1:end,1) = 0; % R
% cmap(1:end,2) = 0; % G
% cmap(1:end,3) = 0; % B
%
% cmap = circshift(cmap,1,2); % color shift
if profile_view
    
    subplot(121);
    
end

ax = gca;

if strcmpi(starsize,'constant')
    
    scatter3(X,Y,Z,1,abs(Z),'.'), hold on;
    
elseif strcmpi(starsize,'adaptated')
    
    scatter3(X,Y,Z,abs(mean_val+sdn-R).^2,abs(Z),'*'), hold on;
    
end

axis equal, axis tight, axis off;
set(ax,'Color',[0 0 0]);
colormap(cmap);
ax.Clipping = 'off';
view(az,el);
zoom(zoom_lvl);

if profile_view
    
    subplot(122);
    ax = gca;
    
    if strcmpi(starsize,'constant')
        
        scatter3(X,Y,Z,1,abs(Z),'.'), hold on;
        
    elseif strcmpi(starsize,'adaptated')
        
        scatter3(X,Y,Z,abs(mean_val+sdn-R).^2,abs(Z),'*'), hold on;
        
    end
    
    axis equal, axis tight, axis off;
    set(ax,'Color',[0 0 0]);
    colormap(cmap);
    ax.Clipping = 'off';
    view(0,0);
    zoom(zoom_lvl);
    
end


end % galaxy_model