function [] = galaxy_model() % M


% Spiral galaxy parameter set example :
%
% [nb_stars = 5e4, amplitude = 0.5, radius = 1, curvature = 9, nb_arms = 2, way = -1, R0 = 3, mean_val = 2, sdn = 4, 
% X_noise_percent = Y_noise_percent = Z_noise_percent = 0.08, a = b = 1, cmap = colormap('bone')]

% Elliptic galaxy parameter set example #1 :
%
% [nb_stars = 1e5, amplitude = 3, radius = 2, curvature = 15, nb_arms = 2,
% way = -1, R0 = 4, mean_val = 1.8, sdn = 1, a = 2, b = 1, noise = 0.15,
% cmap = colormap('autumn')]

% Elliptic galaxy parameter set example #2 :
%
% [nb_stars = 5e4, amplitude = 3, radius = 1, curvature = 360, nb_arms = 2,
% way = -1, R0 = 4, mean_val = 3, sdn = 3, a = 2, b = 1, noise = 0.08,
% colormap('autumn')]


% Spiral geometric parameters
%
nb_stars = 5e4;  % number of stars
amplitude = 0.5; % thickness of the galactic disk
radius = 1;      % basis radius
curvature = 9;   % the higher this value is, the greater the spiral is wrapped around itself
nb_arms = 1;     % number of arms, integer
way = -1;        % -1 / 1 : rotation way
R0 = 3;          % radius translation parameter
mean_val = 2;    % influence stars density at the centre
sdn = 4;         % standard deviations number
k = 1;           % phase at the origin

a = 1;           % ellipsoidal X axis
b = 1;           % ellipsoidal Y axis

% Noise parameters
%
X_noise_percent = 0.08;
Y_noise_percent = 0.08;
Z_noise_percent = 0.08;

% Display parameters
%
profile_view = true;
starsize = 'adaptated'; % 'constant'
az = 0; % -40; % azimut
el = 90; % 78; % elevation
zoom_lvl = 1.5;

% Data distribution
R = radius*(randn(1,nb_stars)) + mean_val;

% Pass band filtering on [mean_val-sdn*radius; mean_val+sdn*radius] where
% radius is the standard deviation of the distribution
%
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

% M = cat(2,X',Y',Z');

h = figure;
set(h,'Position',get(0,'ScreenSize'));
set(gcf,'Color',[0 0 0]);

% Color parameters
%
colormap('bone'); % 'bone', 'hot', 'autumn', etc.
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