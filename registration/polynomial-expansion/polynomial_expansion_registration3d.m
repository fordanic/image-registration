function [update] = ...
    polynomial_expansion_registration3d(moving, fixed, varargin)
% POLYNOMIAL_EXPANSION_REGISTRATION3D Estimates a displacement field using polynomial expansion
%
% INPUT ARGUMENTS
% moving            - Moving image
% fixed             - Fixed image
% 
% OPTIONAL INPUT ARGUMENTS
% 'signalModel'             - Local signal model to use when computing the
%                             polynomial expansion transformation
%                             'linear' (deafult), 'quadratic'
%
% 'transformationModel'     - Transformation model for estimating the
%                             displacement field
%                             translation, affine, non-rigid (default)
%
% 'multiModal'              - Set wheteher to perform multi-modal or
%                             uni-modal image registration
%                             false (default), true
%
% 'numberOfChannels'        - Number of channels to use in when computing
%                             the entropy (based on channel coding). This
%                             is only relevant if multiModal is set to
%                             true.
%                             Default value is 8
%
% OUTPUT ARGUMENTS
% update
%   displacementUpdate      - Estimated update field
%   certaintyUpdate         - Certainty related to the estimated update field
%   transformationMatrix    - Estimate transformation matrix (only if 
%                             transformation model is set to translation or affine)


% Copyright (c) 2012 Daniel Forsberg
% danne.forsberg@outlook.com
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Setup default parameters
% linear, quadratic
signalModel = 'linear';

% translation, affine, non-rigid
transformationModel = 'non-rigid';
multiModal = false;

% Only valid for multi-modal registration
numberOfChannels = 16;

% Only valid for non-rigid registration
sigma = 1.5;
alpha = 0.01;

% Overwrites default parameter
for k=1:2:length(varargin)
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

%% Perform polynomial expansion
if strcmp(signalModel,'quadratic')
    [A_moving, b_moving, c_moving] = make_Abc_fast(moving);
    [A_fixed, b_fixed, c_fixed] = make_Abc_fast(fixed);
    
    A = (A_moving + A_fixed)/2;
    delta_b = b_fixed - b_moving;
else
    [b_moving, c_moving] = make_bc_fast(moving);
    [b_fixed, c_fixed] = make_bc_fast(fixed);
    if multiModal
        b = b_fixed;
        
        [delta_c mask] = estimate_delta_c(c_fixed,c_moving,numberOfChannels);
        delta_c(mask ~= 1) = 0;
        mask = repmat(mask,[1 1 1 3]);
        b(mask ~= 1) = 0;
    else
        b = (b_moving + b_fixed)/2;
        
        delta_c = c_fixed - c_moving;
    end 
end

if strcmp(signalModel,'quadratic')
    switch transformationModel
        case 'translation'
            [G, h] = build_G_h_quadratic3d(A, delta_b, transformationModel);
            d = G \ h;
            update.transformationMatrix = eye(4);
            update.transformationMatrix(1:3,4) = d;
        case {'rigid','affine'}
            [G, h] = build_G_h_quadratic3d(A, delta_b, transformationModel);
            p = G \ h;
            
            update.transformationMatrix = eye(4);
            update.transformationMatrix(1:3,1:4) = [1+p(1) p(2) p(3) p(10);...
                p(4) 1+p(5) p(6) p(11);...
                p(7) p(8) 1+p(9) p(12)];
        case 'non-rigid'
            A11 = A(:,:,:,5);
            A12 = 2*A(:,:,:,4);
            A13 = 2*A(:,:,:,8);
            A21 = 2*A(:,:,:,2);
            A22 = A(:,:,:,1);
            A23 = 2*A(:,:,:,7);
            A31 = 2*A(:,:,:,6);
            A32 = 2*A(:,:,:,3);
            A33 = A(:,:,:,9);
            B1 = delta_b(:,:,:,2);
            B2 = delta_b(:,:,:,1);
            B3 = delta_b(:,:,:,3);
            
            % Set the elements of the equation system
            G11 = A11.^2 + A21.^2 + A31.^2;
            G12 = A11.*A12 + A21.*A22 + A31.*A32;
            G13 = A11.*A13 + A21.*A23 + A31.*A33;
            G22 = A12.^2 + A22.^2 + A32.^2;
            G23 = A12.*A13 + A22.*A23 + A32.*A33;
            G33 = A13.^2 + A23.^2 + A33.^2;
            H1 = A11.*B1 + A21.*B2 + A31.*B3;
            H2 = A12.*B1 + A22.*B2 + A32.*B3;
            H3 = A13.*B1 + A23.*B2 + A33.*B3;
            
            % Smooth the elements of the equation system
            G11 = gauss_smoothing(G11, sigma);
            G12 = gauss_smoothing(G12, sigma);
            G13 = gauss_smoothing(G13, sigma);
            G22 = gauss_smoothing(G22, sigma);
            G23 = gauss_smoothing(G23, sigma);
            G33 = gauss_smoothing(G33, sigma);
            H1 = gauss_smoothing(H1, sigma);
            H2 = gauss_smoothing(H2, sigma);
            H3 = gauss_smoothing(H3, sigma);
            
            % Add alpha*identity matrix for stability
            scaleFactor = max([G11(:); G22(:); G33(:)]);
            G11 = G11 + alpha * scaleFactor;
            G22 = G22 + alpha * scaleFactor;
            G33 = G33 + alpha * scaleFactor;
            
            det = G33.*G12.^2 - 2*G12.*G13.*G23 + G22.*G13.^2 + G11.*G23.^2 - G11.*G22.*G33;
            
            % Estimate the displacement field
            update.displacement = cell(3,1);
            update.displacement{1} = removenan((H1.*(G23.^2 - G22.*G33) - ...
                H3.*(G12.*G23 - G13.*G22) - ...
                H2.*(G13.*G23 - G12.*G33))./(det + eps));
            update.displacement{2} = removenan((H2.*(G13.^2 - G11.*G33) - ...
                H3.*(G12.*G13 - G11.*G23) - ...
                H1.*(G13.*G23 - G12.*G33))./(det + eps));
            update.displacement{3} = removenan((H3.*(G12.^2 - G11.*G22) - ...
                H2.*(G12.*G13 - G11.*G23) - ...
                H1.*(G12.*G23 - G13.*G22))./(det + eps));
            
            % Estimate a certainty
            update.certainty = sqrt(b_fixed(:,:,:,2).^2 + b_fixed(:,:,:,1).^2 + b_fixed(:,:,:,3).^2);
    end
else
    switch transformationModel
        case 'translation'
            [G, h] = build_G_h_linear3d(b, delta_c, transformationModel);
            d = G \ h;
            update.transformationMatrix = eye(4);
            update.transformationMatrix(1:3,4) = d;
        case {'rigid','affine'}
            [G, h] = build_G_h_linear3d(b, delta_c, transformationModel);
            p = G \ h;
            
            update.transformationMatrix = eye(4);
            update.transformationMatrix(1:3,1:4) = [1+p(1) p(2) p(3) p(10);...
                p(4) 1+p(5) p(6) p(11);...
                p(7) p(8) 1+p(9) p(12)];
        case 'non-rigid'
            % Set the elements of the equation system
            % Set the elements of the equation system
            G11 = b(:,:,:,2) .* b(:,:,:,2);
            G12 = b(:,:,:,2) .* b(:,:,:,1);
            G13 = b(:,:,:,2) .* b(:,:,:,3);
            G22 = b(:,:,:,1) .* b(:,:,:,1);
            G23 = b(:,:,:,1) .* b(:,:,:,3);
            G33 = b(:,:,:,3) .* b(:,:,:,3);
            H1 = b(:,:,:,2) .* delta_c;
            H2 = b(:,:,:,1) .* delta_c;
            H3 = b(:,:,:,3) .* delta_c;
            
            % Smooth the elements of the equation system
            G11 = gauss_smoothing(G11, sigma);
            G12 = gauss_smoothing(G12, sigma);
            G13 = gauss_smoothing(G13, sigma);
            G22 = gauss_smoothing(G22, sigma);
            G23 = gauss_smoothing(G23, sigma);
            G33 = gauss_smoothing(G33, sigma);
            H1 = gauss_smoothing(H1, sigma);
            H2 = gauss_smoothing(H2, sigma);
            H3 = gauss_smoothing(H3, sigma);
            
            % Add alpha*identity matrix for stability
            scaleFactor = max([G11(:); G22(:); G33(:)]);
            G11 = G11 + alpha * scaleFactor;
            G22 = G22 + alpha * scaleFactor;
            G33 = G33 + alpha * scaleFactor;
            
            det = G33.*G12.^2 - 2*G12.*G13.*G23 + G22.*G13.^2 + G11.*G23.^2 - G11.*G22.*G33;
            
            % Estimate the displacement field
            update.displacement = cell(3,1);
            update.displacement{1} = removenan((H1.*(G23.^2 - G22.*G33) - ...
                H3.*(G12.*G23 - G13.*G22) - ...
                H2.*(G13.*G23 - G12.*G33))./(det + eps));
            update.displacement{2} = removenan((H2.*(G13.^2 - G11.*G33) - ...
                H3.*(G12.*G13 - G11.*G23) - ...
                H1.*(G13.*G23 - G12.*G33))./(det + eps));
            update.displacement{3} = removenan((H3.*(G12.^2 - G11.*G22) - ...
                H2.*(G12.*G13 - G11.*G23) - ...
                H1.*(G12.*G23 - G13.*G22))./(det + eps));
            
            % Estimate a certainty
            update.certainty = sqrt(b(:,:,:,2).^2 + b(:,:,:,1).^2 + b(:,:,:,3).^2);
    end
end
