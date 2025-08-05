% This function takes an input of a rotation matrix, R (3x3), and outputs
% the equivalent axis of rotation, u (3x1), and angle of rotation, theta.
% Date: 2/9/2025
function [u,theta] = EquivalentAxis3(R)
            zeroThreshold = 1e-6;
            if abs(trace(R)-3)<zeroThreshold % R=I, theta=0, u = NaN
                theta = 0;
                u = nan;
            elseif abs(trace(R)+1)<zeroThreshold %theta=pi
                theta = pi;
                % [u]^2 = u*u' - I
                % diagonals of ([u]^2+I) are ui^2
                uuTranspose = (R+eye(3))/2;
                % find the maximum value of ui^2
                [mm,ii] = max(diag(uuTranspose)); % maximum ui^2
                % choose the column with maximum ui and get u
                u = uuTranspose(:,ii)/sqrt(mm);  
            else %theta is neither pi or 0
                % 
                theta = acos((trace(R)-1)/2);
                u = so32vec(1/(2*sin(theta))*(R-R'));
            end
        end