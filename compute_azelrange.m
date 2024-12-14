function [Az,El,Range] = compute_azelrange(userECEF, satECEF)

% userECEF to be 1 by 3
% satECEF to be M by 3
% Y. Wang edited Sep 2024

nsats = size(satECEF,1);

userLLA = ecef2lla(userECEF);
C = Cecef2enu(userLLA(1),userLLA(2));

LOS_ECEF = satECEF - repmat(userECEF,nsats,1);
[Range,LOS_ECEF] = unitvec(LOS_ECEF,2);

LOS_ENU = (C*LOS_ECEF')';

Az = atan2(LOS_ENU(:,1),LOS_ENU(:,2))*180/pi;
El = asind(LOS_ENU(:,3));