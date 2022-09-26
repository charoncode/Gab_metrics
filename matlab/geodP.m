function geo = geodP(q1d,q2d,t,a,dime)

if dime>2
%     % d'abord on regarde le plan commun à q1d et q2d
    [r1, theta1, phi1] = cart3sph(q1d);
    [r2, theta2, phi2] = cart3sph(q2d);

    q12 = [r2*cos(theta2-theta1)*cos(phi2), r2*sin(theta2-theta1)*cos(phi2), r2*sin(phi2)];
    Ry(1, :) = [cos(phi1) 0 sin(phi1)];
    Ry(2, :) = [0 1 0];
    Ry(3, :) = [-sin(phi1) 0 cos(phi1)];
    q12 = (Ry*q12.').';
    
    [theta12, r12] = cart2pol(q12(2),q12(3));
    Rx(1, :) = [1 0 0];
    Rx(2, :) = [0 cos(-theta12) -sin(-theta12)];
    Rx(3, :) = [0 sin(-theta12) cos(-theta12)];
    q12 = (Rx*q12.').';
    

    q1=[r1,0];
    q2=[q12(1),q12(2)];
else
    q1=q1d;
    q2=q2d;
end


    [theta1, r1] = cart2pol(q1(1),q1(2));
    [theta2, r2] = cart2pol(q2(1),q2(2));


    if theta2-theta1 > pi
        theta2 = theta2-2*pi;
    end
    if theta2-theta1 < -pi
        theta2 = theta2+2*pi;
    end

    if theta2-theta1>0
        if a*(theta2-theta1)<pi
            q1temp = [r1,0];
            q2temp = [r2*cos(a*(theta2-theta1)), r2*sin(a*(theta2-theta1))];
            x = (1-t)*q1temp + t*q2temp;
            [thetat, rt] = cart2pol(x(1),x(2));
            geo = [rt*cos(1/a*thetat+theta1),rt*sin(1/a*thetat+theta1)];
        else
            if t<=r1/(r1+r2)
                geo = (1-(r1+r2)/r1*t)*q1;
            else
                geo = (t-r1/(r1+r2))*(r1+r2)/r2*q2;
            end
        end
    else
        if a*(theta2-theta1)>-pi
            q2temp = [r2,0];
            q1temp = [r1*cos(a*(theta1-theta2)), r1*sin(a*(theta1-theta2))];
            x = (1-t)*q1temp + t*q2temp;
            [thetat, rt] = cart2pol(x(1),x(2));
            geo = [rt*cos(1/a*thetat+theta2),rt*sin(1/a*thetat+theta2)];
        else
            if t<=r1/(r1+r2)
                geo = (1-(r1+r2)/r1*t)*q1;
            else
                geo = (t-r1/(r1+r2))*(r1+r2)/r2*q2;
            end
        end 
    end

if dime>2
    [r1, theta1, phi1] = cart3sph(q1d);
    [r2, theta2, phi2] = cart3sph(q2d);

    Rx(1, :) = [1 0 0];
    Rx(2, :) = [0 cos(theta12) -sin(theta12)];
    Rx(3, :) = [0 sin(theta12) cos(theta12)];
    geo = ([geo(1) geo(2) 0]).';
    geo = Rx*geo;

    Ry(1, :) = [cos(-phi1) 0 sin(-phi1)];
    Ry(2, :) = [0 1 0];
    Ry(3, :) = [-sin(-phi1) 0 cos(-phi1)];
    geo = Ry*geo;
    geo = geo.';
    [rg, thetag, phig] = cart3sph(geo);
    geo = [rg*cos(thetag+theta1)*cos(phig), rg*sin(thetag+theta1)*cos(phig), rg*sin(phig)];
end

end