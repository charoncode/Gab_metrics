function geo_deriv = geodP_deriv(q1,q2,t,a)
    [theta1, r1] = cart2pol(q1(1),q1(2));
    [theta2, r2] = cart2pol(q2(1),q2(2));
    
    if t==0 | t==1
        geo_deriv = [0,0];
    else
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
            rtderiv = (-2*t*(1-t)*r1*r2*(theta2-theta1)*sin(a*(theta2-theta1)))/(2*rt);
            xtemp= ((1-t)*r1+t*r2*cos(a*(theta2-theta1)))/rt;
            thetatderiv = t*(r2*rt*(theta2-theta1)*sin(a*(theta2-theta1))+r2*cos(a*(theta2-theta1))*rtderiv)/(rt^2)*1/sqrt(1-xtemp^2);
            geo_deriv = [rtderiv*cos(1/a*thetat+theta1)-rt*((thetatderiv*a-thetat)/a^2)*sin(1/a*thetat+theta1),rtderiv*sin(1/a*thetat+theta1)+rt*((thetatderiv*a-thetat)/a^2)*cos(1/a*thetat+theta1)];
        else
            geo_deriv = [0,0];
        end
    else
        if a*(theta2-theta1)>-pi
            q2temp = [r2,0];
            q1temp = [r1*cos(a*(theta1-theta2)), r1*sin(a*(theta1-theta2))];
            x = (1-t)*q1temp + t*q2temp;
            [thetat, rt] = cart2pol(x(1),x(2));
            rtderiv = (-2*t*(1-t)*r1*r2*(theta1-theta2)*sin(a*(theta1-theta2)))/(2*rt);
            xtemp= (t*r2+(1-t)*r1*cos(a*(theta1-theta2)))/rt;
            thetatderiv = (1-t)*(r1*rt*(theta1-theta2)*sin(a*(theta1-theta1))+r1*cos(a*(theta1-theta2))*rtderiv)/(rt^2)*1/sqrt(1-xtemp^2);
            geo_deriv = [rtderiv*cos(1/a*thetat+theta2)-rt*((thetatderiv*a-thetat)/a^2)*sin(1/a*thetat+theta2),rtderiv*sin(1/a*thetat+theta2)+rt*((thetatderiv*a-thetat)/a^2)*cos(1/a*thetat+theta2)];
        else
            geo_deriv = [0,0];
        end 
    end
    end
geo_deriv = geo_deriv.';
end