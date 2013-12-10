function Ermito_Splainai
x =[-0.96 -0.3 0.36 0.09 0.09 0.36 -0.3 -0.96];
y =[-0.28 0.65 0.23 -0.11 0.11 -0.23 -0.65 0.28];

t =[0];
for i = 2:length(x)
    t(end+1) = t(end) + ((x(i) - x(i - 1))^2 + (y(i) - y(i - 1))^2)^0.5;
end
t
hold on;
dx = Akima(t, x);
dy = Akima(t, y);
xx = [];
yy = [];
tt = [];
for i = 2:length(t)

    t0 = t(i-1);
    t1 = t(i);
    fx = hermite_spline(t0, t1, x(i-1), x(i), dx(i-1), dx(i));
    fy = hermite_spline(t0, t1, y(i-1), y(i), dy(i-1), dy(i));
    for j = t0:0.01:t1
        xx(end+1) = fx(j);
        yy(end+1) = fy(j);
        tt(end+1) = j;
    end
end

figure(1);
hold on;
grid on;
plot(tt, xx, 'r');
plot(t, x, 'rd');
plot(tt, yy, 'b');
plot(t, y, 'bd');
legend('X(t) grafikas', 'Interpoliavimo taðkai', 'Y(t) grafikas','Interpoliavimo taðkai');
hold off;

figure(2);
hold on;
grid on;
plot(x, y, 'o');
plot(xx, yy, 'b-');
legend('Interpoliavimo taðkai', 'F(x,y)=0 grafikas');

hold off;
end

function f = hermite_spline(x0, x1, y0, y1, dy0, dy1)
    function h=h00(t)
        h = 2 * t^3 - 3 * t^2 + 1;
    end
    function h=h10(t)
        h = t^3 - 2 * t^2 + t;
    end
    function h=h01(t)
        h = -2 * t^3 + 3 * t^2;
    end
    function h=h11(t)
        h = t^3 - t^2;
    end

    %t = (x - x0)/(x1 - x0);
    %f = h00(t) * y0 + h10(t) * (x1 - x0) * dy0 + h01(t) * y1 + h11(t) * (x1 - x0) * dy1
    f = @(x) h00((x - x0)/(x1 - x0)) * y0 + h10((x - x0)/(x1 - x0)) * (x1 - x0) * dy0 + h01((x - x0)/(x1 - x0)) * y1 + h11((x - x0)/(x1 - x0)) * (x1 - x0) * dy1;
end

function DY=Akima(X,Y) % Isvestiniu reiksmiu interpoliavimo taskuose nustatymas pagal skaitinio integravimo formules
    n=length(X);
    f=inline('(2*x-xi-xip1)/((xim1-xi)*(xim1-xip1))*yim1+(2*x-xim1-xip1)/((xi-xim1)*(xi-xip1))*yi+(2*x-xim1-xi)/((xip1-xim1)*(xip1-xi))*yip1');
    for i=1:n
        if i == 1,xim1=X(1);xi=X(2);xip1=X(3); yim1=Y(1);yi=Y(2);yip1=Y(3);DY(i)=f(xim1,xi,xim1,xip1,yi,yim1,yip1);
        elseif i == n, xim1=X(n-2);xi=X(n-1);xip1=X(n); yim1=Y(n-2);yi=Y(n-1);yip1=Y(n); DY(n)=f(xip1,xi,xim1,xip1,yi,yim1,yip1);
        else, xim1=X(i-1);xi=X(i);xip1=X(i+1); yim1=Y(i-1);yi=Y(i);yip1=Y(i+1); DY(i)=f(xi,xi,xim1,xip1,yi,yim1,yip1);
        end
    end
return
end
