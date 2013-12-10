function Aproksimavimas_Vienanariu_Bazeje
format long;
gx = @(x) exp(-x) * cos(x)/(x - 6);
px = []
for i = -5:0.1:5
    px(end + 1) = gx(i);
end
plot(-5:0.1:5, px);
hold on;

x = -5:0.5:5;
y = [];
for i = x
    y(end + 1) = gx(i);
end
t =[0];
for i = 2:length(x)
    t(end+1) = t(end) + ((x(i) - x(i - 1))^2 + (y(i) - y(i - 1))^2)^0.5;
end
Gt = base(6, t);
coeff_x = ((Gt'*Gt)\(Gt'*x'))';
coeff_y = ((Gt'*Gt)\(Gt'*y'))';
coeff_x = fliplr(coeff_x);
coeff_y = fliplr(coeff_y);

fx = matlabFunction(poly2sym(coeff_x))
fy = matlabFunction(poly2sym(coeff_y));
xx = [];
yy = [];
for i = (t(1)-0.5):0.01:(t(end)+0.5)
    xx(end+1) = fx(i);
    yy(end+1) = fy(i);
end
 paklaidos_x = [];
 for i = 2:10
     paklaidos_x(end+1) = tikrinti(gx, x, y, i);
 end 
plot(x, y, 'o');
plot(xx, yy, 'r-');
legend('aproksimuojama funkcija', 'aproksimavimo taskai', 'aproksimuota funkcija');
hold off;
figure(2)
title('Paklaidø priklausomybës nuo aproksimavimo bazës grafikas');
plot(2:10, paklaidos_x);
legend('paiklaidos dydis kai aproksimuojamø taðkø skaièius yra 21');
xlabel('Baziniø funkcijø skaièius');
ylabel('paklaida');
    
end

function G=base(m, x)
    G = [];
    for i = 1:m, G(:, i) = x.^(i-1); end
    return
end

function c=paklaida(f, x, y)
    c = 0;
    for i = 1:length(x)
        c = c + (f(x(i)) - y(i))^2;
    end
    c = c/2
end

function p=tikrinti(funkcija, x, y, laipsnis)
    Gx = base(laipsnis, x);
    coeff = ((Gx'*Gx)\(Gx'*y'))';
    coeff = fliplr(coeff);    
    fx = matlabFunction(poly2sym(coeff));
    p = paklaida(fx, x, y);
end