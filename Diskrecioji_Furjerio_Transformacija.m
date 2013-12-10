function Diskrecioji_Furjerio_Transformacija
T = 4;
F =@(t) sign(sin(2 * pi  * t/T)) * cos(2 * pi * 2 * t/T);
R =@(t) F(t) + 0.15 * cos(2*pi * 110 * t/T) + 0.13 * cos(2 * pi * 100 * t/T)
%ezplot(R);
t = 0:4/1000:4; 
y = [];
for i = t
    y(end+1) = R(i);
end

plot(t,y) 
grid on 
title('Duota funkcija su triukðmu') 

Y= slow_dft(y); 
spektras = abs(Y);
hold off;
figure(6);
title('Furje aproksimacijos harmonikø amplitudës su paþymëtomis slenksèiø reikðmëmis');
grid on;
hold on;
bar([0:length(spektras)-1], spektras, 0.01);

Pyy = sqrt(Y.* conj(Y)); 

filtras = 0.285;
nuo = round(filtras * 250);
iki = round((4 - filtras)* 250);

ind = nuo:iki;
%ind = Pyy<73.5;% tinkamas slenksèio dydis 
line([nuo, nuo], [0, 73.5]);
line([-10, 1010], [73.5, 73.5],'Color', 'r');
line([iki, iki], [0, 73.5]);
legend('Furje amplitudþiø spektras pagal modulá', 'filtravimo pagal daþná riba', 'filtravimo pagal amplitudæ riba');
hold off;
Y(ind)=0; 
yt=slow_idft(Y); 

figure(4) 

title('funkcija be triukðmo') 
plot(t,yt) 
grid on;
end

function y=slow_dft(x)
    n = length(x);
    root = exp(-2 * pi * 1i / n);
    y = zeros(n, 1);
    for p = 1:n
        y(p) = 0;
        for j = 1:n
            y(p) = y(p) + x(j) * root ^ (j * p);
        end        
    end
    
end

function y=slow_idft(x)
    n = length(x);
    y = zeros(n, 1);    
    root = exp(2 * pi * 1i / n);
    
    for p = 1:n
        y(p) = 0;
        for k = 1:n
            y(p) = y(p) + x(k) * root ^ (k * p); 
        end
        y(p) = y(p)/n;
    end
end
