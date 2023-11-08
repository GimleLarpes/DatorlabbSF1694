% %-----DEL A-----%
% 
% clear all
% hold off
% 
% inp = input("Skriv in vattentornstryck: ");
% 
% startP = inp;
% 
% p = [startP,1,1,1,0,0];
% q = [0,0,0,0,0,0];
% kIO = [1,2; 2,1; 1,2];
% L = [400, 600, 600, 500, 700, 700];
% k = [0.005, 0.003, 0.003, 0.004, 0.001, 0.001];
% 
% %q(1) = q(2) + q(3)
% %q(2) + q(4) = q(5)
% %q(3) = q(4) + q(6)
% 
% %0 = (startP-p(2)) *k(1)*L(1) - (p(2)-p(4)) *k(2)*L(2) + (p(2)-p(3)) *k(3)*L(3);
% %0 = startP*k(1)*L(1) - p(2)*k(1)*L(1) - p(2)*k(2)*L(2) + p(4)*k(2)*L(2) + p(2)*k(3)*L(3) - p(3)*k(3)*L(3)
% %startP*k(1)*L(1) = p(2)*(k(1)*L(1) + k(2)*L(2) - k(3)*L(3)) - p(4)*k(2)*L(2) + p(3)*k(3)*L(3)
% %-startP = 1/(k(1)*L(1)) * (-p(2)*(k(1)*L(1) + k(2)*L(2) + k(3)*L(3))+p(4)*k(2)*L(2) +p(3)*k(3)*L(3))
% %p(4) = 1/(k(2)*L(2)) * (p(2)*(k(1)*L(1) + k(2)*L(2) + k(3)*L(3)) -p(3)*k(3)*L(3) -startP*k(1)*L(1))
% %Nod2
% 
% %0 = (p(2)-p(3))*k(3)*L(3) -(p(3)-p(4))*k(4)*L(4) -(p(3)*k(6)*L(6))
% %0 = p(2)*k(3)*L(3)-p(3)*(k(3)*L(3) + k(4)*L(4) + k(6)*L(6)) +p(4)*k(4)*L(4)
% %p(2) = 1/(k(3)*L(3)) * (p(3)*(k(3)*L(3) + k(4)*L(4) + k(6)*L(6)) -p(4)*k(4)*L(4))
% %Nod3
% 
% %0 = (p(2)-p(4))*k(2)*L(2) + (p(3)-p(4))*k(4)*L(4) - p(4)*k(5)*L(5)
% %0 = p(2)*k(2)*L(2)-p(4)*(k(2)*L(2)+k(4)*L(4)+k(5)*L(5))+p(3)*k(4)*L(4)
% %p(3) = 1/(k(4)*L(4)) * (p(4)*(k(2)*L(2) + k(4)*L(4) + k(5)*L(5)) - p(2)*k(2)*L(2))
% %Nod4
% 
% 
% b = [startP*k(1)*L(1);0;0];
% 
% A = [k(1)*L(1)+k(2)*L(2)-k(3)*L(3) , k(3)*L(3) , -k(2)*L(2);
%      k(3)*L(3) , -(k(3)*L(3)+k(4)*L(4)+k(6)*L(6)) , k(4)*L(4);
%      k(2)*L(2) , k(4)*L(4) , -(k(2)*L(2)+k(4)*L(4)+k(5)*L(5))];
% 
% p(2:4) = (A\b)';
% 
% 
% T = table([1,2,3,4,5,6]',p','VariableNames',["Nod Nr:", "Tryck: "]);
% disp(T)
% 
% plot(p)
% hold on
% xlabel("Nod Nr")
% ylabel("Tryck")
% 
% disp(["Genomsnittstryck: ",mean(p(2:4))])

%-----DEL B----%

clear all
hold off
tiledlayout(2,1)

inp = input("\nVilken fil vill du ladda in? ","s");
data = load(inp);
A = data.A;
X = data.X;
Y = data.Y;

k = 0.01;
l = 200;

[L, U, P] = lu(A);

startP = [];

cont = true;

while cont == true
    disp(" ")
    disp("Skriv tryck för resorvar...")
    for i = [1:8]
        startP(i) = input(append("Nr. ", string(i), ": "));
    end
    
    V = [-k*l.*startP]';
    
    b = [V;zeros(6400-8,1)];
    
    y = L\(P*b);
    res = U\y;
    %res = A\b;
     
    disp(" ")
 
    T = table(min(res),max(res),mean(res),'VariableNames',["Min", "Max", "Mean"]);
    disp(T)
    
    inp = input("För nod-tryck graf, skriv 1\nFör 3d-graf, skriv 2\nFör bägge, skriv 3\n");
    
    if inp == 1
        plot(res)
        hold on
        xlabel("Nod Nr")
        ylabel("Tryck")
    elseif inp == 2
        coords = reshape(res,size(X));
        surf(X,Y,coords);
    else
        nexttile()
    
        plot(res)
        hold on
        xlabel("Nod Nr")
        ylabel("Tryck")
    
        nexttile()
        
        coords = reshape(res,size(X));
        surf(X,Y,coords);
    end

    inp = input("\nOm du vill fortsätta, skriv 1: ","s");
    if inp == "1"
        cont = true;
    else
        cont = false;
    end
end