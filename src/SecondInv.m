
function  [aII]  =  SecondInv(a)


aII  =  sqrt(1/4.*(a(:,1) - a(:,2)).^2 + a(:,3).^2);

end