function E=calc_TI(E1,E2)
E1_norm=norm(E1);E2_norm=norm(E2);
cos_phi=dot(E1,E2)/(E1_norm*E2_norm);
if E1_norm<E2_norm
    temp=E1_norm;
    E1_norm=E2_norm;
    E2_norm=temp;
end

if cos_phi+1<0.00001
    E=0;
elseif E1_norm*cos_phi<E2_norm 
    sin_phi=sqrt(1-cos_phi*cos_phi);
    tan_theta=(E1_norm-E2_norm*cos_phi)/(E2_norm*sin_phi);
    cos_theta=sqrt(1/(1+tan_theta^2));
    sin_theta=sqrt(1-cos_theta^2);
    E=min(E1_norm*cos_theta,E2_norm*(cos_theta*cos_phi+sin_theta*sin_phi));
else
    E=E2_norm;
end
end