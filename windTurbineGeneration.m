function [Power_wind]=windTurbineGeneration(windSpeed, windTurbinePara)
v_in = windTurbinePara.cutinSpeed;
v_rate = windTurbinePara.ratedSpeed;
v_out = windTurbinePara.cutoutSpeed;
ratePower = windTurbinePara.ratedPower;

eta = 1;
rho = 1.225;
C_p = 0.4;

Tower_height = 10;
Power_rated_wind = 1;


Wind_tower=windSpeed*(Tower_height/10).^0.142857;
A_swept=2*Power_rated_wind/(eta*rho*C_p*v_rate^3);
Power_wind=zeros(size(windSpeed));

ind=find(Wind_tower<=v_in);
Power_wind(ind)=0;
ind=find(Wind_tower<=v_rate & Wind_tower>v_in);
Power_wind(ind)=1/2*C_p*rho*eta*A_swept*Wind_tower(ind).^3;
ind=find(Wind_tower<v_out & Wind_tower>v_rate);
Power_wind(ind)=Power_rated_wind*eta;
ind=find(Wind_tower>=v_out);
Power_wind(ind)=0;
Power_wind = Power_wind*ratePower;
end