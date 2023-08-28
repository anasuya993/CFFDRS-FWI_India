%% Code to compute Fire Weather Index (FWI) for India
% Anasuya Barik, IIT Delhi
% contact the author at anasuya993@gmail.com

% Ths code computes the CFFDRS-FWI for 1 year at a daily temporal
% resolution.  
% This index involves three moisture codes, the Fine Fuel Moisture Code (FFMC), the Duff Moisture Code (DMC), and the Drought Code (DC). 
% These codes quantify moisture content at depths of 0-1 inch, 2-4 inches, and 4-8 inches, respectively, considering both rainfall and drying phases. 
% The algorithm then proceeds to derive two intermediary indices, the Initial Spread Index (ISI) and the Buildup Index (BUI). 
% The ISI integrates wind speed and the top layer moisture from FFMC to determine the fire's rate of spread. Meanwhile, BUI, computed from DMC and DC, reflects available fuel for combustion. 
% Ultimately, the FWI is obtained from the combination of ISI and BUI.  


% The inputs to the CFFDRS FWI are 2m temperature (K or Degree Celcius), relative
% humidity (%), past 24hr accumulated precipitation (mm), and wind speed (m/s) recorded at 12noon of
% the local time. In the absence of such sub-daily data, daily gridded
% data of mean/maximum temperature, minimum relative humidity, mean wind speed, and daily
% accumulated total rainfall can also be used.

% Few IMPORTANT notes for using the code
% ix is the no of grids in the x direction 
% iy is the no of grids in the y direction 
% The input variables should ideally be in the dimension ix x iy x 365
% (considering a 365-day calendar)for gridded FWI computation. For point
% locations the input files must be 1D array of size 1 x 365.

% NOTE 
% The code is primarily developed for India. However, it is suitable for usage in any tropical region of interest. 
% This code has been developed using CFFDRS-FWI equations as described in Wagner, 1987. Refer to README.doc for a detailed description of the adjustments done.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Read latitude longitude information
% For point location read the respective location values
lon=ncread('filename.nc','longitude'); % 1D matrix
lat=ncread('filename.nc','latitude');   % 1D matrix

%% Input data read 
temp=ncread('filename.nc','temp_variable');
temp_K=temp_111-273.15;

rh_input=ncread('filename.nc','rh_variable');

prec_input=ncread('filename.nc','pr_variable');

ws_input=ncread('filename.nc','uv_magnitude_variable'); % computed from the zonal and meridional components

%% FFMC, DMC, DC initialization (Refer README doc for spin up details)

% initialization for spin up
ffmc_yda(1:ix,1:iy) = 85;       % standard initialization values for day 1 
dmc_yda(1:ix,1:iy)=6;
dc_yda(1:ix,1:iy)=150;

% load previous year ffmc, dmc, dc     (for warm start of the model)             
ffmc_yda=ffmc_prev_year(:,:,365);     %december 31st value of previous year
dmc_yda=dmc_prev_year(:,:,365);
dc_yda=dc_prev_year(:,:,365);

% Initialize current run year indices                               %Change
ffmc(1:ix,1:iy,1:365)=NaN;
dmc(1:ix,1:iy,1:365)=NaN;
dc(1:ix,1:iy,1:365)=NaN;
fwi(1:ix,1:iy,1:365)=NaN;

%% FFMC, DMC, DC, FWI computation
for ii=1:365
   % Provide meteorological data for current run year
   temp= temp_K(:,:,ii);                                     %Change
   rh= rh_input(:,:,ii);
   rh(rh>=100)=100;
   ws=ws_input(:,:,ii);         % Bias Corrected with GSOD
   ws(ws<=0)=0;
   prec=prec_input(:,:,ii);
  
wmo=147.2.*(101-ffmc_yda)./(59.5+ffmc_yda);   % eq 2b
WMO=wmo;                %accounting wmo for later use
%rain reduction to allow for loss in canopy head
ra1=prec;               %prec in mm
ra1(ra1<=0.5)=NaN;      
ra1=ra1-0.5;            %rainfall upto 0.5mm are to be ignored (eq 14 in wagner et al, 1987)
ra11=ra1;               %accounting ra1 for later use
ra2=prec;               %avoiding negative values in rainfall loss calculation
ra2(ra2>0.5)=NaN;
ra1(isnan(ra1))=ra2(isnan(ra1));
ra=ra1;                 %rainfall loss corrected precipitation

%masking values in wmo with ra1 and ra2
wmo(isnan(ra1))= NaN;   
wmo1=wmo;
wmo=WMO;
wmo(isnan(ra2))= NaN;
wmo2=wmo;
%effect of rainfall
%when initial moisture content is more than 150
wmo11=wmo1;
wmo11(wmo11<=150)=NaN;
wmo11=wmo11 + 0.0015.*(wmo11 - 150).*(wmo11 - 150).* sqrt(ra11) + 42.5.*ra11.*exp(-100./(251 - wmo11)).*(1 - exp(-6.93 ./ ra11)); %eq 12 and 13
%when initial moisture content is less than 150
wmo12=wmo1;
wmo12(wmo11>150)=NaN;
wmo12=wmo12 + 42.5.*ra11.*exp(-100./(251 - wmo12)).*(1 - exp(-6.93 ./ ra11)); %eq 12 and 13 in wagner et al, 1987
wmo_1=wmo11;
wmo_1(isnan(wmo_1))=wmo12(isnan(wmo_1));    %rainfall corrected increase in moisture content
wmo1(isnan(wmo1))=wmo2(isnan(wmo1));        %actual moisture content
%capping the actual moisture content to 250 percent
wmo1(wmo1>250)=250;
% Equilibrium moisture content from drying
ed=0.942.*(rh.^0.679)+(11.*exp((rh - 100)./10))+0.18.*(21.1 - temp).*(1 - 1 ./ exp(rh .* 0.115));   %eq 8a in wagner et al, 1987
% Equilibrium moisture content from drying
ew=0.618.*(rh.^0.753)+(10.*exp((rh - 100)./10))+0.18.*(21.1 - temp).*(1 - 1 ./ exp(rh .* 0.115));   %eq 8b in wagner et al, 1987
% Log drying rate at the normal termperature of 21.1 C
for i=1:ix
    for j=1:iy
        if wmo1(i,j)<ed(i,j) && wmo1(i,j)<ew(i,j)
        z(i,j)=0.424*(1-(((100 - rh(i,j))/100)^1.7))+0.0694*(ws(i,j)^0.5)*(1-((100 - rh(i,j))/100)^8);     %eq 4 in wagner et al, 1987
        x(i,j)=z(i,j)*0.581*exp(0.0365*temp(i,j));              %effect of temperature on  drying rate
        wm(i,j)=ew(i,j)-(ew(i,j)-wmo1(i,j))/(10^x(i,j));
        else z(i,j)=0; x(i,j)=0; wm(i,j)=wmo1(i,j);
        end
    end
end
% Log wetting rate at the normal termperature of 21.1 C
for i=1:ix
    for j=1:iy
        if wmo1(i,j)>ed(i,j) 
        z(i,j)=0.424*(1-(rh(i,j)/100)^1.7)+0.0694*sqrt(ws(i,j))*(1-(rh(i,j)/100)^8);     %eq 5 in wagner et al, 1987
        wm(i,j)=ed(i,j)+(wmo1(i,j)-ed(i,j))/(10^x(i,j));                                    
        else z(i,j)=z(i,j);  wm(i,j)=wm(i,j);
        end
        x(i,j)=z(i,j)*0.581*exp(0.0365*temp(i,j));              %Effect of temperature on  drying rate
    end
end 
ffmc_cal=(59.5 .* (250 - wm))./(147.2 + wm);                                               %eq 2a
ffmc_cal(ffmc_cal<0)=0;
ffmc_yda=ffmc_cal;              % Saving ffmc_yda for next day initialization
ffmc(:,:,ii)=ffmc_cal;

%*************************************************************************************************************
% DMC calculation
%Rainfall phase of DMC
ra=prec;                                    %observed rainfall in the open
ra(ra<1.5)=0;                               %rainfall upto 1.5mm are ignored for DMC calc
rw=0.92.*ra-1.27;                           %effective rain
wmi= 20 + 280./exp(0.023.*dmc_yda);
if dmc_yda <= 33
    b=100./(0.5+0.3.*dmc_yda);
elseif dmc_yda <=65 & dmc_yda>33
    b=14-1.3.*log(dmc_yda);
else
    b=6.2.*log(dmc_yda)-17.2;
end
wmr=wmi+1000.*rw./(48.77+b.*rw); %moisture content after rain
wmr(wmr<21)=21;
pr=43.43.*(5.6348-log(wmr-20));
pr(pr<0)=0;         %constrain P
pr(isinf(pr)|isnan(pr)) = 0;
%Drying phase of DMC
%log drying rate (k)
if ii<=31
    index=1;
elseif ii>31 && ii<=59
    index =2
elseif ii>59 && ii<=90
    index =3
elseif ii>90 && ii<=120
    index =4
elseif ii>120 && ii<=151
    index =5
elseif ii>151 && ii<=181
    index =6
elseif ii>181 && ii<=212
    index =7
elseif ii>212 && ii<=243
    index =8
elseif ii>243 && ii<=273
    index =9
elseif ii>273 && ii<=304
    index =10
elseif ii>304 && ii<=335
    index =11
elseif ii>335 && ii<=365
    index =12
end
for jj=1:ix
    if lat(jj)<30  
        le=[7.9;8.4;8.9;9.5;9.9;10.2;10.1;9.7;9.1;8.6;8.1;7.8]; %DMC day length adjustment
        k=1.894.*(temp+1.1).*(100-rh).*le(index).*10^-6;        %Eq. 20
    elseif lat(jj)>=10
        le=[7.9;8.4;8.9;9.5;9.9;10.2;10.1;9.7;9.1;8.6;8.1;7.8]; %DMC day length adjustment
        k=1.894.*(temp+1.1).*(100-rh).*le(index).*10^-6;        %Eq. 20
    elseif lat(jj)>=30 
        le=[6.5;7.5;9;12.8;13.9;13.9;12.4;10.9;9.4;8;7;6];      %DMC day length adjustment
        k=1.894.*(temp+1.1).*(100-rh).*le(index).*10^-6;        %Eq. 20
    else
        k=1.894.*(temp+1.1).*(100-rh).*9.*10^-6; % for latitudes close to equator
    end
end
dmc_today=pr+100.*k;
dmc_yda=dmc_today;
dmc(:,:,ii)=real(dmc_today);

%*************************************************************************************************************
% DC calculation
t0=temp;
t0(t0<-2.8)=-2.8;                           %constrain temperature
%DC day length adjustment
for jj=1:ix
if lat(jj)>20 
    f1=[-1.6,-1.6,-1.6,0.9,3.8,5.8,6.4,5,2.4,0.4,-1.6,-1.6]; 
    pe=(0.36.*(temp+2.8)+f1(index))/2;      %potential evapotranspiration
else
    pe=(0.36.*(temp+2.8)+1.4)./2;           %for latitudes close to equator we consider constant index
end
end
smi=800.*exp(-1.*dc_yda./400);
pe(pe<0)=0;                                 %constrain pe
ra=prec;
rw=0.83.*ra-1.27;                           %effective rainfall
for i=1:ix
    for j=1:iy
        if ra(i,j)<=2.8
        dr(i,j)=dc_yda(i,j);
        else
        dr(i,j)=real(dc_yda(i,j)-400.*log(1+3.937.*rw(i,j)./smi(i,j)));
        end
    end
end
dc1=dr+pe;
dc1(dc1<0)=0;
dc_yda=dc1;
dc(:,:,ii)=dc1;

%**********************************************************************************************************
% Initial Spread Index (ISI)
fw=exp(0.05039.*ws);                  %wind effect
fm=147.2.*(101-ffmc)./(59.5+ffmc);    %moisture effect
ff=91.9.*exp(-0.1386.*fm).*(1+(fm.^5.31)./49300000); %Fine Fuel Moisture
isi=0.208.*fw.*ff;                    %spread index equation

%**********************************************************************************************************
% Buildup Index (BUI)
if any(dmc_today<=0.4*dc1)
    bui=0.8.*dc1.*dmc_today./(dmc_today+0.4.*dc1);
else
    bui=dmc_today-(1-0.8.*dc1./(dmc_today+0.4.*dc1)).*(0.92+(0.0114.*dmc_today).^1.7);
end
bui(dmc_today==0)=0;

%**********************************************************************************************************
% Fire Weather Index (FWI)
for i=1:369
    for j=1:369
if bui(i,j)<=80
    fd1(i,j)=real(0.626.*(bui(i,j)).^(0.809) + 2);
else
    fd1(i,j)=real(1000./(25 + 108.64.*exp(-0.023.*bui(i,j))));
end
    end
end
bb=0.1.*isi.*fd1;
if any(bb>1)
    fwi_calc=real(exp(2.72.*(0.434.*log(bb)).^0.647));
else
    fwi_calc=bb;
end
fwi(:,:,ii)=fwi;

disp(ii)
end

% SAVE THE VARIABLES ffmc, dmc, dc, fwi from the workspace into .mat format
% or find the code to create nc files in https://github.com/anasuya993/postprocessing_DSCESM/blob/main/DSCESM_postprocessing_ncwrite.m


%% Plot and Check the output
pcolor(lon,lat,(squeeze(nanmean(fwi,3)))'); shading interp; 
hold on;
% plot([ss.X],[ss.Y],'-k') % Use shapefile if required
colorbar;
title('FWI gridded - Annual mean ')
set(gca,'FontSize',16, 'FontName','Times New Roman','FontWeight','Bold','Layer','top')


