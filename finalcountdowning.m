clc; clear all; close all; format long

%definition of steps

%Maximum Stress Failure Theory

    
%% Given values
E_1 =181*10^6;
E_2 = 10.30 *10^6;
V_12 =0.28;
G_12= 7.17;
angle=0:90;
theta= angle*pi/180;
%% positve values of those are tensile negative values are compresssion values.Those values are taken from data sheet adviced on HomeWork sheet.
sigma_1_ult_t = 1500;
sigma_2_ult_t =40;
tao_12_ult_t= 68;

sigma_1_ult_c = -1500;
sigma_2_ult_c =-246;
tao_12_ult_c= -68;

sigma_1_ult_t_2 = 1500/2;
sigma_2_ult_t_2 =40/2;
tao_12_ult_t_2= 68/2;

sigma_1_ult_c_2 = 1500/2;
sigma_2_ult_c_2 =246/2;
tao_12_ult_c_2 = 68/2;

sigma_x_r =4;
sigma_y_r =-3;
tou_xy_r =5;

%%
%matrix definitions

    m_1_t = [sigma_1_ult_t;sigma_2_ult_t;tao_12_ult_t]; %tensile 
    m_1_c = [sigma_1_ult_c;sigma_2_ult_c;tao_12_ult_c]; %comprassion 
    
    m_2 = [sigma_x_r;sigma_y_r;tou_xy_r]; %R value matrix
  

%matrix  m_3 (is our mechanical prob) = T*given values

for i=1:91
    

        % T is transformation matrix
    T = [cosd(i-1)^2, sind(i-1)^2, 2*sind(i-1)*cosd(i-1); sind(i-1)^2 ,cosd(i-1)^2,-2*sind(i-1)*cosd(i-1);-sind(i-1)*cosd(i-1),sind(i-1)*cosd(i-1), cosd(i-1)^2 - sind(i-1)^2] ;
   
         m_3(:,i) = T * m_2; %we made matrix multiplication for finding R values.
         
           
              
       for j =1:3  
          
           if m_3(j,i)<0
            R(j,i) = m_1_c(j)./m_3(j,i);
            
           end     
     
           if m_3(j,i)>0
           
            R(j,i) = m_1_t(j)./m_3(j,i);
           end
       
       end
end



for k = 1:size(R,2)
    
   R_minimums(k) = min(R(:,k));

    
end

    
sigma_x = 4* [R_minimums];
sigma_y = -3* [R_minimums];
tou_xy = 5* [R_minimums];

m_4= [sigma_x;sigma_y;tou_xy];

for k = 1:size(m_4,2)
    
   Stress_max(k) = max(m_4(:,k));
 
end

    
figure
plot(angle,R_minimums,'-g','LineWidth',2);

%%
% In this section we are calculate Tsai-Wu failure theory



%Tsai-Wu Failure Theory

H_1 = 1/sigma_1_ult_t_2 -1/sigma_1_ult_c_2;
H_2= 1/sigma_2_ult_t_2 - 1/sigma_2_ult_c_2;
H_11 = 1/sigma_1_ult_t_2 *1/sigma_1_ult_c_2;
H_22 = 1/sigma_2_ult_t_2 * 1/sigma_2_ult_c_2 ;
H_66 =1/(tao_12_ult_c_2)^2 ;
H_12= -0.5*sqrt(1/(sigma_1_ult_t_2*sigma_1_ult_c_2*sigma_2_ult_t_2*sigma_2_ult_c_2));
 

 % we define syms x because we need to solve this equation in order to find
 % out roots of equations.
 syms x
 
 for f = 1:size(m_3,2)
 sigma_1=x* m_3(1,f); 
 sigma_2 =x* m_3(2,f);
 tao_12 = x*m_3(3,f);
 
 eqn = H_1*sigma_1+H_2*sigma_2+H_11*sigma_1^2+H_22*sigma_2^2+H_66*tao_12^2+2*H_12*sigma_1*sigma_2 == 1;

 S = solve(eqn); % Solver funtion 
 Solution(f) = double(max(S));
 
 end
 
 figure
 plot(angle,Solution,'-r','LineWidth',2)

   