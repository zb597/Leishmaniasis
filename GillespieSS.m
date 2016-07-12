function [new_state,Use]=GillespieSS(old, Parameters)
%Uses the tau leaping gillespie algorithm to calculate populations of S,I,R at
%t+timestep given populations S,I,R at time t.

%initialise values

AA=Parameters(1);             %total infection rate
BB=Parameters(2);             %recovery rate
CC=Parameters(3);             %birth rate & death rate
N=Parameters(4);              %total population
timestep=Parameters(5);       %timestep
super=Parameters(6);          %proportion of superspreaders
r2=Parameters(7);             %avg number of secondardy infections due to superspreaders
X=old(1); Y=old(2); Z=old(3); %inital sus, infec, recov

%Define system of ODE's for populations rate of change
r1=AA-super*r2;               %calculate 'normal' infection rate s.t AA=r1+super*r2

PopRate(1) = CC*N;            %births
RateOfChange(1,:)=[+1 0 0];   %births applied to S class
PopRate(2) = r1*X*Y/N;        %normal infections 
RateOfChange(2,:)=[-1 +1 0];  %infections removed from S and added to I
PopRate(3) = BB*Y;            %recovery
RateOfChange(3,:)=[0 -1 +1];  %removed from I added to R
PopRate(4) = CC*X;            %death of S
RateOfChange(4,:)=[-1 0 0];   %removed from S
PopRate(5) = CC*Y;            %death of I
RateOfChange(5,:)=[0 -1 0];   %removed from I
PopRate(6) = CC*Z;            %death of R
RateOfChange(6,:)=[0 0 -1]; %removed from R
PopRate(7) = super*X*Y/N;  %infections from SSE       
RateOfChange(7,:)=[-1 +1 0];  %infections removed from S and added to I

new_state=old; %store current population

%calculate new population using Gillespie, number of objects added/removed 
%     from a given class is a poisson r.v with mean = rate of change of pop.*timestep 

for i=1:7
    Num=poissrnd(PopRate(i)*timestep);                    %pick number of objects to change state from poission distribution with rate=pop rate*dt  
    
    if(i==7 && Num>0) %SSE case
        inf=0; %reset number of SS infections
        for j=1:Num %for each SSE
            Num2=poissrnd(r2); %number of infections caused by jth SS
            inf=inf+Num2; %keep track of total num of new infections 
        end
        Num=inf;
    end
    
    Use(i)=min([Num new_state(find(RateOfChange(i,:)<0))]);  %use = number to add or remove from class, 
                                                          %if use>number objects in class when removing,  remove remaining number of objects
                                                          
                                                       
    new_state=new_state+RateOfChange(i,:)*Use(i);            %add or remove 'use' number of objects from class S, I or R    
end


