function [T,P, infect, spread, hittime]=SIRmodelSS(Time,Initial,Parameters)

%Inputs values of population size and parameter values from Epidemic and performs
%stochasic integration using the Gillespie algorithm for t=0 to t=Tmax with
%timestep given in Epidemic. 
%Returns a matrix of each population size at each timestep

Xstart=Initial(1);      %initial S
Ystart=Initial(2);      %initial I
Zstart=Initial(3);      %initial R
timestep=Parameters(5); %timestep
infect=0; %number of normal infections
spread=0; %number of SS infections
hittime=Time(2);

T=[0:timestep:Time(2)];        %time vector from 0 to Total time in timestep
P(1,:)=[Xstart Ystart Zstart]; %inital population values
old=[Xstart Ystart Zstart];    %inital population values

loopcount=1;

while (T(loopcount)<Time(2))         %while time<Total time
    [new,obj]=GillespieSS(old,Parameters); %calculate new pop values using Gillespie algorithm
    loopcount=loopcount+1;           %next time value
    P(loopcount,:)=new;              %assign new pop values to vector P
                           
    infect=infect+obj(2);       %increase number of normal infections
    spread=spread+obj(7);       %increase number of SS infections

    if(new(2)==0)hittime=T(loopcount); end %calculate time until infection dies out (if it doesn't then hitime=Tmax)
    
    old=new;  %re-initialise starting value with new pop values
end

