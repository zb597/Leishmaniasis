%function [] = EpidemicSS(AA,BB,CC,X0,Y0,N0,r2,super ,timestep, Tmax)

%Function epidemicSS is a stochastic model of a disease outbreak with three classes Suseptibles(S),
%Infected(I) and recovered(R). Objects move between classes through normal infection, superspreading(SSE), recovery, 
%birth & death. Note that birth & death rate are equal for all classes, births are added to the S class, 
%only S objects can become infected(I). 
%
%Input values: Total Infection rate AA(/day),
%              Recovery rate BB(/day) 
%              Death rate CC(/day)
%              Initial S populations X0 
%              Initial I population Y0
%              Initial Total population N0
%              Intensity of SSE's
%              mean number of infections per SSE
%              Integration timestep(days)
%              Total run time(days)
%
%  If user does not define initial values, function will start from pre-defined values
clear all;

store=0;
n=5;
for i=1:n
hold on;
%if nargin == 0 %if no input arguments given, initialise rate parameters and inital S, I and R populations
   AA=2.5;       %avg infection rate(/day)
   BB=1/10.0;  %recovery rate
   CC=5e-4;    %birth/death rate
   N0=5000;    %total pop
   Y0=ceil(CC*N0/BB);  %initial infectious
   X0=floor(BB*N0/AA); %initial suseptiable
   timestep=1;
   Tmax=2*365; %run for 2 years
   super=0.2; %proportion of infected which are superspreaders
   r2=8;   %mean number of infections caused by each superspreader
%end

Z0=N0-X0-Y0;  %recovered = total - suseptiable - infectious
[t, pop, infect, spread, hittime]=SIRmodelSS([0 Tmax],[X0 Y0 Z0],[AA BB CC N0 timestep, super, r2]);
T=t/365; %convert time to years
XX=pop(:,1);  %susceptible
YY=pop(:,2);  %infectious
ZZ=pop(:,3);  %recovered
hold on;
subplot(3,1,1);
h=plot(T,XX,'-g');
xlabel 'Time in years'; %plot susceptible population XX
ylabel 'Susceptible';
hold on;
subplot(3,1,2);
h=plot(T,YY,'-r');
xlabel 'Time in years'; %plot infectious population YY
ylabel 'Infectious';
hold on;
subplot(3,1,3);
h=plot(T,ZZ,'-k');    %plot recovered population ZZ
xlabel 'Time in years';
ylabel 'Recovered';



(spread+infect)/hittime;

end

