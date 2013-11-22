
% The function Gillespie caclulates a stochastic trajectory according to
% the Gillespie algorithm
% Parameters:
%   W:  Transition matrix
%   n0: initial state
% Return values:
%   T:  Vector of time values
%   N:  Vector of states at the times in T

function [T,N] = Gillespie(W, n0)

%Initialize
t = 0;
i = 0;
n = n0;
finished = false;
    while ~finished
    %Determine next reaction to occur
    %First get the non-zero elements (possible destination states)
    %correponding to current state n
    ns = find(W(n,:)>0);
    rates = W(n, ns);
    sumrates = sum(rates);
    if sumrates == 0
        finished = true;
        break;
    end;
    cumrates = cumsum(rates);
    %Select a reaction randomly proportional to its rate
    rn = rand(1);
    rk = find(rn <= cumrates/sumrates,1, 'first'); %find index in rates
    nn = ns(rk);                            %corresponding change in particle number

    n = nn;

    %time interval from exponential distribution
    rt = rand(1);
    dt = -log(rt)/sumrates;

    %Update time and state
    t = t+dt;
    i = i + 1;
    T(i) = t;
    N(i) = n;
    
    %Abort if...
    if i > 10000
        finished = true;
    end;
end;