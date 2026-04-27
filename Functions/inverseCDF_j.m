function [AgeTransit] = inverseCDF_j(Cell,Probability)
CellRand = rand(1,Cell);
Age = 0:0.1:(10*Probability.Mean);
Prob =Probability.Scaling.*normpdf(Age,Probability.Mean,Probability.STD);
Cumulative = 1-cumprod(1-Prob);
Result = ~bsxfun(@lt,CellRand,Cumulative');
idx = (sum(Result,1));
idx(idx==0) = 1;
AgeTransit = Age(idx) ;
end