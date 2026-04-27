function [AgeTransit] = inverseCDF(Cell,Probability)
CellRand = rand(1,Cell);
Age = 0:0.1:(10*Probability.Mean);
Prob =Probability.Scaling.*normpdf((Age),Probability.Mean,Probability.STD);
Cumulative = 1-cumprod(1-Prob);
Result = ~bsxfun(@lt,CellRand,Cumulative');
idx = min((sum(Result,1)+1),size(Age,2));
AgeTransit = Age(idx);
end