function [NCelltoadd] = CellIntroduction (introductiontype,IntroductionRate,t,steptime,N0,TotalCelladded)

TotalCellThisTimestep = round((t + steptime)^introductiontype * N0 / (IntroductionRate^introductiontype));
NCelltoadd = max(TotalCellThisTimestep - TotalCelladded,0);
if (t==0 && NCelltoadd == 0)
    NCelltoadd = 1;
end
if (NCelltoadd + TotalCelladded > N0)
    NCelltoadd = N0 - TotalCelladded;
end

end