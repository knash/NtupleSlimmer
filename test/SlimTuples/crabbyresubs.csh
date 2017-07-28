foreach n (crab_*V12)
crab status $n
crab resubmit $n --siteblacklist=T1_RU_JINR,T2_RU_JINR,T3_US_UCR,T2_US_Florida
end

