foreach n (crab_JetHT*V9)
crab status $n
crab resubmit $n --siteblacklist=T1_RU_JINR,T2_RU_JINR,T3_US_UCR,T2_US_Florida
end

