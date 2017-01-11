foreach n (crab_JetHT*V6_v1)
#crab report $n
cp $n/results/processedLumis.json ./$n.json
end

