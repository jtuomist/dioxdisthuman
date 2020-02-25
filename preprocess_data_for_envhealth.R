library(thlVerse)
library(thlConnect)

dbi <- thlDbInfo("pubhealth", dbengine="postgresql")
tst <- thlDbQuery(dbi,"SELECT * FROM envhealth.valvontayksikko_kunta")

thlDbConnect("pubhealth",dbengine = "postgresql")
