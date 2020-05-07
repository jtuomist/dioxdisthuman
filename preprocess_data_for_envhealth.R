# This code is github.com/jtuomist/dioxdisthuman.git/preprocess_data_for_envhealth.R
# It takes data from three sources and creates a coherent data table for further analysis.
##* PCB concentrations in blood in children from LASERI and LUKAS projects.
##* PCB and dioxin concentrations in mother's milk from WHO dioxin studies
##* PCB and dioxin concentrations in blood of donors of Finnish Red Cross.

library(OpasnetUtils)
library(thlVerse)
library(thlConnect)
library(reshape2)

# PCBs of interest (the order of decreasing correlation with SUM-TEQ will be determined automatically later)
pcb9 <- c("PCB118","PCB138","PCB74","PCB156","PCB153","PCB99","PCB187","PCB170","PCB180")
teq3 <- c("PCDDF_TEQ", "PCB_TEQ", "Total_TEQ")

# Data from children in LUKAS and LASERI projects

boys <- read.csv("../Dioxdistboys/NuorDiox-Individual-Data_Jouni_2019-10-16.csv", header=TRUE, sep=";", dec=",")
boys <- boys[boys$Cohort != "", ] # Remove empty rows
boys$Id <- 1:nrow(boys)
colnames(boys)[colnames(boys)=="BirthYear"] <- "Birthyear"
boys$Birthyear <- as.integer(gsub(" ", "", as.character(boys$Birthyear)))
boys <- boys[!(grepl("TEQ", colnames(boys)) | colnames(boys) %in%
                 c("Code","CohortClass","SexClass","FirstBorn","FirstBornClass",
                   "AgeClass","Post_Office","VisitDate",""))]
boys$year <- boys$Age + boys$Birthyear
boys$Age <- factor(boys$Age, levels=c(7,8,9,10,17,18,19),
                   labels=c(rep("7-9",3),"10", rep("17-19",3)))
boys <- melt(boys, id.vars =c("Id","Cohort","Sex","Age","Center","Parity","Birthyear", "year"),
             variable.name="Compound",value.name="Result")
Var <- factor(substr((boys$Compound), nchar(as.character(boys$Compound)),1000)) # F fat or V volume-based?
boys <- boys[Var=="F",]
boys$Compound <- factor(substr(boys$Compound, 1, nchar(as.character(boys$Compound))-2))
boys <- boys[boys$Compound!="PCB9",] # SUM9PCB removed as redundant

# Data from mothers in WHO mothers' milk project 2000-2010

### THERE ARE MORE MOTHERS HERE THAN IN THE NEWER DATA! WHY?
#momsb <- read.csv("Äidinmaidot-primipara_2000-2010_yhteenveto_all.csv", header=TRUE, sep=";", dec=",", stringsAsFactors = FALSE)

moms <- read.csv("../Dioxdistboys/Äidinmaidot-kaikki_yhteenveto-SPSS-Siirto_2019-10-16.csv",
                 skip=5, header=TRUE, sep=";", dec=",", stringsAsFactors = FALSE)
moms <- moms[!is.na(moms$year),] # Remove empty rows from end
colnames(moms)[colnames(moms)=="id"] <- "Id"
colnames(moms)[colnames(moms)=="mage"] <- "Age"
colnames(moms)[colnames(moms)=="parity"] <- "Parity"
moms <- moms[!colnames(moms) %in% c("X","X.1","X.2","X.3","X.4","X.5",  # Remove empty columns
                                    "PCDD.F.2005TEQ","PCB.2005.TEQ","Total.TEQ","Sum_IndPCB","Sum_PCB9")] # Remove non-original columns

moms <- moms[!grepl("(BDE|PCN|DDD|DDE|DDT|PCB80|PCB122|PCB159)", colnames(moms))] # There is no data for PCB80,122,159

moms$Cohort <- "WHO mothers' milk"
moms$Center <- factor(ifelse(grepl("R",moms$Id),"Rovaniemi",ifelse(grepl("K",moms$Id),"Kuopio","Helsinki")))
moms$Birthyear <- moms$year - moms$Age
moms$Age <- factor(moms$Age)
moms$Sex <- "female"
moms$Parity <- as.character(ifelse(moms$Parity<4,moms$Parity,4))

moms <- melt(moms, id.vars =c("Id","Cohort","Sex","Age","Center","Parity","Birthyear","year"),
             variable.name="Compound",value.name="Result")
levels(moms$Compound) <- toupper(gsub("(^X|\\.)","", levels(moms$Compound))) # Remove unnecessary variation from congener names

# Add data from three pools of blood samples from the Finnish Red Cross

men <- read.csv("../Dioxdistboys/SPR-Diox-PCB.csv", header=TRUE, sep=";", dec=",", skip=2, stringsAsFactors = FALSE)
men <- men[men$Helsinki!="" , 1:4] # Remove empty rows and columns
colnames(men)[1] <- "Compound"
men$Compound <- toupper(trimws(gsub("(\\*|\\.|-|,)", "", men$Compound)))
men <- melt(men, id.vars = "Compound", variable.name = "Center", value.name = "Result")

# Adjust values below LOQ

LOQ <- grepl("<", men$Result) 
men$Result <- as.numeric(gsub(",",".",gsub("<","",men$Result)))
men$Result[LOQ] <- men$Result[LOQ]/10 # Take the lower bound but don't use zero
men$Id <- paste0("SPR_", men$Center)
men$year <- 2019
men$Cohort <- "SPR"
men$Sex <- "male"
men$Age <- "18-20"
men$Birthyear <- 2000
men$Parity <- "unknown"

## Combine children, women and donors to a single dataframe

pops <- rbind(
  cbind(boys, Subgroup = "Child"),
  cbind(moms, Subgroup = "Woman"),
  cbind(men, Subgroup = "Donor")
)
colnames(pops)[colnames(pops)=="year"] <- "Year"

pops$Parity[is.na(pops$Parity)] <- "unknown"

pops <- EvalOutput(Ovariable("pops",data=pops))

# Convert units from ng/g to pg/g where needed

conv <- paste0("PCB", c(18,28,33,47,49,51,52,60,66,74,80,99,101,105,110,114,118,122,123,
                        128,138,141,153,156,157,167,170,180,183,187,189,194,206,209)) # ng/g; convert to pg/g

result(pops)[pops$Compound %in% conv] <- result(pops)[pops$Compound %in% conv] * 1000 # ng/g --> pg/g

## Save data to R object file.

save(pops, file="../Dioxdistboys/pops_ovariable")

## Upload data to THL database

#dbi <- thlDbInfo("pubhealth", dbengine="postgresql")
#tst <- thlDbQuery(dbi,"SELECT * FROM envhealth.valvontayksikko_kunta")

con <- thlDbConnect("pubhealth",dbengine = "postgresql")

thlWriteTable(con, dat=pops@output, schema="envhealth", tablename="dioxdisthuman",
              owner="pubhealth_yte_admin",role="pubhealth_yte_user",grant="select")

thlSendStatement(con,"COMMENT ON COLUMN envhealth.dioxdisthuman.id IS 'identifier';")
thlSendStatement(con,"COMMENT ON COLUMN envhealth.dioxdisthuman.cohort IS 'which cohort the person belongs to';")
thlSendStatement(con,"COMMENT ON COLUMN envhealth.dioxdisthuman.popsresult IS 'dioxin or PCB concentration in pg/g';")

