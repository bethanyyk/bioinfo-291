#Code used to analyse the data associated with the paper, "Domestic dogs (Canis familiaris) distinguish conspecific emotional chemosignals."

#Authors: Alice Wang and Alexandra Horowitz

#Part 1: Habituation-discrimination paradigm
#Read data
install.packages("readxl")
library(readxl)
Part1 = 'Part_1.xlsx'

#Comparing investigation time across all trials
all <- read_excel(Part1, sheet = 'all')
#Testing habituation: comparing investigation time between first and third odor presentations during the habituation phase
wilcox.test(all$habit1, all$habit3, paired = TRUE)
#Testing discrimination: comparing investigation time between the novel odor and the fourth presentation of the habituated odor, which was presented simultaneously
wilcox.test(all$habit4, all$novel, paired = TRUE)
#Testing dishabituation: comparing investigation time between the novel odor and the third presentation of the habituated odor
wilcox.test(all$habit3, all$novel, paired = TRUE)

#Comparing investigation time in specific odor pairings
pairing <- read_excel(Part1, sheet = 'emotionalpairings')
#Testing discrimination between joy vs. baseline condition
wilcox.test(pairing$JBhabit4, pairing$JBnovel, paired = TRUE)
#Testing dishabituation between joy vs. baseline condition
wilcox.test(pairing$JBhabit3, pairing$JBnovel, paired = TRUE)
#Testing discrimination between stress vs. baseline
wilcox.test(pairing$SBhabit4, pairing$SBnovel, paired = TRUE)
#Testing dishabituation between stress vs. baseline
wilcox.test(pairing$SBhabit3, pairing$SBnovel, paired = TRUE)
#Testing discrimination between stress vs. joy
wilcox.test(pairing$SJhabit4, pairing$SJnovel, paired = TRUE)
#Testing dishabituation between stress vs. joy
wilcox.test(pairing$SJhabit3, pairing$SJnovel, paired = TRUE)

#Comparing investigation time accounting for presentation order within each emotional pairing
orders <- read_excel(Part1, sheet = 'specificorders')
#Testing discrimination when joy is the habituated odor and baseline is novel
wilcox.test(orders$Jhabit4_vB, orders$Bnovel_vJ, paired = TRUE)
#Testing dishabituation when joy is the habituated odor and baseline is novel
wilcox.test(orders$Jhabit3_vB, orders$Bnovel_vJ, paired = TRUE)
#Testing discrimination when baseline is the habituated odor and joy is novel
wilcox.test(orders$Bhabit4_vJ, orders$Jnovel_vB, paired = TRUE)
#Testing dishabituation when baseline is the habituated odor and joy is novel
wilcox.test(orders$Bhabit3_vJ, orders$Jnovel_vB, paired = TRUE)
#Testing discrimination when baseline is the habituated odor and stress is novel
wilcox.test(orders$Bhabit4_vS, orders$Snovel_vB, paired = TRUE)
#Testing dishabituation when baseline is the habituated odor and stress is novel
wilcox.test(orders$Bhabit3_vS, orders$Snovel_vB, paired = TRUE)
#Testing discrimination when stress is the habituated odor and baseline is novel
wilcox.test(orders$Shabit4_vB, orders$Bnovel_vS, paired = TRUE)
#Testing dishabituation when stress is the habituated odor and baseline is novel
wilcox.test(orders$Shabit3_vB, orders$Bnovel_vS, paired = TRUE)
#Testing discrimination when stress is the habituated odor and joy is novel
wilcox.test(orders$Shabit4_vJ, orders$Jnovel_vS, paired = TRUE)
#Testing dishabituation when stress is the habituated odor and joy is novel
wilcox.test(orders$Shabit3_vJ, orders$Jnovel_vS, paired = TRUE)
#Testing discrimination when joy is the habituated odor and stress is novel
wilcox.test(orders$Jhabit4_vS, orders$Snovel_vJ, paired = TRUE)
#Testing dishabituation when joy is the habituated odor and stress is novel
wilcox.test(orders$Jhabit3_vS, orders$Snovel_vJ, paired = TRUE)

#Testing dishabituation between the third and fourth presentation of habituated odors w/ change in presentation format
wilcox.test(all$habit3, all$habit4, paired = TRUE)
#Testing side bias or significant differences across data collection time periods in any emotional pairing
summary(lm(JBdiscrmdifference ~ JBside + Season, data = pairing))
summary(lm(JBdishabitdifference ~ JBside + Season, data = pairing))
summary(lm(SBdiscrmdifference ~ SBside + Season, data = pairing))
summary(lm(SBdishabitdifference ~ SBside + Season, data = pairing))
summary(lm(SJdiscrmdifference ~ SJside + Season, data = pairing))
summary(lm(SJdishabitdifference ~ SJside + Season, data = pairing))

#Part 2: Testing behavioral reactivity to conspecific emotional odors
#Read data
Part2 = 'Part_2.xlsx'
stressbehav <- read_excel(Part2, sheet = "stressbehav")
joybehav <- read_excel(Part2, sheet = "joybehav")
owner <- read_excel(Part2, sheet = "owner")
box <- read_excel(Part2, sheet = "box")
stranger <- read_excel(Part2, sheet = "stranger")
upright <- read_excel(Part2, sheet = "upright") 
sit <- read_excel(Part2, sheet = "sit")
prone <- read_excel(Part2, sheet = "prone")

#Testing whether dataset is normally distributed to determine whether to use ANOVA w/ Tukey’s pairwise post hoc comparisons (normally distributed) or Kruskal-Wallis test w/ Dunn’s post-hoc comparisons (not normally distributed) 
shapiro.test(stressbehav$Frequency)
shapiro.test(joybehav$Frequency)
shapiro.test(owner$Frequency)
shapiro.test(box$Frequency)
shapiro.test(stranger$Frequency)
shapiro.test(upright$Frequency)
shapiro.test(sit$Frequency)
shapiro.test(prone$Frequency)

#Comparing frequency of behaviors across different emotional odor conditions
install.packages("dunn.test")
library(dunn.test)
#Comparing frequency of stress behaviors across different emotional odor conditions
kruskal.test(Frequency~Condition, data = stressbehav)
dunn.test(stressbehav$Frequency, stressbehav$Condition, kw = TRUE, label = TRUE)
#Comparing frequency of joy behaviors across different emotional odor conditions
kruskal.test(Frequency~Condition, data = joybehav)
dunn.test(joybehav$Frequency, joybehav$Condition, kw = TRUE, label = TRUE)
#Comparing frequency of owner proximity across different emotional odor conditions
kruskal.test(Frequency~Condition, data = owner)
dunn.test(owner$Frequency, owner$Condition, kw = TRUE, label = TRUE)
#Comparing frequency of box proximity cross different emotional odor conditions
summary(aov(Frequency~Condition, data = box))
#Comparing frequency of stranger proximity across different emotional odor conditions
kruskal.test(Frequency~Condition, data = stranger)
dunn.test(stranger$Frequency, stranger$Condition, kw = TRUE, label = TRUE)
#Comparing frequency of upright posture across different emotional odor conditions
kruskal.test(Frequency~Condition, data = upright)
dunn.test(upright$Frequency, upright$Condition, kw = TRUE, label = TRUE)
#Comparing frequency of sit posture across different emotional odor conditions
kruskal.test(Frequency~Condition, data = sit)
#Comparing frequency of prone posture across different emotional odor conditions
kruskal.test(Frequency~Condition, data = prone)
dunn.test(prone$Frequency, prone$Condition, kw = TRUE, label = TRUE)

#The effect of emotional odor condition and sex on frequency of stress behaviors
glm(Frequency_MF~Condition_MF * Sex_MF, data = stressbehav)
summary(glm(Frequency_MF~Condition_MF * Sex_MF, data = stressbehav))
#The effect of emotional odor condition and sex on frequency of joy behaviors
glm(Frequency_MF~Condition_MF * Sex_MF, data = joybehav)
summary(glm(Frequency_MF~Condition_MF * Sex_MF, data = joybehav))
#The effect of emotional odor condition and sex on frequency of owner proximity 
glm(Frequency_MF~Condition_MF * Sex_MF, data = owner)
summary(glm(Frequency_MF~Condition_MF * Sex_MF, data = owner))
#The effect of emotional odor condition and sex on frequency of box proximity 
glm(Frequenc_MFy~Condition_MF * Sex_MF, data = box)
summary(glm(Frequency_MF~Condition_MF * Sex_MF, data = box))
#The effect of emotional odor condition and sex on frequency of stranger proximity 
glm(Frequency_MF~Condition_MF * Sex_MF, data = stranger)
summary(glm(Frequency_MF~Condition_MF * Sex_MF, data = stranger))
