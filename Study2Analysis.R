##########################################
#### LUTTRELL, PETTY, & BRINOL (2016) ####
############# STUDY 2: ALCOHOL ###########
##########################################

#### SET UP ###

# Set working directory
setwd("/Users/zhouleer/Desktop/研究生的文件夹/研一下学期26.3-26.7/R语言/再复现之/rd/")

### Load Functions and Import Data ###
if (!require("tidyverse")) install.packages("tidyverse") #For data manipulation
if (!require("broom")) install.packages("broom") #For presenting simple model outputs
if (!require("psych")) install.packages("psych") #For summary statistics
if (!require("Hmisc")) install.packages("Hmisc") #For correlation tables
if (!require("ggeffects")) install.packages("ggeffects") #For plotting results
if (!require("extrafont")) install.packages("extrafont") #For using alternative fonts in figures
if (!require("Cairo")) install.packages("Cairo") #For using rendering plots
select <- dplyr::select #Ensuring that select() calls the dplyr function

# Load function for displaying linear models
source("display_lm function.R")

# Plotting theme
## (Thanks to John Sakaluk for this particular setup: 
## https://sakaluk.wordpress.com/2015/08/27/6-make-it-pretty-plotting-2-way-interactions-with-ggplot2/#APA)
apatheme <- theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())

# Read in data
S2.full <- read_csv("Study2Data.csv") %>%
  mutate(Date1 = as.Date(Date1, "%m/%d/%Y"),
         Date2 = as.Date(Date2, "%m/%d/%Y"),
         Date3 = as.Date(Date3, "%m/%d/%Y"),
         delay.12 = Date2 - Date1,
         delay = Date3 - Date1,
         gender = factor(gender),
         amb = (neg + pos)/2 - abs(neg - pos),
         attT1 = (att1 + att2 + att3)/3,
         attT2 = (att1_T2 + att2_T2 + att3_T2)/3,
         attT3 = (att1_T3 + att2_T3 + att3_T3)/3,
         change = abs(attT3 - attT1), #Paper focuses on the longer one-year delay
         change.12 = abs(attT1 - attT2), #Footnotes effect when T2 is taken sooner
         samb = (subj1 + subj2 + subj3)/3,
         cert = select(., matches("cert[1-7]")) %>% rowMeans(na.rm = T),
         cert_clar = select(., matches("cert[1-4]")) %>% rowMeans(na.rm = T),
         cert_corr = select(., matches("cert[5-7]")) %>% rowMeans(na.rm = T),
         ex = abs(attT1 - 4))


#### METHOD ####

# Participants
## Time 1
nrow(S2.full) #T1 Sample Size
psych::describe(S2.full$age)
table(S2.full$gender, useNA = "ifany")

## Time 2
filter(S2.full, is.na(attT3) == F) %>% dplyr::summarize(T2Sample = n()) #Sample size for follow-up

## Predicting retention
S2.full <- mutate(S2.full, retention = ifelse(is.na(attT3) == F, 1, 0)) #0 = T1 Only; 1 = Both surveys
glm(retention ~ attT1 + cert + samb + amb + cert:amb, data = S2.full, family = "binomial") %>% summary()

S2 <- filter(S2.full, retention == 1) #Isolating only participants with T1 and T2 responses

# Measured Variables
## Attitudes (T1)
S2.full %>% select(att1, att2, att3) %>% psych::alpha() %>% with(total) %>% round(2)

## Attitudes (T3)
S2.full %>% select(matches("att[1-3]_T3")) %>% psych::alpha() %>% with(total) %>% round(2)

## Certainty
### Clarity
S2.full %>% select(matches("cert[1-4]")) %>% psych::alpha() %>% with(total) %>% round(2)

### Correctness
S2.full %>% select(matches("cert[5-7]")) %>% psych::alpha() %>% with(total) %>% round(2)

### Overall certainty
S2.full %>% select(matches("cert[1-7]")) %>% psych::alpha() %>% with(total) %>% round(2)

## Subjective Ambivalence
S2.full %>% select(matches("subj[1-7]")) %>% psych::alpha() %>% with(total) %>% round(2)


#### RESULTS ####

# Table 1 (Correlations & Summary Statistics)
select(S2, cert, amb, samb, change) %>% as.matrix() %>% rcorr()
select(S2, cert, amb, samb, change) %>% psych::describe() %>% select(mean, sd)

# Ambivalence-Certainty Correlations
cor.test(S2.full$amb, S2.full$cert)
cor.test(S2$amb, S2$cert)

# Stability analyses
## Main Effects
lm(change ~ ex + amb + cert, S2) %>% display.lm()

## Interaction
lm(change ~ ex + amb * cert, S2) %>% display.lm(rnd = 3)

## Simple Slopes
S2 <- S2 %>% mutate(ambLO = amb - mean(amb, na.rm = T) + sd(amb, na.rm = T),
                    ambHI = amb - mean(amb, na.rm = T) - sd(amb, na.rm = T))

lm(change ~ ex + ambLO * cert, S2) %>% display.lm() %>% filter(term == "cert")
lm(change ~ ex + ambHI * cert, S2) %>% display.lm() %>% filter(term == "cert")

## Simple Effects
S2 <- S2 %>% mutate(certLO = cert - mean(cert, na.rm = T) + sd(cert, na.rm = T),
                    certHI = cert - mean(cert, na.rm = T) - sd(cert, na.rm = T))

lm(change ~ ex + certLO * amb, S2) %>% display.lm() %>% filter(term == "amb")
lm(change ~ ex + certHI * amb, S2) %>% display.lm() %>% filter(term == "amb")

## Controlling for subjective ambivalence
lm(change ~ ex + samb + amb + cert, S2) %>% display.lm() #Main Effect Model
lm(change ~ ex + samb + amb * cert, S2) %>% display.lm() #Interaction Model
lm(change ~ ex + samb + ambLO * cert, S2) %>% display.lm() %>% filter(term == "cert") #Low Ambivalence Simple Effect
lm(change ~ ex + samb + ambHI * cert, S2) %>% display.lm() %>% filter(term == "cert") #High Ambivalence Simple Effect
lm(change ~ ex + samb + amb * certLO, S2) %>% display.lm() %>% filter(term == "amb") #Low Certainty Simple Slope
lm(change ~ ex + samb + amb * certHI, S2) %>% display.lm() %>% filter(term == "amb")#High Certainty Simple Slope

## Testing effects using subjective ambivalence instead of objective ambivalence
lm(change ~ ex + samb + cert, S2) %>% display.lm()
lm(change ~ ex + samb * cert, S2) %>% display.lm()


#### GRAPHS ####
# Plot as in paper
int.model <- lm(change ~ ex + amb * cert, S2)
plot.int <- ggpredict(int.model, terms = c("cert [meansd]", "amb [meansd]")) %>%
  as.data.frame() %>%
  mutate(amb = factor(group, labels = c("Low Ambivalence (-1 SD)", "Mean", "High Ambivalence (+1 SD)")),
         cert = factor(x, labels = c("(-1 SD)", "Mean", "(+1 SD)"))) %>%
  filter(amb != "Mean" & cert != "Mean") %>%
  mutate(amb = factor(amb))

ggplot(plot.int, aes(x = cert, y = predicted, group = amb, linetype = amb)) +
  geom_point(stat = 'summary', fun.y = sum) +
  stat_summary(fun.y = sum, geom = "line", size = 1) +
  apatheme + 
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Certainty", y = "Attitude Change", linetype = "") +
  theme(legend.position = 'top') +
  ylim(c(0,1.8))

# Full plot
plot.int.full <- ggpredict(int.model, terms = c("cert", "amb [meansd]")) %>%
  as.data.frame() %>%
  mutate(amb = factor(group, labels = c("Low Ambivalence (-1 SD)", "Mean", "High Ambivalence (+1 SD)"))) %>%
  filter(amb != "Mean")%>%
  mutate(amb = factor(amb))

ggplot(plot.int.full, aes(x = x, y = predicted, linetype = amb)) +
  geom_point(data = S2, aes(x = cert, y = change), alpha = .10,
             position = position_jitter(height = .05), 
             inherit.aes = F) +
  geom_line(size = 2, color = "black") +
  labs(x = "Certainty", y = "Attitude Change", linetype = "") + 
  apatheme +
  theme(legend.position = "top",
        legend.key.width = unit(.5, "inches"),
        text = element_text(family = "Book Antiqua", size = 16),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 12, b = 5)),
        axis.title.y = element_text(margin = margin(r = 8)),
        legend.text = element_text(size = 14)) +
  ylim(c(0,1.8))


#### FOOTNOTES ####

# Footnote 1: Accounting for delay
lm(change ~ delay + ex + amb * cert, S2) %>% display.lm()

# Footnote 2: Clarity vs. Correctness
## Clarity
lm(change ~ ex + cert_clar * amb, S2) %>% display.lm()

## Correctness
lm(change ~ ex + cert_corr * amb, S2) %>% display.lm()

# Footnote 6: Removing extremity covariate
lm(change ~ amb * cert, S2) %>% display.lm()


#### FUTURE FOOTNOTES! ####
# At one point, we had included this in the paper, but we ended up dropping it.
# We had originally fun this study with a 2-month delay with the same goal of showing the ambivalence x
# certainty interaction using a similar delay as in Study 1. We later collected attitudes again at a 
# one-year delay to test the longevity of the effect. Because it extends the effect beyond Study 1, we
# focused on this one-year delay in the paper. However, the effect was there at just a 2-month delay as well.

## Main Effects
lm(change.12 ~ ex + amb + cert, S2) %>% display.lm()

## Interaction
lm(change.12 ~ ex + amb * cert, S2) %>% display.lm(rnd = 3)
