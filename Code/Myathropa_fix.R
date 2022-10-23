MF_df_raw <- as.data.frame(readxl::read_excel("Data/MFlorea_raw.xlsx",
                                              sheet=1, .name_repair="minimal"))

MF_df_raw_fix <- as.data.frame(readxl::read_excel("Data/MFlorea_raw.xlsx",
                                                  sheet=1, .name_repair="minimal"))

names(MF_df_raw_fix) <- make.unique(names(MF_df_raw_fix))

MF_df_raw_fix_mig <- MF_df_raw_fix 

################################################################################
###################### MF239
################################################################################

######## General

# -> change 173 to NA just for LUX_137_1_30-07-21_B
MF_df_raw_fix[MF_df_raw_fix$`Sample Name`=="LUX_137_1_30-07-21_B","MF239.1"] <- NA

######## Migraine

# -> change 173 to NA
MF_df_raw_fix_mig[,"MF239"] <- replace(MF_df_raw_fix_mig[,"MF239"],
                                        which(MF_df_raw_fix_mig[,"MF239"]%in%c(173)),
                                        NA)

MF_df_raw_fix_mig[,"MF239.1"] <- replace(MF_df_raw_fix_mig[,"MF239.1"],
                                          which(MF_df_raw_fix_mig[,"MF239.1"]%in%c(173)),
                                          NA)

################################################################################
###################### MF265
################################################################################

######## General

# -> change 262 to NA
MF_df_raw_fix[,"MF265"] <- replace(MF_df_raw_fix[,"MF265"],
                                       which(MF_df_raw_fix[,"MF265"]%in%c(262)),
                                       NA)

MF_df_raw_fix[,"MF265.1"] <- replace(MF_df_raw_fix[,"MF265.1"],
                                   which(MF_df_raw_fix[,"MF265.1"]%in%c(262)),
                                   NA)

######## Migraine

# -> change 262 to NA
MF_df_raw_fix_mig[,"MF265"] <- replace(MF_df_raw_fix_mig[,"MF265"],
                                   which(MF_df_raw_fix_mig[,"MF265"]%in%c(262)),
                                   NA)

MF_df_raw_fix_mig[,"MF265.1"] <- replace(MF_df_raw_fix_mig[,"MF265.1"],
                                     which(MF_df_raw_fix_mig[,"MF265.1"]%in%c(262)),
                                     NA)

# -> change 244 to 245
MF_df_raw_fix_mig[,"MF265"] <- replace(MF_df_raw_fix_mig[,"MF265"],
                                           which(MF_df_raw_fix_mig[,"MF265"]%in%c(244)),
                                           245)

MF_df_raw_fix_mig[,"MF265.1"] <- replace(MF_df_raw_fix_mig[,"MF265.1"],
                                         which(MF_df_raw_fix_mig[,"MF265.1"]%in%c(244)),
                                         245)

# -> change 248 to 247
MF_df_raw_fix_mig[,"MF265"] <- replace(MF_df_raw_fix_mig[,"MF265"],
                                       which(MF_df_raw_fix_mig[,"MF265"]%in%c(248)),
                                       247)

MF_df_raw_fix_mig[,"MF265.1"] <- replace(MF_df_raw_fix_mig[,"MF265.1"],
                                         which(MF_df_raw_fix_mig[,"MF265.1"]%in%c(248)),
                                         247)

# -> change 250 and 252 to NA
MF_df_raw_fix_mig[,"MF265"] <- replace(MF_df_raw_fix_mig[,"MF265"],
                                       which(MF_df_raw_fix_mig[,"MF265"]%in%c(250,252)),
                                       NA)

MF_df_raw_fix_mig[,"MF265.1"] <- replace(MF_df_raw_fix_mig[,"MF265.1"],
                                         which(MF_df_raw_fix_mig[,"MF265.1"]%in%c(250,252)),
                                         NA)

################################################################################
###################### MF270
################################################################################

######## General

######## Migraine

# -> change 210 to NA
MF_df_raw_fix_mig[,"MF270"] <- replace(MF_df_raw_fix_mig[,"MF270"],
                                       which(MF_df_raw_fix_mig[,"MF270"]%in%c(210)),
                                       NA)

MF_df_raw_fix_mig[,"MF270.1"] <- replace(MF_df_raw_fix_mig[,"MF270.1"],
                                       which(MF_df_raw_fix_mig[,"MF270.1"]%in%c(210)),
                                       NA)

# -> change 208 to 209
MF_df_raw_fix_mig[,"MF270"] <- replace(MF_df_raw_fix_mig[,"MF270"],
                                       which(MF_df_raw_fix_mig[,"MF270"]%in%c(208)),
                                       209)

MF_df_raw_fix_mig[,"MF270.1"] <- replace(MF_df_raw_fix_mig[,"MF270.1"],
                                         which(MF_df_raw_fix_mig[,"MF270.1"]%in%c(208)),
                                         209)

################################################################################
###################### MF59
################################################################################

######## General

# -> change 178 to NA
MF_df_raw_fix[,"MF59"] <- replace(MF_df_raw_fix[,"MF59"],
                                      which(MF_df_raw_fix[,"MF59"]%in%c(178)),
                                      NA)

MF_df_raw_fix[,"MF59.1"] <- replace(MF_df_raw_fix[,"MF59.1"],
                                        which(MF_df_raw_fix[,"MF59.1"]%in%c(178)),
                                        NA)

######## Migraine

# -> change 178 to NA
MF_df_raw_fix_mig[,"MF59"] <- replace(MF_df_raw_fix_mig[,"MF59"],
                                       which(MF_df_raw_fix_mig[,"MF59"]%in%c(178)),
                                       NA)

MF_df_raw_fix_mig[,"MF59.1"] <- replace(MF_df_raw_fix_mig[,"MF59.1"],
                                         which(MF_df_raw_fix_mig[,"MF59.1"]%in%c(178)),
                                         NA)

# -> change 180 to NA
MF_df_raw_fix_mig[,"MF59"] <- replace(MF_df_raw_fix_mig[,"MF59"],
                                       which(MF_df_raw_fix_mig[,"MF59"]%in%c(180)),
                                       NA)

MF_df_raw_fix_mig[,"MF59.1"] <- replace(MF_df_raw_fix_mig[,"MF59.1"],
                                         which(MF_df_raw_fix_mig[,"MF59.1"]%in%c(180)),
                                         NA)

# -> change 175 to NA
MF_df_raw_fix_mig[,"MF59"] <- replace(MF_df_raw_fix_mig[,"MF59"],
                                      which(MF_df_raw_fix_mig[,"MF59"]%in%c(175)),
                                      NA)

MF_df_raw_fix_mig[,"MF59.1"] <- replace(MF_df_raw_fix_mig[,"MF59.1"],
                                        which(MF_df_raw_fix_mig[,"MF59.1"]%in%c(175)),
                                        NA)

# -> change 156 to 155
MF_df_raw_fix_mig[,"MF59"] <- replace(MF_df_raw_fix_mig[,"MF59"],
                                      which(MF_df_raw_fix_mig[,"MF59"]%in%c(156)),
                                      155)

MF_df_raw_fix_mig[,"MF59.1"] <- replace(MF_df_raw_fix_mig[,"MF59.1"],
                                        which(MF_df_raw_fix_mig[,"MF59.1"]%in%c(156)),
                                        155)

# -> change 159 to 158
MF_df_raw_fix_mig[,"MF59"] <- replace(MF_df_raw_fix_mig[,"MF59"],
                                      which(MF_df_raw_fix_mig[,"MF59"]%in%c(159)),
                                      158)

MF_df_raw_fix_mig[,"MF59.1"] <- replace(MF_df_raw_fix_mig[,"MF59.1"],
                                        which(MF_df_raw_fix_mig[,"MF59.1"]%in%c(159)),
                                        158)

# -> change 168 to NA and 169

MF_df_raw_fix_mig[MF_df_raw_fix_mig$`Sample Name`=="CO_304_1_01-09-21","MF59"] <- 169
MF_df_raw_fix_mig[MF_df_raw_fix_mig$`Sample Name`=="CO_304_1_01-09-21","MF59.1"] <- 169
MF_df_raw_fix_mig[MF_df_raw_fix_mig$`Sample Name`=="CO_184_2_06-08-21","MF59"] <- NA

################################################################################
###################### MF197
################################################################################

######## General

######## Migraine

# -> change 139 to 140
MF_df_raw_fix_mig[,"MF197"] <- replace(MF_df_raw_fix_mig[,"MF197"],
                                      which(MF_df_raw_fix_mig[,"MF197"]%in%c(139)),
                                      140)

MF_df_raw_fix_mig[,"MF197.1"] <- replace(MF_df_raw_fix_mig[,"MF197.1"],
                                        which(MF_df_raw_fix_mig[,"MF197.1"]%in%c(139)),
                                        140)

# -> change 141 to NA
MF_df_raw_fix_mig[,"MF197"] <- replace(MF_df_raw_fix_mig[,"MF197"],
                                       which(MF_df_raw_fix_mig[,"MF197"]%in%c(141)),
                                       NA)

MF_df_raw_fix_mig[,"MF197.1"] <- replace(MF_df_raw_fix_mig[,"MF197.1"],
                                         which(MF_df_raw_fix_mig[,"MF197.1"]%in%c(141)),
                                         NA)

# -> change 143 to NA
MF_df_raw_fix_mig[,"MF197"] <- replace(MF_df_raw_fix_mig[,"MF197"],
                                       which(MF_df_raw_fix_mig[,"MF197"]%in%c(143)),
                                       NA)

MF_df_raw_fix_mig[,"MF197.1"] <- replace(MF_df_raw_fix_mig[,"MF197.1"],
                                         which(MF_df_raw_fix_mig[,"MF197.1"]%in%c(143)),
                                         NA)

# -> change 149 to NA
MF_df_raw_fix_mig[,"MF197"] <- replace(MF_df_raw_fix_mig[,"MF197"],
                                       which(MF_df_raw_fix_mig[,"MF197"]%in%c(149)),
                                       NA)

MF_df_raw_fix_mig[,"MF197.1"] <- replace(MF_df_raw_fix_mig[,"MF197.1"],
                                         which(MF_df_raw_fix_mig[,"MF197.1"]%in%c(149)),
                                         NA)

################################################################################
###################### MF36
################################################################################

######## General

######## Migraine

# -> change 89 to 90

MF_df_raw_fix_mig[,"MF36"] <- replace(MF_df_raw_fix_mig[,"MF36"],
                                       which(MF_df_raw_fix_mig[,"MF36"]%in%c(89)),
                                       90)

MF_df_raw_fix_mig[,"MF36.1"] <- replace(MF_df_raw_fix_mig[,"MF36.1"],
                                         which(MF_df_raw_fix_mig[,"MF36.1"]%in%c(89)),
                                         90)
                                        
# -> change 102 to NA

MF_df_raw_fix_mig[,"MF36"] <- replace(MF_df_raw_fix_mig[,"MF36"],
                                      which(MF_df_raw_fix_mig[,"MF36"]%in%c(102)),
                                      NA)

MF_df_raw_fix_mig[,"MF36.1"] <- replace(MF_df_raw_fix_mig[,"MF36.1"],
                                        which(MF_df_raw_fix_mig[,"MF36.1"]%in%c(102)),
                                              NA)

# -> change 109 to 108

MF_df_raw_fix_mig[,"MF36"] <- replace(MF_df_raw_fix_mig[,"MF36"],
                                      which(MF_df_raw_fix_mig[,"MF36"]%in%c(109)),
                                      108)

MF_df_raw_fix_mig[,"MF36.1"] <- replace(MF_df_raw_fix_mig[,"MF36.1"],
                                        which(MF_df_raw_fix_mig[,"MF36.1"]%in%c(109)),
                                              108)                                        

# -> change 103 to 102

MF_df_raw_fix_mig[,"MF36"] <- replace(MF_df_raw_fix_mig[,"MF36"],
                                       which(MF_df_raw_fix_mig[,"MF36"]%in%c(103)),
                                       102)

MF_df_raw_fix_mig[,"MF36.1"] <- replace(MF_df_raw_fix_mig[,"MF36.1"],
                                        which(MF_df_raw_fix_mig[,"MF36.1"]%in%c(103)),
                                              102)
                                        
# -> change 98 to NA

MF_df_raw_fix_mig[,"MF36"] <- replace(MF_df_raw_fix_mig[,"MF36"],
                                      which(MF_df_raw_fix_mig[,"MF36"]%in%c(98)),
                                      NA)

MF_df_raw_fix_mig[,"MF36.1"] <- replace(MF_df_raw_fix_mig[,"MF36.1"],
                                        which(MF_df_raw_fix_mig[,"MF36.1"]%in%c(98)),
                                              NA)
                                        
# -> change 92 to NA

MF_df_raw_fix_mig[,"MF36"] <- replace(MF_df_raw_fix_mig[,"MF36"],
                                      which(MF_df_raw_fix_mig[,"MF36"]%in%c(92)),
                                      NA)

MF_df_raw_fix_mig[,"MF36.1"] <- replace(MF_df_raw_fix_mig[,"MF36.1"],
                                        which(MF_df_raw_fix_mig[,"MF36.1"]%in%c(92)),
                                              NA)

################################################################################
###################### MF492
################################################################################

######## General

######## Migraine                                        
                               
# -> change 151 to 150

MF_df_raw_fix_mig[,"MF492"] <- replace(MF_df_raw_fix_mig[,"MF492"],
                                      which(MF_df_raw_fix_mig[,"MF492"]%in%c(151)),
                                      150)

MF_df_raw_fix_mig[,"MF492.1"] <- replace(MF_df_raw_fix_mig[,"MF492.1"],
                                        which(MF_df_raw_fix_mig[,"MF492.1"]%in%c(151)),
                                              150)
                                        
# -> change 152 to 151

MF_df_raw_fix_mig[,"MF492"] <- replace(MF_df_raw_fix_mig[,"MF492"],
                                       which(MF_df_raw_fix_mig[,"MF492"]%in%c(152)),
                                       151)

MF_df_raw_fix_mig[,"MF492.1"] <- replace(MF_df_raw_fix_mig[,"MF492.1"],
                                         which(MF_df_raw_fix_mig[,"MF492.1"]%in%c(152)),
                                               151)
                                         
# -> change 153 to 152

MF_df_raw_fix_mig[,"MF492"] <- replace(MF_df_raw_fix_mig[,"MF492"],
                                       which(MF_df_raw_fix_mig[,"MF492"]%in%c(153)),
                                       152)

MF_df_raw_fix_mig[,"MF492.1"] <- replace(MF_df_raw_fix_mig[,"MF492.1"],
                                         which(MF_df_raw_fix_mig[,"MF492.1"]%in%c(153)),
                                         152)
                                         
# -> change 147 to NA

MF_df_raw_fix_mig[,"MF492"] <- replace(MF_df_raw_fix_mig[,"MF492"],
                                       which(MF_df_raw_fix_mig[,"MF492"]%in%c(147)),
                                       NA)

MF_df_raw_fix_mig[,"MF492.1"] <- replace(MF_df_raw_fix_mig[,"MF492.1"],
                                         which(MF_df_raw_fix_mig[,"MF492.1"]%in%c(147)),
                                               NA)
                                         
# -> change 145 to NA

MF_df_raw_fix_mig[,"MF492"] <- replace(MF_df_raw_fix_mig[,"MF492"],
                                       which(MF_df_raw_fix_mig[,"MF492"]%in%c(145)),
                                       NA)

MF_df_raw_fix_mig[,"MF492.1"] <- replace(MF_df_raw_fix_mig[,"MF492.1"],
                                         which(MF_df_raw_fix_mig[,"MF492.1"]%in%c(145)),
                                               NA)
                                         
# -> change 141 to 140

MF_df_raw_fix_mig[,"MF492"] <- replace(MF_df_raw_fix_mig[,"MF492"],
                                       which(MF_df_raw_fix_mig[,"MF492"]%in%c(141)),
                                       140)

MF_df_raw_fix_mig[,"MF492.1"] <- replace(MF_df_raw_fix_mig[,"MF492.1"],
                                         which(MF_df_raw_fix_mig[,"MF492.1"]%in%c(141)),
                                               140)
                                         
# -> change 139 to 138

MF_df_raw_fix_mig[,"MF492"] <- replace(MF_df_raw_fix_mig[,"MF492"],
                                       which(MF_df_raw_fix_mig[,"MF492"]%in%c(139)),
                                       138)

MF_df_raw_fix_mig[,"MF492.1"] <- replace(MF_df_raw_fix_mig[,"MF492.1"],
                                         which(MF_df_raw_fix_mig[,"MF492.1"]%in%c(139)),
                                               138)                                          

################################################################################
###################### MF103
################################################################################

######## General

######## Migraine

# -> change 150 and 152 to NA

MF_df_raw_fix_mig[,"MF103"] <- replace(MF_df_raw_fix_mig[,"MF103"],
                                       which(MF_df_raw_fix_mig[,"MF103"]%in%c(150, 152)),
                                       NA)

MF_df_raw_fix_mig[,"MF103.1"] <- replace(MF_df_raw_fix_mig[,"MF103.1"],
                                         which(MF_df_raw_fix_mig[,"MF103.1"]%in%c(150, 152)),
                                               NA)
                                         
################################################################################
###################### MF26
################################################################################

######## General

######## Migraine

# -> change 147 to 148

MF_df_raw_fix_mig[,"MF26"] <- replace(MF_df_raw_fix_mig[,"MF26"],
                                       which(MF_df_raw_fix_mig[,"MF26"]%in%c(147)),
                                       148)

MF_df_raw_fix_mig[,"MF26.1"] <- replace(MF_df_raw_fix_mig[,"MF26.1"],
                                         which(MF_df_raw_fix_mig[,"MF26.1"]%in%c(147)),
                                               148)                                          

# -> change 129 to 128

MF_df_raw_fix_mig[,"MF26"] <- replace(MF_df_raw_fix_mig[,"MF26"],
                                      which(MF_df_raw_fix_mig[,"MF26"]%in%c(129)),
                                      128)

MF_df_raw_fix_mig[,"MF26.1"] <- replace(MF_df_raw_fix_mig[,"MF26.1"],
                                        which(MF_df_raw_fix_mig[,"MF26.1"]%in%c(129)),
                                        128)
                                        
# -> change 141 to 140

MF_df_raw_fix_mig[,"MF26"] <- replace(MF_df_raw_fix_mig[,"MF26"],
                                      which(MF_df_raw_fix_mig[,"MF26"]%in%c(141)),
                                      140)

MF_df_raw_fix_mig[,"MF26.1"] <- replace(MF_df_raw_fix_mig[,"MF26.1"],
                                        which(MF_df_raw_fix_mig[,"MF26.1"]%in%c(141)),
                                              140)
                                        
# -> change 139 and 143 to NA

MF_df_raw_fix_mig[,"MF26"] <- replace(MF_df_raw_fix_mig[,"MF26"],
                                      which(MF_df_raw_fix_mig[,"MF26"]%in%c(139,143)),
                                      NA)

MF_df_raw_fix_mig[,"MF26.1"] <- replace(MF_df_raw_fix_mig[,"MF26.1"],
                                        which(MF_df_raw_fix_mig[,"MF26.1"]%in%c(139,143)),
                                              NA)
################################################################################
###################### MF263
################################################################################

######## General

######## Migraine

# -> change 87 to 88

MF_df_raw_fix_mig[,"MF263"] <- replace(MF_df_raw_fix_mig[,"MF263"],
                                      which(MF_df_raw_fix_mig[,"MF263"]%in%c(87)),
                                      88)

MF_df_raw_fix_mig[,"MF263.1"] <- replace(MF_df_raw_fix_mig[,"MF263.1"],
                                        which(MF_df_raw_fix_mig[,"MF263.1"]%in%c(87)),
                                              88)
                    
# -> change 97 to NA

MF_df_raw_fix_mig[,"MF263"] <- replace(MF_df_raw_fix_mig[,"MF263"],
                                       which(MF_df_raw_fix_mig[,"MF263"]%in%c(97)),
                                       NA)

MF_df_raw_fix_mig[,"MF263.1"] <- replace(MF_df_raw_fix_mig[,"MF263.1"],
                                         which(MF_df_raw_fix_mig[,"MF263.1"]%in%c(97)),
                                         NA)

################################################################################
###################### MF269
################################################################################

######## General

######## Migraine

# -> change 245 to 246

MF_df_raw_fix_mig[,"MF269"] <- replace(MF_df_raw_fix_mig[,"MF269"],
                                       which(MF_df_raw_fix_mig[,"MF269"]%in%c(245)),
                                       246)

MF_df_raw_fix_mig[,"MF269.1"] <- replace(MF_df_raw_fix_mig[,"MF269.1"],
                                         which(MF_df_raw_fix_mig[,"MF269.1"]%in%c(245)),
                                               246)
                                         
################################################################################
###################### MF323
################################################################################

######## General

######## Migraine

# -> change 166 to 167

MF_df_raw_fix_mig[,"MF323"] <- replace(MF_df_raw_fix_mig[,"MF323"],
                                       which(MF_df_raw_fix_mig[,"MF323"]%in%c(166)),
                                       167)

MF_df_raw_fix_mig[,"MF323.1"] <- replace(MF_df_raw_fix_mig[,"MF323.1"],
                                         which(MF_df_raw_fix_mig[,"MF323.1"]%in%c(166)),
                                               167)
                                         
# -> change 178 to 179

MF_df_raw_fix_mig[,"MF323"] <- replace(MF_df_raw_fix_mig[,"MF323"],
                                       which(MF_df_raw_fix_mig[,"MF323"]%in%c(178)),
                                       179)

MF_df_raw_fix_mig[,"MF323.1"] <- replace(MF_df_raw_fix_mig[,"MF323.1"],
                                         which(MF_df_raw_fix_mig[,"MF323.1"]%in%c(178)),
                                               179)
                                         
# -> change 175 to NA

MF_df_raw_fix_mig[,"MF323"] <- replace(MF_df_raw_fix_mig[,"MF323"],
                                       which(MF_df_raw_fix_mig[,"MF323"]%in%c(175)),
                                       NA)

MF_df_raw_fix_mig[,"MF323.1"] <- replace(MF_df_raw_fix_mig[,"MF323.1"],
                                         which(MF_df_raw_fix_mig[,"MF323.1"]%in%c(175)),
                                               NA)

################################################################################
###################### MF457
################################################################################

######## General

######## Migraine

# -> change 176 to 175

MF_df_raw_fix_mig[,"MF457"] <- replace(MF_df_raw_fix_mig[,"MF457"],
                                       which(MF_df_raw_fix_mig[,"MF457"]%in%c(176)),
                                       175)

MF_df_raw_fix_mig[,"MF457.1"] <- replace(MF_df_raw_fix_mig[,"MF457.1"],
                                         which(MF_df_raw_fix_mig[,"MF457.1"]%in%c(176)),
                                         175)

# -> change 192 to 191

MF_df_raw_fix_mig[,"MF457"] <- replace(MF_df_raw_fix_mig[,"MF457"],
                                       which(MF_df_raw_fix_mig[,"MF457"]%in%c(192)),
                                       191)

MF_df_raw_fix_mig[,"MF457.1"] <- replace(MF_df_raw_fix_mig[,"MF457.1"],
                                         which(MF_df_raw_fix_mig[,"MF457.1"]%in%c(192)),
                                         191)

# -> change 170 to NA

MF_df_raw_fix_mig[,"MF457"] <- replace(MF_df_raw_fix_mig[,"MF457"],
                                       which(MF_df_raw_fix_mig[,"MF457"]%in%c(170)),
                                       NA)

MF_df_raw_fix_mig[,"MF457.1"] <- replace(MF_df_raw_fix_mig[,"MF457.1"],
                                         which(MF_df_raw_fix_mig[,"MF457.1"]%in%c(170)),
                                         NA)

################################################################################
###################### MF491
################################################################################

######## General

######## Migraine

# -> change 225 to 226

MF_df_raw_fix_mig[,"MF491"] <- replace(MF_df_raw_fix_mig[,"MF491"],
                                       which(MF_df_raw_fix_mig[,"MF491"]%in%c(225)),
                                       226)

MF_df_raw_fix_mig[,"MF491.1"] <- replace(MF_df_raw_fix_mig[,"MF491.1"],
                                         which(MF_df_raw_fix_mig[,"MF491.1"]%in%c(225)),
                                         226)

# -> change 200 to 199

MF_df_raw_fix_mig[,"MF491"] <- replace(MF_df_raw_fix_mig[,"MF491"],
                                       which(MF_df_raw_fix_mig[,"MF491"]%in%c(200)),
                                       199)

MF_df_raw_fix_mig[,"MF491.1"] <- replace(MF_df_raw_fix_mig[,"MF491.1"],
                                         which(MF_df_raw_fix_mig[,"MF491.1"]%in%c(200)),
                                         199)

################################################################################
###################### MF70
################################################################################

######## General

# -> change 169 to NA for Cologne individuals

MF_df_raw_fix[MF_df_raw_fix$`Sample Name`=="CO_143_1_01-08-21_B","MF70"] <- NA
MF_df_raw_fix[MF_df_raw_fix$`Sample Name`=="CO_146_1_04-09-21_A","MF70"] <- NA
MF_df_raw_fix[MF_df_raw_fix$`Sample Name`=="CO_266_1_26-07-21_A","MF70"] <- NA

######## Migraine

# -> change 169 to NA

MF_df_raw_fix_mig[,"MF70"] <- replace(MF_df_raw_fix_mig[,"MF70"],
                                       which(MF_df_raw_fix_mig[,"MF70"]%in%c(169)),
                                       NA)

MF_df_raw_fix_mig[,"MF70.1"] <- replace(MF_df_raw_fix_mig[,"MF70.1"],
                                         which(MF_df_raw_fix_mig[,"MF70.1"]%in%c(169)),
                                         NA)

################################################################################
# Get proper names back
names(MF_df_raw_fix_mig) <- names(MF_df_raw_fix)

################################################################################
# GET RID OF LOCI WITH IRREPARABLE WEIRD STEP SIZES (base disruption etc)

MF_df_raw_fix_mig[,"MF130"] <- NULL
MF_df_raw_fix_mig[,"MF130.1"] <- NULL
MF_df_raw_fix_mig[,"MF419"] <- NULL
MF_df_raw_fix_mig[,"MF419.1"] <- NULL
MF_df_raw_fix_mig[,"MF261"] <- NULL
MF_df_raw_fix_mig[,"MF261.1"] <- NULL
MF_df_raw_fix_mig[,"MF28"] <- NULL
MF_df_raw_fix_mig[,"MF28.1"] <- NULL