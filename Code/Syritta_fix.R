SP_df_raw <- as.data.frame(readxl::read_excel("Data/SPipiens_raw.xlsx",
                                                  sheet=1, .name_repair="minimal"))

SP_df_raw_fix <- as.data.frame(readxl::read_excel("Data/SPipiens_raw.xlsx",
                                              sheet=1, .name_repair="minimal"))

names(SP_df_raw_fix) <- make.unique(names(SP_df_raw_fix))

SP_df_raw_fix_mig <- SP_df_raw_fix 

################################################################################
###################### Spp010
################################################################################

######## General

# -> change 150 and 151 to 153
SP_df_raw_fix[,"Spp010"] <- replace(SP_df_raw_fix[,"Spp010"],
                                    which(SP_df_raw_fix[,"Spp010"]%in%c(150,151)),
                                    153)

SP_df_raw_fix[,"Spp010.1"] <- replace(SP_df_raw_fix[,"Spp010.1"],
                                    which(SP_df_raw_fix[,"Spp010.1"]%in%c(150,151)),
                                    153)
######## Migraine

# -> change 150 and 151 to 153
SP_df_raw_fix_mig[,"Spp010"] <- replace(SP_df_raw_fix_mig[,"Spp010"],
                                    which(SP_df_raw_fix_mig[,"Spp010"]%in%c(150,151)),
                                    153)

SP_df_raw_fix_mig[,"Spp010.1"] <- replace(SP_df_raw_fix_mig[,"Spp010.1"],
                                      which(SP_df_raw_fix_mig[,"Spp010.1"]%in%c(150,151)),
                                      153)

# -> change 160 to NA
SP_df_raw_fix_mig[,"Spp010"] <- replace(SP_df_raw_fix_mig[,"Spp010"],
                                    which(SP_df_raw_fix_mig[,"Spp010"]%in%c(160)),
                                    NA)

SP_df_raw_fix_mig[,"Spp010.1"] <- replace(SP_df_raw_fix_mig[,"Spp010.1"],
                                    which(SP_df_raw_fix_mig[,"Spp010.1"]%in%c(160)),
                                    NA)

# -> change 166 to 165, 170 to 169, and 174 to 173
SP_df_raw_fix_mig[,"Spp010"] <- replace(SP_df_raw_fix_mig[,"Spp010"],
                                        which(SP_df_raw_fix_mig[,"Spp010"]%in%c(166)),
                                        165)

SP_df_raw_fix_mig[,"Spp010.1"] <- replace(SP_df_raw_fix_mig[,"Spp010.1"],
                                          which(SP_df_raw_fix_mig[,"Spp010.1"]%in%c(166)),
                                          165)

SP_df_raw_fix_mig[,"Spp010"] <- replace(SP_df_raw_fix_mig[,"Spp010"],
                                        which(SP_df_raw_fix_mig[,"Spp010"]%in%c(170)),
                                        169)

SP_df_raw_fix_mig[,"Spp010.1"] <- replace(SP_df_raw_fix_mig[,"Spp010.1"],
                                          which(SP_df_raw_fix_mig[,"Spp010.1"]%in%c(170)),
                                          169)

SP_df_raw_fix_mig[,"Spp010"] <- replace(SP_df_raw_fix_mig[,"Spp010"],
                                        which(SP_df_raw_fix_mig[,"Spp010"]%in%c(174)),
                                        173)

SP_df_raw_fix_mig[,"Spp010.1"] <- replace(SP_df_raw_fix_mig[,"Spp010.1"],
                                          which(SP_df_raw_fix_mig[,"Spp010.1"]%in%c(174)),
                                          173)

################################################################################
###################### Spp053
################################################################################

######## General

# no changes

######## Migraine

# -> change 137 to 136
SP_df_raw_fix_mig[,"Spp053"] <- replace(SP_df_raw_fix_mig[,"Spp053"],
                                        which(SP_df_raw_fix_mig[,"Spp053"]%in%c(137)),
                                        136)

SP_df_raw_fix_mig[,"Spp053.1"] <- replace(SP_df_raw_fix_mig[,"Spp053.1"],
                                          which(SP_df_raw_fix_mig[,"Spp053.1"]%in%c(137)),
                                          136)

# -> change 147 to 146
SP_df_raw_fix_mig[,"Spp053"] <- replace(SP_df_raw_fix_mig[,"Spp053"],
                                        which(SP_df_raw_fix_mig[,"Spp053"]%in%c(147)),
                                        146)

SP_df_raw_fix_mig[,"Spp053.1"] <- replace(SP_df_raw_fix_mig[,"Spp053.1"],
                                          which(SP_df_raw_fix_mig[,"Spp053.1"]%in%c(147)),
                                          146)

# -> change 155 to NA
SP_df_raw_fix_mig[,"Spp053"] <- replace(SP_df_raw_fix_mig[,"Spp053"],
                                        which(SP_df_raw_fix_mig[,"Spp053"]%in%c(155)),
                                        NA)

SP_df_raw_fix_mig[,"Spp053.1"] <- replace(SP_df_raw_fix_mig[,"Spp053.1"],
                                          which(SP_df_raw_fix_mig[,"Spp053.1"]%in%c(155)),
                                          NA)

# -> change 159 to NA
SP_df_raw_fix_mig[,"Spp053"] <- replace(SP_df_raw_fix_mig[,"Spp053"],
                                        which(SP_df_raw_fix_mig[,"Spp053"]%in%c(159)),
                                        NA)

SP_df_raw_fix_mig[,"Spp053.1"] <- replace(SP_df_raw_fix_mig[,"Spp053.1"],
                                          which(SP_df_raw_fix_mig[,"Spp053.1"]%in%c(159)),
                                          NA)

################################################################################
###################### Spp080
################################################################################

######## General

# -> change 163 to 164
SP_df_raw_fix[,"Spp080"] <- replace(SP_df_raw_fix[,"Spp080"],
                                    which(SP_df_raw_fix[,"Spp080"]%in%c(163)),
                                    164)

SP_df_raw_fix[,"Spp080.1"] <- replace(SP_df_raw_fix[,"Spp080.1"],
                                      which(SP_df_raw_fix[,"Spp080.1"]%in%c(163)),
                                      164)

######## Migraine

# -> change 163 to 164
SP_df_raw_fix_mig[,"Spp080"] <- replace(SP_df_raw_fix_mig[,"Spp080"],
                                    which(SP_df_raw_fix_mig[,"Spp080"]%in%c(163)),
                                    164)

SP_df_raw_fix_mig[,"Spp080.1"] <- replace(SP_df_raw_fix_mig[,"Spp080.1"],
                                      which(SP_df_raw_fix_mig[,"Spp080.1"]%in%c(163)),
                                      164)

# -> change 147 to 146
SP_df_raw_fix_mig[,"Spp080"] <- replace(SP_df_raw_fix_mig[,"Spp080"],
                                        which(SP_df_raw_fix_mig[,"Spp080"]%in%c(147)),
                                        146)

SP_df_raw_fix_mig[,"Spp080.1"] <- replace(SP_df_raw_fix_mig[,"Spp080.1"],
                                          which(SP_df_raw_fix_mig[,"Spp080.1"]%in%c(147)),
                                          146)

# -> change 165 to NA
SP_df_raw_fix_mig[,"Spp080"] <- replace(SP_df_raw_fix_mig[,"Spp080"],
                                        which(SP_df_raw_fix_mig[,"Spp080"]%in%c(165)),
                                        NA)

SP_df_raw_fix_mig[,"Spp080.1"] <- replace(SP_df_raw_fix_mig[,"Spp080.1"],
                                          which(SP_df_raw_fix_mig[,"Spp080.1"]%in%c(165)),
                                          NA)

# -> change 167 to NA
SP_df_raw_fix_mig[,"Spp080"] <- replace(SP_df_raw_fix_mig[,"Spp080"],
                                        which(SP_df_raw_fix_mig[,"Spp080"]%in%c(167)),
                                        NA)

SP_df_raw_fix_mig[,"Spp080.1"] <- replace(SP_df_raw_fix_mig[,"Spp080.1"],
                                          which(SP_df_raw_fix_mig[,"Spp080.1"]%in%c(167)),
                                          NA)

################################################################################
###################### Spp142
################################################################################

######## General

######## Migraine

# -> change 75 to 76
SP_df_raw_fix_mig[,"Spp142"] <- replace(SP_df_raw_fix_mig[,"Spp142"],
                                        which(SP_df_raw_fix_mig[,"Spp142"]%in%c(75)),
                                        76)

SP_df_raw_fix_mig[,"Spp142.1"] <- replace(SP_df_raw_fix_mig[,"Spp142.1"],
                                          which(SP_df_raw_fix_mig[,"Spp142.1"]%in%c(75)),
                                          76)

# -> change 91 to NA
SP_df_raw_fix_mig[,"Spp142"] <- replace(SP_df_raw_fix_mig[,"Spp142"],
                                        which(SP_df_raw_fix_mig[,"Spp142"]%in%c(91)),
                                        NA)

SP_df_raw_fix_mig[,"Spp142.1"] <- replace(SP_df_raw_fix_mig[,"Spp142.1"],
                                          which(SP_df_raw_fix_mig[,"Spp142.1"]%in%c(91)),
                                          NA)

# -> change 93 to NA
SP_df_raw_fix_mig[,"Spp142"] <- replace(SP_df_raw_fix_mig[,"Spp142"],
                                        which(SP_df_raw_fix_mig[,"Spp142"]%in%c(93)),
                                        NA)

SP_df_raw_fix_mig[,"Spp142.1"] <- replace(SP_df_raw_fix_mig[,"Spp142.1"],
                                          which(SP_df_raw_fix_mig[,"Spp142.1"]%in%c(93)),
                                          NA)


# -> change 97 to NA
SP_df_raw_fix_mig[,"Spp142"] <- replace(SP_df_raw_fix_mig[,"Spp142"],
                                        which(SP_df_raw_fix_mig[,"Spp142"]%in%c(97)),
                                        NA)

SP_df_raw_fix_mig[,"Spp142.1"] <- replace(SP_df_raw_fix_mig[,"Spp142.1"],
                                          which(SP_df_raw_fix_mig[,"Spp142.1"]%in%c(97)),
                                          NA)

# -> change change 99 to 100
SP_df_raw_fix_mig[,"Spp142"] <- replace(SP_df_raw_fix_mig[,"Spp142"],
                                        which(SP_df_raw_fix_mig[,"Spp142"]%in%c(99)),
                                        100)

SP_df_raw_fix_mig[,"Spp142.1"] <- replace(SP_df_raw_fix_mig[,"Spp142.1"],
                                          which(SP_df_raw_fix_mig[,"Spp142.1"]%in%c(99)),
                                          100)

# -> change change 101 to 102
SP_df_raw_fix_mig[,"Spp142"] <- replace(SP_df_raw_fix_mig[,"Spp142"],
                                        which(SP_df_raw_fix_mig[,"Spp142"]%in%c(101)),
                                        102)

SP_df_raw_fix_mig[,"Spp142.1"] <- replace(SP_df_raw_fix_mig[,"Spp142.1"],
                                          which(SP_df_raw_fix_mig[,"Spp142.1"]%in%c(101)),
                                          102)

# -> change change 105 to 106
SP_df_raw_fix_mig[,"Spp142"] <- replace(SP_df_raw_fix_mig[,"Spp142"],
                                        which(SP_df_raw_fix_mig[,"Spp142"]%in%c(105)),
                                        106)

SP_df_raw_fix_mig[,"Spp142.1"] <- replace(SP_df_raw_fix_mig[,"Spp142.1"],
                                          which(SP_df_raw_fix_mig[,"Spp142.1"]%in%c(105)),
                                          106)

################################################################################
###################### Spp231
################################################################################

######## General

######## Migraine

# -> change 124 to 125
SP_df_raw_fix_mig[,"Spp231"] <- replace(SP_df_raw_fix_mig[,"Spp231"],
                                        which(SP_df_raw_fix_mig[,"Spp231"]%in%c(124)),
                                        125)

SP_df_raw_fix_mig[,"Spp231.1"] <- replace(SP_df_raw_fix_mig[,"Spp231.1"],
                                          which(SP_df_raw_fix_mig[,"Spp231.1"]%in%c(124)),
                                          125)

################################################################################
###################### Spp273
################################################################################

######## General

######## Migraine

# -> change 228 to 229
SP_df_raw_fix_mig[,"Spp273"] <- replace(SP_df_raw_fix_mig[,"Spp273"],
                                        which(SP_df_raw_fix_mig[,"Spp273"]%in%c(228)),
                                        229)

SP_df_raw_fix_mig[,"Spp273.1"] <- replace(SP_df_raw_fix_mig[,"Spp273.1"],
                                          which(SP_df_raw_fix_mig[,"Spp273.1"]%in%c(228)),
                                          229)

# -> change 236 to 235
SP_df_raw_fix_mig[,"Spp273"] <- replace(SP_df_raw_fix_mig[,"Spp273"],
                                        which(SP_df_raw_fix_mig[,"Spp273"]%in%c(236)),
                                        235)

SP_df_raw_fix_mig[,"Spp273.1"] <- replace(SP_df_raw_fix_mig[,"Spp273.1"],
                                          which(SP_df_raw_fix_mig[,"Spp273.1"]%in%c(236)),
                                          235)

# -> change 242 to 241
SP_df_raw_fix_mig[,"Spp273"] <- replace(SP_df_raw_fix_mig[,"Spp273"],
                                        which(SP_df_raw_fix_mig[,"Spp273"]%in%c(242)),
                                        241)

SP_df_raw_fix_mig[,"Spp273.1"] <- replace(SP_df_raw_fix_mig[,"Spp273.1"],
                                          which(SP_df_raw_fix_mig[,"Spp273.1"]%in%c(242)),
                                          241)

# -> change 246 to NA
SP_df_raw_fix_mig[,"Spp273"] <- replace(SP_df_raw_fix_mig[,"Spp273"],
                                        which(SP_df_raw_fix_mig[,"Spp273"]%in%c(246)),
                                        NA)

SP_df_raw_fix_mig[,"Spp273.1"] <- replace(SP_df_raw_fix_mig[,"Spp273.1"],
                                          which(SP_df_raw_fix_mig[,"Spp273.1"]%in%c(246)),
                                          NA)

# -> change 254 to 253
SP_df_raw_fix_mig[,"Spp273"] <- replace(SP_df_raw_fix_mig[,"Spp273"],
                                        which(SP_df_raw_fix_mig[,"Spp273"]%in%c(254)),
                                        253)

SP_df_raw_fix_mig[,"Spp273.1"] <- replace(SP_df_raw_fix_mig[,"Spp273.1"],
                                          which(SP_df_raw_fix_mig[,"Spp273.1"]%in%c(254)),
                                          253)

# -> change 260 to 259
SP_df_raw_fix_mig[,"Spp273"] <- replace(SP_df_raw_fix_mig[,"Spp273"],
                                        which(SP_df_raw_fix_mig[,"Spp273"]%in%c(260)),
                                        259)

SP_df_raw_fix_mig[,"Spp273.1"] <- replace(SP_df_raw_fix_mig[,"Spp273.1"],
                                          which(SP_df_raw_fix_mig[,"Spp273.1"]%in%c(260)),
                                          259)

# -> change 266 to 265
SP_df_raw_fix_mig[,"Spp273"] <- replace(SP_df_raw_fix_mig[,"Spp273"],
                                        which(SP_df_raw_fix_mig[,"Spp273"]%in%c(266)),
                                        265)

SP_df_raw_fix_mig[,"Spp273.1"] <- replace(SP_df_raw_fix_mig[,"Spp273.1"],
                                          which(SP_df_raw_fix_mig[,"Spp273.1"]%in%c(266)),
                                          265)

# -> change 272 to 271
SP_df_raw_fix_mig[,"Spp273"] <- replace(SP_df_raw_fix_mig[,"Spp273"],
                                        which(SP_df_raw_fix_mig[,"Spp273"]%in%c(272)),
                                        271)

SP_df_raw_fix_mig[,"Spp273.1"] <- replace(SP_df_raw_fix_mig[,"Spp273.1"],
                                          which(SP_df_raw_fix_mig[,"Spp273.1"]%in%c(272)),
                                          271)

# -> change 278 to 277
SP_df_raw_fix_mig[,"Spp273"] <- replace(SP_df_raw_fix_mig[,"Spp273"],
                                        which(SP_df_raw_fix_mig[,"Spp273"]%in%c(278)),
                                        277)

SP_df_raw_fix_mig[,"Spp273.1"] <- replace(SP_df_raw_fix_mig[,"Spp273.1"],
                                          which(SP_df_raw_fix_mig[,"Spp273.1"]%in%c(278)),
                                          277)

# -> change 284 to 283
SP_df_raw_fix_mig[,"Spp273"] <- replace(SP_df_raw_fix_mig[,"Spp273"],
                                        which(SP_df_raw_fix_mig[,"Spp273"]%in%c(284)),
                                        283)

SP_df_raw_fix_mig[,"Spp273.1"] <- replace(SP_df_raw_fix_mig[,"Spp273.1"],
                                          which(SP_df_raw_fix_mig[,"Spp273.1"]%in%c(284)),
                                          283)

# -> change 287 to 286
SP_df_raw_fix_mig[,"Spp273"] <- replace(SP_df_raw_fix_mig[,"Spp273"],
                                        which(SP_df_raw_fix_mig[,"Spp273"]%in%c(287)),
                                        286)

SP_df_raw_fix_mig[,"Spp273.1"] <- replace(SP_df_raw_fix_mig[,"Spp273.1"],
                                          which(SP_df_raw_fix_mig[,"Spp273.1"]%in%c(287)),
                                          286)

# -> change 290 to 289
SP_df_raw_fix_mig[,"Spp273"] <- replace(SP_df_raw_fix_mig[,"Spp273"],
                                        which(SP_df_raw_fix_mig[,"Spp273"]%in%c(290)),
                                        289)

SP_df_raw_fix_mig[,"Spp273.1"] <- replace(SP_df_raw_fix_mig[,"Spp273.1"],
                                          which(SP_df_raw_fix_mig[,"Spp273.1"]%in%c(290)),
                                          289)
# -> change 279 to NA
SP_df_raw_fix_mig[,"Spp273"] <- replace(SP_df_raw_fix_mig[,"Spp273"],
                                        which(SP_df_raw_fix_mig[,"Spp273"]%in%c(279)),
                                        NA)

SP_df_raw_fix_mig[,"Spp273.1"] <- replace(SP_df_raw_fix_mig[,"Spp273.1"],
                                          which(SP_df_raw_fix_mig[,"Spp273.1"]%in%c(279)),
                                          NA)

################################################################################
###################### Spp476
################################################################################

######## General

######## Migraine

# -> change 104 to 105
SP_df_raw_fix_mig[,"Spp476"] <- replace(SP_df_raw_fix_mig[,"Spp476"],
                                        which(SP_df_raw_fix_mig[,"Spp476"]%in%c(104)),
                                        105)

SP_df_raw_fix_mig[,"Spp476.1"] <- replace(SP_df_raw_fix_mig[,"Spp476.1"],
                                          which(SP_df_raw_fix_mig[,"Spp476.1"]%in%c(104)),
                                          105)

# -> change 100 to 101
SP_df_raw_fix_mig[,"Spp476"] <- replace(SP_df_raw_fix_mig[,"Spp476"],
                                        which(SP_df_raw_fix_mig[,"Spp476"]%in%c(100)),
                                        101)

SP_df_raw_fix_mig[,"Spp476.1"] <- replace(SP_df_raw_fix_mig[,"Spp476.1"],
                                          which(SP_df_raw_fix_mig[,"Spp476.1"]%in%c(100)),
                                          101)

# -> change 102 to 103
SP_df_raw_fix_mig[,"Spp476"] <- replace(SP_df_raw_fix_mig[,"Spp476"],
                                        which(SP_df_raw_fix_mig[,"Spp476"]%in%c(102)),
                                        103)

SP_df_raw_fix_mig[,"Spp476.1"] <- replace(SP_df_raw_fix_mig[,"Spp476.1"],
                                          which(SP_df_raw_fix_mig[,"Spp476.1"]%in%c(102)),
                                          103)

# -> change 110 to NA
SP_df_raw_fix_mig[,"Spp476"] <- replace(SP_df_raw_fix_mig[,"Spp476"],
                                        which(SP_df_raw_fix_mig[,"Spp476"]%in%c(110)),
                                        NA)

SP_df_raw_fix_mig[,"Spp476.1"] <- replace(SP_df_raw_fix_mig[,"Spp476.1"],
                                          which(SP_df_raw_fix_mig[,"Spp476.1"]%in%c(110)),
                                          NA)

################################################################################
###################### Spp051
################################################################################

######## General

######## Migraine

# -> change 120 to 119
SP_df_raw_fix_mig[,"Spp051"] <- replace(SP_df_raw_fix_mig[,"Spp051"],
                                        which(SP_df_raw_fix_mig[,"Spp051"]%in%c(120)),
                                        119)

SP_df_raw_fix_mig[,"Spp051.1"] <- replace(SP_df_raw_fix_mig[,"Spp051.1"],
                                          which(SP_df_raw_fix_mig[,"Spp051.1"]%in%c(120)),
                                          119)

# -> change 124 to NA
SP_df_raw_fix_mig[,"Spp051"] <- replace(SP_df_raw_fix_mig[,"Spp051"],
                                        which(SP_df_raw_fix_mig[,"Spp051"]%in%c(124)),
                                        NA)

SP_df_raw_fix_mig[,"Spp051.1"] <- replace(SP_df_raw_fix_mig[,"Spp051.1"],
                                          which(SP_df_raw_fix_mig[,"Spp051.1"]%in%c(124)),
                                          NA)

################################################################################
###################### Spp108
################################################################################

######## General

# -> change 191 to NA
SP_df_raw_fix[,"Spp108"] <- replace(SP_df_raw_fix[,"Spp108"],
                                    which(SP_df_raw_fix[,"Spp108"]%in%c(191)),
                                    NA)

SP_df_raw_fix[,"Spp108.1"] <- replace(SP_df_raw_fix[,"Spp0108.1"],
                                      which(SP_df_raw_fix[,"Spp0108.1"]%in%c(191)),
                                      NA)

######## Migraine NTO FINISHED BECAUSE WE GET RID OF IT ANYWAY

# -> change 191 to NA
SP_df_raw_fix_mig[,"Spp108"] <- replace(SP_df_raw_fix_mig[,"Spp108"],
                                    which(SP_df_raw_fix_mig[,"Spp108"]%in%c(191)),
                                    NA)

SP_df_raw_fix_mig[,"Spp108.1"] <- replace(SP_df_raw_fix_mig[,"Spp108.1"],
                                      which(SP_df_raw_fix_mig[,"Spp108.1"]%in%c(191)),
                                      NA)

################################################################################
###################### Spp141
################################################################################

######## General

######## Migraine

# -> change 123 to 124
SP_df_raw_fix_mig[,"Spp141"] <- replace(SP_df_raw_fix_mig[,"Spp141"],
                                        which(SP_df_raw_fix_mig[,"Spp141"]%in%c(123)),
                                        124)

SP_df_raw_fix_mig[,"Spp141.1"] <- replace(SP_df_raw_fix_mig[,"Spp141.1"],
                                          which(SP_df_raw_fix_mig[,"Spp141.1"]%in%c(123)),
                                          124)

################################################################################
###################### Spp313
################################################################################

######## General

# -> change 233 to 232
SP_df_raw_fix[,"Spp313"] <- replace(SP_df_raw_fix[,"Spp313"],
                                    which(SP_df_raw_fix[,"Spp313"]%in%c(233)),
                                    232)

SP_df_raw_fix[,"Spp313.1"] <- replace(SP_df_raw_fix[,"Spp313.1"],
                                      which(SP_df_raw_fix[,"Spp313.1"]%in%c(233)),
                                      232)

######## Migraine

# -> change 215 to 216
SP_df_raw_fix_mig[,"Spp313"] <- replace(SP_df_raw_fix_mig[,"Spp313"],
                                        which(SP_df_raw_fix_mig[,"Spp313"]%in%c(215)),
                                        216)

SP_df_raw_fix_mig[,"Spp313.1"] <- replace(SP_df_raw_fix_mig[,"Spp313.1"],
                                          which(SP_df_raw_fix_mig[,"Spp313.1"]%in%c(215)),
                                          216)

# -> change 233 to 232
SP_df_raw_fix_mig[,"Spp313"] <- replace(SP_df_raw_fix_mig[,"Spp313"],
                                        which(SP_df_raw_fix_mig[,"Spp313"]%in%c(233)),
                                        232)

SP_df_raw_fix_mig[,"Spp313.1"] <- replace(SP_df_raw_fix_mig[,"Spp313.1"],
                                          which(SP_df_raw_fix_mig[,"Spp313.1"]%in%c(233)),
                                          232)

# -> change 237 to NA
SP_df_raw_fix_mig[,"Spp313"] <- replace(SP_df_raw_fix_mig[,"Spp313"],
                                        which(SP_df_raw_fix_mig[,"Spp313"]%in%c(237)),
                                        NA)

SP_df_raw_fix_mig[,"Spp313.1"] <- replace(SP_df_raw_fix_mig[,"Spp313.1"],
                                          which(SP_df_raw_fix_mig[,"Spp313.1"]%in%c(237)),
                                          NA)

# -> change 229 to NA
SP_df_raw_fix_mig[,"Spp313"] <- replace(SP_df_raw_fix_mig[,"Spp313"],
                                        which(SP_df_raw_fix_mig[,"Spp313"]%in%c(229)),
                                        NA)

SP_df_raw_fix_mig[,"Spp313.1"] <- replace(SP_df_raw_fix_mig[,"Spp313.1"],
                                          which(SP_df_raw_fix_mig[,"Spp313.1"]%in%c(229)),
                                          NA)

################################################################################
###################### Spp360
################################################################################

######## General

# -> change 129 to 130
SP_df_raw_fix[,"Spp360"] <- replace(SP_df_raw_fix[,"Spp360"],
                                    which(SP_df_raw_fix[,"Spp360"]%in%c(129)),
                                    130)

SP_df_raw_fix[,"Spp360.1"] <- replace(SP_df_raw_fix[,"Spp360.1"],
                                      which(SP_df_raw_fix[,"Spp360.1"]%in%c(129)),
                                      130)

# -> change 119 to 120
SP_df_raw_fix[,"Spp360"] <- replace(SP_df_raw_fix[,"Spp360"],
                                    which(SP_df_raw_fix[,"Spp360"]%in%c(119)),
                                    120)

SP_df_raw_fix[,"Spp360.1"] <- replace(SP_df_raw_fix[,"Spp360.1"],
                                      which(SP_df_raw_fix[,"Spp360.1"]%in%c(119)),
                                      120)

######## Migraine

# -> change 129 to 130
SP_df_raw_fix_mig[,"Spp360"] <- replace(SP_df_raw_fix_mig[,"Spp360"],
                                    which(SP_df_raw_fix_mig[,"Spp360"]%in%c(129)),
                                    130)

SP_df_raw_fix_mig[,"Spp360.1"] <- replace(SP_df_raw_fix_mig[,"Spp360.1"],
                                      which(SP_df_raw_fix_mig[,"Spp360.1"]%in%c(129)),
                                      130)

# -> change 119 to 120
SP_df_raw_fix_mig[,"Spp360"] <- replace(SP_df_raw_fix_mig[,"Spp360"],
                                    which(SP_df_raw_fix_mig[,"Spp360"]%in%c(119)),
                                    120)

SP_df_raw_fix_mig[,"Spp360.1"] <- replace(SP_df_raw_fix_mig[,"Spp360.1"],
                                      which(SP_df_raw_fix_mig[,"Spp360.1"]%in%c(119)),
                                      120)

# -> change 123 to NA
SP_df_raw_fix_mig[,"Spp360"] <- replace(SP_df_raw_fix_mig[,"Spp360"],
                                        which(SP_df_raw_fix_mig[,"Spp360"]%in%c(123)),
                                        NA)

SP_df_raw_fix_mig[,"Spp360.1"] <- replace(SP_df_raw_fix_mig[,"Spp360.1"],
                                          which(SP_df_raw_fix_mig[,"Spp360.1"]%in%c(123)),
                                          NA)

# -> change 124 to NA
SP_df_raw_fix_mig[,"Spp360"] <- replace(SP_df_raw_fix_mig[,"Spp360"],
                                        which(SP_df_raw_fix_mig[,"Spp360"]%in%c(124)),
                                        NA)

SP_df_raw_fix_mig[,"Spp360.1"] <- replace(SP_df_raw_fix_mig[,"Spp360.1"],
                                          which(SP_df_raw_fix_mig[,"Spp360.1"]%in%c(124)),
                                          NA)

################################################################################
###################### Spp391
################################################################################

######## General

######## Migraine

# -> change 152 to NA
SP_df_raw_fix_mig[,"Spp391"] <- replace(SP_df_raw_fix_mig[,"Spp391"],
                                        which(SP_df_raw_fix_mig[,"Spp391"]%in%c(152)),
                                        NA)

SP_df_raw_fix_mig[,"Spp391.1"] <- replace(SP_df_raw_fix_mig[,"Spp391.1"],
                                          which(SP_df_raw_fix_mig[,"Spp391.1"]%in%c(152)),
                                          NA)

# -> change 166 to NA
SP_df_raw_fix_mig[,"Spp391"] <- replace(SP_df_raw_fix_mig[,"Spp391"],
                                        which(SP_df_raw_fix_mig[,"Spp391"]%in%c(166)),
                                        NA)

SP_df_raw_fix_mig[,"Spp391.1"] <- replace(SP_df_raw_fix_mig[,"Spp391.1"],
                                          which(SP_df_raw_fix_mig[,"Spp391.1"]%in%c(166)),
                                          NA)
# Get proper names back
names(SP_df_raw_fix_mig) <- names(SP_df_raw_fix)

################################################################################
# GET RID OF LOCI WITH IRREPARABLE WEIRD STEP SIZES (base disruption etc)

SP_df_raw_fix_mig[,"Spp273"] <- NULL
SP_df_raw_fix_mig[,"Spp273.1"] <- NULL
SP_df_raw_fix_mig[,"Spp108"] <- NULL
SP_df_raw_fix_mig[,"Spp108.1"] <- NULL
