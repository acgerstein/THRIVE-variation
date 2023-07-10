library(diskImageR)

#Input directory: /Users/acgerstein/Nextcloud/Umanitoba/Research/18Rebekah/Research/THRIVE/Research/Pilot-VariationStudy-YST/DDA/Photos/pH4_48h/ImageJ_Adjusted_DDA1_AG/
IJMacro("pH4_DDA1_t")
maxLik("pH4_DDA1_t", 1)
createDataframe("pH4_DDA1_t", 1)

#Input directory: /Users/acgerstein/Nextcloud/Umanitoba/Research/18Rebekah/Research/THRIVE/Research/Pilot-VariationStudy-YST/DDA/Photos/pH4_48h/ImageJ_Adjusted_DDA2_AG/
IJMacro("pH4_DDA2_t")
maxLik("pH4_DDA2_t", 1)
createDataframe("pH4_DDA2_t", 1)

#Input directory: /Users/acgerstein/Nextcloud/Umanitoba/Research/18Rebekah/Research/THRIVE/Research/Pilot-VariationStudy-YST/DDA/Photos/pH4_48h/ImageJ_Adjusted_DDA3_AG/
IJMacro("pH4_DDA3_t")
maxLik("pH4_DDA3_t", 1)
createDataframe("pH4_DDA3_t", 1)

#Input directory: /Users/acgerstein/Nextcloud/Umanitoba/Research/18Rebekah/Research/THRIVE/Research/Pilot-VariationStudy-YST/DDA/Photos/pH4_48h/ImageJ_Adjusted_DDA4_AG/
IJMacro("pH4_DDA4_t")
maxLik("pH4_DDA4_t", 1)
createDataframe("pH4_DDA4_t", 1)

#Input directory: /Users/acgerstein/Nextcloud/Umanitoba/Research/18Rebekah/Research/THRIVE/Research/Pilot-VariationStudy-YST/DDA/Photos/pH4_48h/ImageJ_Adjusted_DDA5_AG/
IJMacro("pH4_DDA5_t")
maxLik("pH4_DDA5_t", 1)
createDataframe("pH4_DDA5_t", 1)

#Input directory: /Users/acgerstein/Nextcloud/Umanitoba/Research/18Rebekah/Research/THRIVE/Research/Pilot-VariationStudy-YST/DDA/Photos/pH7_48h/210515_DDA_YST6_V7-18_R12-24_AG/
IJMacro("pH7_DDA1_t")
maxLik("pH7_DDA1_t", 1)
createDataframe("pH7_DDA1", 1)

#Input directory: /Users/acgerstein/Nextcloud/Umanitoba/Research/18Rebekah/Research/THRIVE/Research/Pilot-VariationStudy-YST/DDA/Photos/pH7_48h/210520_DDA_YST6_V1-6,19-24_R1-12_AG/
IJMacro("pH7_DDA2")
maxLik("pH7_DDA2", 1, maxDist = 20)
createDataframe("pH7_DDA2", 1)
createDataframe("pH7_DDA2", 1)

#Input directory: /Users/acgerstein/Nextcloud/Umanitoba/Research/18Rebekah/Research/THRIVE/Research/Pilot-VariationStudy-YST/DDA/Photos/pH7_48h/210603_DDA_YST6_V7-24_R13-24_AG/
IJMacro("pH7_DDA3")
maxLik("pH7_DDA3", 1)
createDataframe("pH7_DDA3", 1)

#Input directory: /Users/acgerstein/Nextcloud/Umanitoba/Research/18Rebekah/Research/THRIVE/Research/Pilot-VariationStudy-YST/DDA/Photos/pH7_48h/210610_DDA_YST7_AG/
IJMacro("pH7_DDA4")
maxLik("pH7_DDA4", 1)
createDataframe("pH7_DDA4", 1)

#Input directory: /Users/acgerstein/Nextcloud/Umanitoba/Research/18Rebekah/Research/THRIVE/Research/Pilot-VariationStudy-YST/DDA/Photos/pH7_48h/210619_DDA_YST7_AG/
IJMacro("pH7_DDA5")
maxLik("pH7_DDA5", 1)
createDataframe("pH7_DDA5", 1)

#Input directory: /Users/acgerstein/Nextcloud/Umanitoba/Research/18Rebekah/Research/THRIVE/Research/Pilot-VariationStudy-YST/DDA/Photos/pH7_48h/210701_DDA_YST7_AG/
IJMacro("pH7_DDA6")
maxLik("pH7_DDA6", 1)
createDataframe("pH7_DDA6", 1)

IJMacro("combined")
maxLik("combined", 1)
createDataframe("combined", 1)

IJMacro("pH7_YST6_FLC")
readInExistingIJ("pH7_YST6_FLC")
maxLik("pH7_YST6_FLC", 1)

maxLik("combined", 1)
createDataframe("combined", 1)


