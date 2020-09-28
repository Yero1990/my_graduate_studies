#codee to aacount for the difference in thickness
#between Al. dummy an LH2 or LD2

Al_density = 2.7 #g / cm^3

#Aluminum wall of LH2 / LD2 cryotargets
LH2_ent = 0.0150 #cm
LH2_ent_err = 0.0011 #cm
LH2_ext = (0.0191+0.0219)/2.  #cm: (tip + wall)/2

LD2_ent = 0.0130 #cm
LD2_ent_err = 0.0012 #cm
LD2_ext = (0.0188+0.0184)/2.  #cm: (tip + wall)/2

LH2_thick = Al_density * (LH2_ent + LH2_ext)/2. # g / cm^2
LD2_thick = Al_density * (LD2_ent + LD2_ext)/2. # g / cm^2


#Aluminum Dummy (al. foils at z = +/- 5 cm) to mimic cryotarget walls
Al_ent = 0.1816 # g / cm^2
Al_ext = 0.1815 # g / cm^2
Al_thick = (Al_ent+Al_ext)/2.

#RATIO
LH2_factor =  LH2_thick / Al_thick
LD2_factor =  LD2_thick / Al_thick

#The ratio represents the thickness ratio of cryogenic as compared to Al. dummy target
#This factor should be multiplied by the Al. dummy yield to convert the thickness from
#Al DUMMY to CRYOTARGET END-CAPS THICKNESS

print('LH2_thick = ',LH2_thick,' g/cm^2')
print('LD2_thick = ',LD2_thick,' g/cm^2')
print('Al_thick = ',Al_thick,' g/cm^2')

print('LH2_factor = ', LH2_factor)
print('LD2_factor = ', LD2_factor)
