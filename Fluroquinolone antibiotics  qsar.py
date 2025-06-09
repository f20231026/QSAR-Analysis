from matplotlib import pyplot as plt 


from rdkit import Chem
from rdkit.Chem import Descriptors



# drugs in market against DNA gyrase 
inhibitors = [('norfloxacin', 'CCN1C=C(C(=O)C2=CC(=C(C=C21)N3CCNCC3)F)C(=O)O'),
('ciprofloxacin', 'C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCNCC4)F)C(=O)O'),
('gatifloxacin', 'CC1CN(CCN1)C2=C(C=C3C(=C2OC)N(C=C(C3=O)C(=O)O)C4CC4)F'),
('moxifloxacin', 'COC1=C2C(=CC(=C1N3C[C@@H]4CCCN[C@@H]4C3)F)C(=O)C(=CN2C5CC5)C(=O)O')
        
         

]
mahi=[] # mahi is the name of the list which contains logp values of our drugs
i=0         
for name,smiles in inhibitors: 
  mol = Chem.MolFromSmiles(smiles)
  if mol :
    logp=Descriptors.MolLogP(mol) 
    
    mahi.append(logp)
    
  i=i+1 

print(mahi)

sri=[] # sri is the name of the list which contains tpsa of our drugs 
i=0         
for name,smiles in inhibitors: 
  mol = Chem.MolFromSmiles(smiles)
  if mol :
    tpsa=Descriptors.TPSA(mol) 
    
    sri.append(tpsa)
  i=i+1    
print(sri)
    
    

x1=mahi
x2=sri
y=[2.28,1.1,0.109,1]
plt.scatter(x1,y,label='logp vs logic50')
plt.scatter(x2,y,label='tpsa vs loic50 ')
plt.title('logp vs ic50')
plt.ylabel('ic50')
plt.xlabel('logp')
plt.legend()
plt.show()

from scipy.stats import linregress
slope, intercept, r_value, p_value, std_err = linregress(x1,y)
print(f"R²: {r_value**2}")


from scipy.stats import linregress
slope, intercept, r_value, p_value, std_err = linregress(x2,y)
print(f"R²: {r_value**2}")






