
from matplotlib import pyplot as plt 


from rdkit import Chem
from rdkit.Chem import Descriptors






inhibitors = [('Methotrexate','CN(CC1=CN=C2C(=N1)C(=NC(=N2)N)N)C3=CC=C(C=C3)C(=O)N[C@@H](CCC(=O)O)C(=O)O')
,('Trimethoprim','COC1=CC(=CC(=C1OC)OC)CC2=CN=C(N=C2N)N'),
('Aminopterin', 'C1=CC(=CC=C1C(=O)N[C@@H](CCC(=O)O)C(=O)O)NCC2=CN=C3C(=N2)C(=NC(=N3)N)N'),

('Pralatrexate', 'CCCC(CC1=CN=C2C(=N1)C(=NC(=N2)N)N)C3=CC=C(C=C3)C(=O)N[C@@H](CCC(=O)O)C(=O)O'),

('Metoprine' ,'CC1=C(C(=NC(=N1)N)N)C2=CC(=C(C=C2)Cl)Cl'),
('WR99210', 'CC1(N=C(N=C(N1OCCCOC2=CC(=C(C=C2Cl)Cl)Cl)N)N)C')
        
         

]
mahi=[] #mahi is the name of the list which contains logp values of our drugs 
i=0         
for name,smiles in inhibitors: 
  mol = Chem.MolFromSmiles(smiles)
  if mol :
    logp=Descriptors.MolLogP(mol) 
    
    mahi.append(logp)
    
  i=i+1 

print(mahi)

sri=[] #sri is the name of the list which contains tpsa values
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
y=[-1.34,-1.732,-0.875,-0.857,-2,1.124]
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






