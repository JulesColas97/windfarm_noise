from prepost import Source 
import matplotlib.pyplot as plt 

src = Source()
src.load("./swl.dat")
# src.load("./swl.dat")
src.plot_all_u()
# plt.savefig('flow.png')
plt.show()


quit()

