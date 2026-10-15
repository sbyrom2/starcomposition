 ## Tracking stellar composition

 Converting between [Fe/H], Z, X, Y, [C/Fe], [alpha/Fe], $A_i$ = log10($\frac{N_i}{N_H}$) + 12, etc is the bane of my existence

 So I made a class that will keep track of everything for me. Upon creation, an instance of the StarComp class has the Grevesse, Asplund, and Sauval 2007 solar mixture. One can procceed to change many of the values of the mixture, then get other values associated with the mixture (X,Y,Z, [Fe/H], etc).
 
 I have yet to add anything to change which solar mixture is used or to allow the user to change isotope ratios, but the scaffolding for it is there.


 Notes from creation:

 ## Abundances Object

What do I want this to be able to do?

Keep track of the composition of a "star" by tracking the number of each element


Class Attributes:
- $N_{tot}$ = number of particles (should be of order $10^12$, I think. To avoid losing precision for elements like Hf)
- $N_i$ = number of particles of element i
- $N_{i,\odot}$ = number of particles of element i for a solar mixture
- $R_{ij1j2}$ = for element i, mass ratio of isotope Mj1/Mj2 (used to calculate)
- $N_{ij,\odot}$ = number of particles of element i, isotope j for a solar mixture
- $m_i$ = atomic mass of element i
- $m_{ij}$ = atomic mass of element i, isotope j (such that $\Sigma_j m_{ij} \frac{N_{ij}}{N_i}  = m_i$)
- $m_{ij,\odot}$ = atomic mass of element i, isotope j for a solar mixture (such that $\Sigma_j m_{ij} \frac{N_{ij}}{N_i}  = m_i$)
- elem_dict = mapping of atomic Z to element abbreviation

Class SubAttributes (not a real term) (calculated when called for):
- $X_i$ = mass fraction of element i ($X_i = M_i / \Sigma_i M_i $)
- [Fe/H] = log($\frac{N_{Fe}/N_H}{N_{Fe,\odot}/N_{H,\odot}}$)

- XYZ = tuple of X,Y,Z: ($N_H m_H, N_{He} m_{He}, \Sigma_{i=Li} N_i m_i$) $/ \Sigma_I N_i m_i$
- $\mu_Z$ = average mass of metals ($\mu_Z = M_Z/N_Z$)
- $R_{i,j1,j2}$ = ratio of isotopes j1/j2 for element i: $N_{ij1}/N_{ij2}$
- $A$(elem $i$) = log($N_i/N_{H}$) + 12
- FRAC_C (N,O) = mass fraction of metals in C (N,O) ($M_C / M_{Z}$)

Class Methods:
- rescale_to_Ntot(): $N_{tot,new} = \Sigma_i N_i$, then for each elem i $N_i = N_i \frac{N_{tot,new}}{N_{tot}}$
	- also rescales isotopes: $N_{ij} = N_{ij}  \frac{N_{tot,new}}{N_{tot}}$
- change_feh(A): $N_{Fe} = N_{H} \frac{N_{Fe,\odot}}{N_{H,\odot}}10^{A}$, rescale_to_Ntot()
- change_X_fixedZ(Xnew): Ynew = 1-Z-Xnew, then do $N_H = X_{new} M_{tot} / m_H$, and $N_{He} = Y_{new} M_{tot} / m_{He}$, rescale_to_Ntot()
- change_Y_fixedZ(Ynew): Xnew = 1-Z-Ynew, then do $N_H = X_{new} M_{tot} / m_H$, and $N_{He} = Y_{new} M_{tot} / m_{He}$, rescale_to_Ntot()
- change_X_fixedY(Xnew): Znew = 1-Y-Xnew, then do $N_H = X_{new} M_{tot} / m_H$, $N_{Z,new} = Z_{new} M_{tot} / \mu_Z$, for k in metals $N_k = N_k \frac{N_{Z,new}}{N_{Z}}$, rescale_to_Ntot()
- change_Y_fixedX(Ynew): Znew = 1-X-Ynew, then do $N_{He} = Y_{new} M_{tot} / m_{He}$, $N_{Z,new} = Z_{new} M_{tot} / \mu_Z$, for k in metals $N_k = N_k \frac{N_{Z,new}}{N_{Z}}$, rescale_to_Ntot()
- set_XYZ(X,Y,Z): check X+Y+Z=1, then do change_X_fixedZ(X), change_Z_fixedX(Z)

- change_isotope_ratio(elem,iso1,iso2,r12): $i$ = elem_dict[elem]; let $N_{i,both,old}$ = $N_{ij1}+N_{ij2}$. $N_{ij1} = N_{ij1}/r_{12}$ and $N_{ij2} = N_{ij2} r_{12}$. Then do $N_{i,both,new} = N_{ij1} + N_{ij2}$, then $N_{ij1} = N_{ij1} \frac{N_{i,both,old}}{N_{i,both,new}}$ and then $N_{ij2} = N_{ij2} \frac{N_{i,both,old}}{N_{i,both,new}}$, rescale_to_Ntot()


TBD additions:
- add_mass(elem, amount): where amount is in units of $M_{tot}$, assuming the same isotope ratios as currently present
- switch_solarmix(newsolarmix): change to a new solar mixture. Note: might be safer to set the solar mixture when the object is instantiated and not allow changes later.
