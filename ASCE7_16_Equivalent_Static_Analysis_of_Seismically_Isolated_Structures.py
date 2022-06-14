"""
* Open source (MIT license) python package for conducting equivalent static analysis on seismically isolated structures.
* The procedure is according to ASCE7-16 design requirements.
* Developed by Amin Moghaddasi SEI member of ASCE. Member I.D. #: 1221458
* I would be happy if you contribute further to this project.
* Reach me on Github: https://github.com/Amin-Moghaddasi 
* Also pardon me explaining basic things in comments and every docstring. Those explanations are meant to be for the beginer designers using this code.
"""
from math import pow as power # Function returning x to the power of y (power(x,y)). Example: power(3,4)= 3^4=81
from math import sqrt as square_root # Function returning the square root of a given number. Example: square_root(4)=2
from math import pi # Function returning value of the number Pi
class Initials: # Initial necessary calculations
    @staticmethod
    def Response_modification_coefficient(seismic_force_resisting_system_classification,type_number):
        """Returns the response modification coefficient for the chosen type of force-resisting system 
           present in the given structure (R)

        Args:
            seismic_force_resisting_system_classification (string): Designation pointed out in Table 12.2-1, ASCE7-16
                                                                    ('A','B','C','D','E','F','G','H')
            type_number (int): Row number under the chosen designation (for 'F' and 'H' this number is equal to 0)
        """
        a=[5,4,2,1.5,4,3,5,3.5,2,2,1.5,1.5,2,1.5,6.5,6.5,2,4]
        b=[8,6,3.25,6,5,2,1.5,5,4,8,5,3,6.5,6,5,5.5,4,2,2,1.5,1.5,7,7,2.5,8,7]
        c=[8,7,4.5,3.5,8,5,3,8,5,6,3,3.5]
        d=[8,7,7,6,8,6,7.5,7,6,5.5,4,8,8]
        e=[6,6.5,3,3.5,5.5,3.5,5,5.5]
        f=4.5
        g=[2.5,1.25,2.5,1.5,1,1.5]
        h=3
        if seismic_force_resisting_system_classification=='A' and 1<=type_number and type_number<=18:
            return a[type_number-1]
        elif seismic_force_resisting_system_classification=='B' and 1<=type_number and type_number<=26:
            return b[type_number-1]
        elif seismic_force_resisting_system_classification=='C' and 1<=type_number and type_number<=12:
            return c[type_number-1]
        elif seismic_force_resisting_system_classification=='D' and 1<=type_number and type_number<=13:
            return d[type_number-1]
        elif seismic_force_resisting_system_classification=='E' and 1<=type_number and type_number<=8:
            return e[type_number-1]
        elif seismic_force_resisting_system_classification=='F':
            return f
        elif seismic_force_resisting_system_classification=='G'and 1<=type_number and type_number<=6:
            return g[type_number-1]
        elif seismic_force_resisting_system_classification=='H':
            return h
        else:
            print("*Error! the designation, or row number, or both do not exist in the list! revise your input!")
            return 0
    @staticmethod
    def Any_B_Value(Any_beta):
        """Calculates desired B value according to entered beta value and Table 17.5-1

        Args:
            Any_beta (double): Any beta (damping) value
        """        
        if Any_beta>=0 and Any_beta<=0.02:
            Any_B = 0.8
        elif Any_beta>0.02 and Any_beta<=0.05:
            Any_B = (6.667*Any_beta)+0.667
        elif Any_beta>0.05 and Any_beta<=0.1:
            Any_B = (3.999*Any_beta)+0.8
        elif Any_beta>0.1 and Any_beta<=0.2:
            Any_B = (3*Any_beta)+0.899
        elif Any_beta>0.2 and Any_beta<=0.3:
            Any_B = (2*Any_beta)+1.1
        elif Any_beta>0.3 and Any_beta<=0.4:
            Any_B = (1.999*Any_beta)+1.1
        elif Any_beta>0.4 and Any_beta<=0.5:
            Any_B = (0.510*Any_beta)+1.745
        elif Any_beta>0.5:
            Any_B = 2
        else:
            print('**Error! Entered beta value is out of range! Revise the input!')
            return 0
        return round(Any_B,3)    
class Isolation_System: # Calculations related to the isolation system
    @staticmethod
    def Initial_Stiffness(Dy,Fy):
        """Calculates the initial stiffness value for the isolator device assuming its stress-strain curve as a bilinear model (Ki)

        Args:
            Dy (double): Yield displacement/strain
            Fy (double): Yield force/stress
        """        
        K1 = Fy/Dy
        return round(K1,3)
    @staticmethod
    def Secondary_Stiffness(Dy,Fy,DM,FM):
        """Calculates the secondary stiffness value for the isolator device assuming its stress-strain curve as a bilinear model (Ks)

        Args:
            Dy (double): Yield displacement/strain
            Fy (double): Yield force/stress
            DM (double): Maximum (ultimate) displacement/strain
            FM (double): Maximum (ultimate) force/stress
        """        
        Ks = (FM-Fy)/(DM-Dy)
        return round(Ks,3)
    @staticmethod
    def Effective_Stiffness_Individual_Device(F_positive,F_negative,Delta_positive,Delta_negative):
        """Calculates the effective stiffness of a given isolator device (K_eff)

        Args:
            F_positive (double): Maximum (positive) force/stress
            F_negative (double): Minimum (negative) force/stress
            Delta_positive (double): Maximum (positive) displacement/strain 
            Delta_negative (double): Minimum (negative) displacement/strain
        """                
        K_eff = (abs(F_positive)+abs(F_negative))/(abs(Delta_positive)+abs(Delta_negative))
        return round(K_eff,3)
    @staticmethod
    def Effective_Stiffness_Whole_Isolation_System(F_positive_list,F_negative_list,DM):
        """Calculates the effective stiffness of the whole isolation system at maximum displacements (kM)

        Args:
            F_positive_list (list): List containing positive force values for each isolator device at maximum displacement (DM)
            F_negative_list (list): List containing negative force values for each isolator device at maximum displacement (DM)
            DM (doubel): Maximum (ultimate) displacement/strain
        """        
        numerator = 0
        denominator = 2 * DM
        for i in range (len(F_positive_list)):
            numerator += abs(F_negative_list[i])+abs(F_negative_list[i])
        kM = numerator/denominator
        return round(kM,3)
    @staticmethod
    def Hysteresis_Loop_Dissipated_Energy(Fy,DM):
        """Calculates the total energy dissipated in a single cycle at displacement Dm by an isolator device (E_loop)

        Args:
            Fy (double): Yield force/stress
            DM (double): Maximum displacement/strain
        """        
        E_loop = 4*Fy*DM
        return round(E_loop,3)
    @staticmethod
    def Effective_Damping_of_the_Isolation_System(E_M_list,kM,DM):
        """Calculates the effective damping of the isolation system at maximum displacement (beta_M)_summary_

        Args:
            E_M_list (list): List containing energy dissipation values for each isolator device at maximum displacement
            kM (double): The effective stiffness value of the whole isolation system
            DM (double): Maximum (ultimate) displacement/strain
        """        
        numerator = 0
        denominator = 2*pi*kM*power(DM,2)
        for i in range(len(E_M_list)):
            numerator += E_M_list[i]
        beta_M = numerator/denominator
        return round(beta_M,3)
    @staticmethod
    def Maximum_displacement(SM1,TM,BM,g):
        """Calculates the maximum displacment (DM)

        Args:
            SM1 (double):Maximum Common Earthquake (MCE) 5% damped spectral acceleration value at 1 second period
            TM (double): Effective period of the seismically isolated structure at the maximum displacement (DM)
            BM (double): Numerical coeeficient as set forth in Table 17.5-1 for the effective damping (beta_M)
            g (double): Gravitational acceleration value (units must comply)
        """        
        DM = (g*SM1*TM)/(4*power(pi,2)*BM)
        return round(DM,3)
    @staticmethod
    def PT_Ratio(Plan_width,Plan_length,N,X_list,Y_list):
        """Calculates the PT ratio (effective translational period of th isolation system / effective torsional period of the isolation system)

        Args:
            Plan_width (double): Width of the structure plan
            Plan_length (double): Length of the structure plan
            N (double): Number of isolator devices/units
            X_list (list):List containing values of distances between the center of mass of the structure and each isolator device/unit in X direction
            Y_list (list):List containing values of distances between the center of mass of the structure and each isolator device/unit in Y direction
        """        
        r_I = square_root((power(Plan_length,2)+power(Plan_width,2))/12)
        numerator = 0
        for i in range (N):
            numerator += power(X_list[i],2)+power(Y_list[i],2)
        PT = (square_root(numerator/N))/r_I
        return round(PT,3)
class Structure: #Calculations related to the structure itself
    @staticmethod
    def Effective_Period_at_Maximum_Displacement(W_sesimic,kM,g):
        """Calculates the effective period of the isolated structure at maximum displacement (TM)
        
        Args:
            W_sesimic (_type_): Effective seismic eight of the structure above the isolation interface
            kM (double): The effective stiffness value of the whole isolation system
            g (double): Gravitational acceleration value (units must comply)
        """        
        TM = (2*pi)*(square_root(W_sesimic/(kM*g)))
        return round(TM,3)
