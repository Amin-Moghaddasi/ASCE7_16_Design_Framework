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
    def R_I(R):
        """Calculates the isolation system response modification factor (R_I)

        Args:
            R (double): Response modification factor
        """        
        R_I = (R*3)/(8)
        if R_I<=1:
            R_I=1
        elif R_I>=2:
            R_I=2
        return round(R_I,3)
    @staticmethod
    def Story_Weight_Values(combined_load_list,story_area_list):
        """Calculates weight values for each story (wi)

        Args:
            combined_load (list): List containing combined vertical load values (per unit area) from first floor to top roof respectively
            story_area_list (list): List containing story area values from the first floor to the top roof, respectively 
        """       
        story_weight_list=[0]*len(story_area_list)
        if len(combined_load_list)==1: # If loads on each story have equal values
            for i in range(len(story_area_list)):
                story_weight_list[i]=combined_load_list[0]*story_area_list[i]
        else: # If each story has a separate load value
            for i in range (len(story_area_list)):
                story_weight_list[i]=combined_load_list[i]*story_area_list[i]
        return story_weight_list
    @staticmethod
    def Total_Weight_of_the_Structure(story_weight_values):
        """Calculates the total mass of the structure given the values of each story (W_bar)

        Args:
            story_weight_values (list): List containing story weight values from the first floor to the top roof,respectively
        """                
        total_mass=0
        for i in range(len(story_weight_values)):
            total_mass+=story_weight_values[i]
        return round(total_mass,3)
    @staticmethod
    def Seismic_Weight(DL,LL,story_area_list):
        """Calculates the effective seismic weigth of the structure (W_seismic)

        Args:
            DL (double): Dead load per unit area/length
            LL (double): Live load per unit area/length
            story_area_list (list): List containing story area values from the first floor to the top roof, respectively 

        Returns:
            _type_: _description_
        """        
        W_sesimic = 0
        w = DL + (0.2*LL)
        for i in range(len(story_area_list)):
            W_sesimic+= w*story_area_list[i]
        return round(W_sesimic,3)
    @staticmethod
    def Average_Vertical_Load(DL,LL):
        """Calcultes the avaerage combined load per unit area/length

        Args:
            DL (double): Dead load per unit area/length
            LL (double): Live load per unit area/length
        """        
        w = DL + (0.5*LL)
        return round(w,3)
    @staticmethod
    def Maximum_Vertical_Load(DL,LL,S,SDS,Eh):
        """Calcultes the avaerage combined load per unit area/length

        Args:
            DL (double): Dead load per unit area/length
            LL (double): Live load per unit area/length
            S (double): Snow load per unit area/length
            SDS (double): DBE 5% damped spectral acceleration value at 1 second period
            Eh (double): Horizontal earthquake load per unit area/length
        """
        Ev = 0.2*SDS*DL         
        w = (1.2*DL)+Ev+Eh+LL+(0.2*S)
        return round(w,3)
    @staticmethod
    def Minimum_Vertical_Load(DL,SDS,Eh):
        """Calcultes the avaerage combined load per unit area/length

        Args:
            DL (double): Dead load per unit area/length
            SDS (double): DBE 5% damped spectral acceleration value at 1 second period
            Eh (double): Horizontal earthquake load per unit area/length
        """
        Ev = 0.2*SDS*DL         
        w = (0.9*DL)-Ev+Eh
        return round(w,3)
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
    def Isolation_System_Lateral_Distributed_Force(Vb,Vst,R_I):
        """Calculates the lateral seismic force induced at base level (F_I)

        Args:
            Vb (double): Base shear
            Vst (double): Total unreduced lateral seismic design shear
            R_I (double): Isolation system response modification factor
        """        
        F_I = abs((Vb-Vst)/(R_I))
        return round(F_I,3)
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
    def Secondary_Stiffness(F1,F2,d1,d2):
        """Calculates the secondary stiffness value for the isolator device assuming its stress-strain curve as a bilinear model (K2)

        Args:
            F1 (double): Yield force
            F2 (double): Ultimate force
            d1 (double): Yield displacement
            d2 (double): Ultimate displacement

        Returns:
            _type_: _description_
        """        
        k2 = (F2-F1)/(d2-d1)
        return round(k2,3)
    @staticmethod
    def Effective_Stiffness_Individual_Device(F1,K2,d1,DM):
        """Calculates the effective stiffness of a given isolator device (K_M_eff)

        Args:
            F1 (double): Yield force
            K2 (double): Secondary stiffness
            d1 (double): Yield displacement
            DM (double): Maximum displacement

        Returns:
            _type_: _description_
        """                       
        K_eff = (F1+((DM-d1)*K2))/DM
        return round(K_eff,3)
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
    def Effective_Damping_of_the_Isolation_System(E_M,kM,DM):
        """Calculates the effective damping of the isolation system at maximum displacement (beta_M)

        Args:
            E_M (double): Energy dissipation value for each isolator device at maximum displacement
            kM (double): The effective stiffness value of the whole isolation system
            DM (double): Maximum (ultimate) displacement/strain
        """        
        denominator = 2*pi*kM*(power(DM,2))
        beta_M = E_M/denominator
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
    @staticmethod
    def Total_Maximum_Displacement(DM,y,PT,Plan_length,Plan_width,e):
        """Calculates the total maximum displacement of elements present in the isolation system (D_TM)

        Args:
            DM (double): Maximum displacement
            y (double): The distance between the centers of rigidity of the isolation system and the element of interest measured perpendicular to the direction of interest
            PT (double): Ratio of the effective translational period of the isolation system to the effective torsional period of the isolation system
            Plan_length (double): Length of the structure plan
            Plan_width (double): Width of the structure plan
            e (double): Actual eccentrity
        """        
        numerator = y*12*e
        denominator = power(PT,2)+(power(Plan_length,2)+power(Plan_width,2))
        fraction = numerator/denominator
        D_TM = DM * (1+fraction)
        return round(D_TM,3)
    @staticmethod
    def Isolation_System_Yield_Shear(Fy,R_I,N,Ws,W,beta_M,hysteresis_condition):
        """Calculates the isolation system shearing force (Vs_iso)

        Args:
            Fy (double): Isolator device yield force
            R_I (double): Isolation system response modification factor
            N (double): Number of Isolator devices
            Ws (double): Weigth of superstructure
            W (double): Total weight of the structure
            beta_M (double): Damping of isolation system
            hysteresis_condition (string): Enter "abrupt_transition" or "smooth_transition" depending on the isolator hysteresis type
        """        
        if hysteresis_condition == 'abrupt_transition':
            constant = 1-(3.5*beta_M)
        elif hysteresis_condition == 'smooth_transition':
            constant = 1-(2.5*beta_M)
        else:
            print('**Error! Check the input on hysteresis condition!')
            return 0
        constant1=Ws/W
        Vs_iso = (1.5*Fy*N*(power(constant1,constant)))/R_I
        return round(Vs_iso,3)
class Structure: #Calculations related to the structure itself
    @staticmethod
    def Effective_Period_at_Maximum_Displacement(W_sesimic,kM,g,N):
        """Calculates the effective period of the isolated structure at maximum displacement (TM)
        
        Args:
            W_sesimic (double): Effective seismic eight of the structure above the isolation interface
            kM (double): The effective stiffness value of the whole isolation system
            g (double): Gravitational acceleration value (units must comply)
            N (double): Number of columns with isolator devices
        """        
        TM = (2*pi)*(square_root(W_sesimic/(kM*g*N)))
        return round(TM,3)
    @staticmethod
    def Eccentrity(Plane_length,Plan_Width):
        """Calculates the eccentrity

        Args:
            Plane_length (double): Length of the structure plan
            Plan_Width (double): Width of the structure plan

        Returns:
            _type_: _description_
        """        
        if Plan_Width>Plane_length:
            a = Plan_Width
        else:
            a = Plane_length
        e = 0.05*a
        return round(e,3)
    @staticmethod
    def Minimum_Substructure_Shear(kM,DM):
        """Calculates the minimum shear imposed on substructure (Vb) 

        Args:
            kM (double): Effective stiffness of a given isolator device
            DM (double): Maximum displacement
        """        
        Vb = kM*DM
        return round(Vb,3)
    @staticmethod
    def SuperStructure_Weight(Total_Weight,Ground_Floor_Weight):
        """Calculates the weight of the superstructure

        Args:
            Total_Weight (double): Total weight of the structure
            Ground_Floor_Weight (double): Weight of the ground floor (directly above the base level)
        """        
        Ws = Total_Weight - Ground_Floor_Weight
        return round(Ws,3)
    @staticmethod
    def Total_Lateral_Seismic_Force(Vb,Ws,W,beta_M,hysteresis_condition):
        """Calculates the total unreduced lateral seismic force (Vst)

        Args:
            Vb (double): Substructure shear
            Ws (double): Weigth of superstructure
            W (double): Total weight of the structure
            beta_M (double): Damping of isolation system
            hysteresis_condition (string): Enter "abrupt_transition" or "smooth_transition" depending on the isolator hysteresis type
        """        
        if hysteresis_condition == 'abrupt_transition':
            constant = 1-(3.5*beta_M)
        elif hysteresis_condition == 'smooth_transition':
            constant = 1-(2.5*beta_M)
        else:
            print('**Error! Check the input on hysteresis condition!')
            return 0
        constant1=Ws/W
        Vst = Vb*(power(constant1,constant))
        if Vst<=Vb:
            Vst=Vb
        return round(Vst,3)
    @staticmethod
    def Fundamental_Approxiamte_Period(structure_type,total_structure_height):
        """"Calculates the approximate value of the period of fundamental mode for the given structure
            based on Table 12.8-2 of ASCE7-16 (Tfb)

        Args:
            structure_type (string): type 'MRFS', or 'MRFRC' for steel, and concrete moment resisting frames,respectively
                                    type 'SBRB' for steel eccentrically braced frames (utilizing BRBs)
                                    type 'other' if structure type is none of the previously described 
            total_structure_height (double): total height of the given structure in meters (base to top roof)
        """        
        if structure_type == 'MRFS':
            ct=0.0724
            x=0.8
        elif structure_type=='MRFC':
            ct=0.0466
            x=0.9
        elif structure_type=='SBRB':
            ct=0.0731
            x=0.75
        elif structure_type=='other':
            ct=0.0488
            x=0.75
        else:
            print("*Error! selected type of structure is not included in the list!")
            return 0
        T=ct*(power(total_structure_height,x))
        return round(T,3)
    @staticmethod
    def Fixed_Base_Shear(R,SDS,SD1,Ie,T1,W):
        """Calculates the fixed base shearing force (V_Fix)

        Args:
            R (double): Response modification factor
            SDS (double): Design, 5% damped, spectral response acceleration parameter at short periods
            SD1 (double): Design, 5% damped, spectral response acceleration parameter at a period of 1 s
            Ie (double): Importance factor
            T1 (double): Fundamental period of the structure
            W (double): Total seismic mass of the structure
        """        
        a = (SDS*W*Ie)/R
        b = (SD1*W*Ie)/(T1*R)
        if a<=b:
            V_Fix = a
        else:
            V_Fix = b
        return round(V_Fix,3)
    @staticmethod
    def Vertical_Force_Distribution_Factor(k,sup_story_weight_list,sup_story_height_list):
        """Calculates the vertical force distribution factor list (Cv)

        Args:
            k (double): Is equal to 14*beta_m*T_fixed_base
            sup_story_weight_list (list): List containing superstructure weight values (for each story)
            sup_story_height_list (_type_): List containing superstructure height values (for each story)
        """        
        Cv = [0]*len(sup_story_weight_list)
        denominator =0
        for i in range(len(sup_story_height_list)):
            denominator+= sup_story_weight_list[i]*power(sup_story_height_list[i],k)
        for i in range(len(sup_story_height_list)):
            Cv[i]=round(((sup_story_weight_list[i]*power(sup_story_height_list[i],k))/denominator),3)
        return Cv
    @staticmethod
    def Lateral_Seismic_Distributed_Force_List(Cv,Vs_min):
        """Calculates the list containing distributed lateral force values (F)

        Args:
            Cv (list): Vertical force distribution factor list
            Vs_min (double): Minimum design shearing force 
        """        
        F = [0]*len(Cv)
        for i in range(len(Cv)):
            F[i]= round((Cv[i]*Vs_min),3)
        return F
