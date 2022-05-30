"""
* Open source (MIT license) python package for conducting equivalent static analysis on structures with dampers as vibration control devices.
* The procedure is according to ASCE7-16 design requirements.
* Developed by Amin Moghaddasi SEI member of ASCE. Member I.D. #: 1221458
* I would be happy if you contribute further to this project.
* Reach me on Github: https://github.com/Amin-Moghaddasi 
* Also pardon me explaining basic things in comments and every docstring. Those explanations are meant to be for the beginer designers using this code.
"""
from math import pow as power # Function returning x to the power of y (power(x,y)). Example: power(3,4)= 3^4=81
from math import sqrt as square_root # Function returning the square root of a given number. Example: square_root(4)=2
from math import pi # Function returning value of the number Pi
class Initials: # A class for calculating the initial analysis values
    @staticmethod
    def Any_B_value(any_beta):
        """Calculates any B coefficient values (Table 18-7.2, ASCE7-16) given its corresponding beta value
             (B_v+1,B1_D,B1_M,Bm_D,Bm_M)

        Args:
            any_beta (double): Corresponding beta value for the desired B value
        """        
        if any_beta<=0.02:
            any_B=0.8
        elif any_beta>0.02 and any_beta<=0.05:
            any_B=(6.67*any_beta)+.667
        elif 0.05<any_beta and any_beta<=0.1:
            any_B=(3.99*any_beta)+0.8
        elif 0.1<any_beta and any_beta<=0.2:
            any_B=(3*any_beta)+0.899
        elif 0.2<any_beta and any_beta<=0.3:
            any_B=(3*any_beta)+0.899
        elif 0.3<any_beta and any_beta<=0.4:
            any_B=(2.99*any_beta)+0.9
        elif 0.4<any_beta and any_beta<=0.5:
            any_B=(2.99*any_beta)+0.9
        elif 0.5<any_beta and any_beta<=0.6:
            any_B=(3*any_beta)+0.899
        elif 0.6<any_beta and any_beta<=0.7:
            any_B=(2.99*any_beta)+0.9
        elif 0.7<any_beta and any_beta<=0.8:
            any_B=(2.99*any_beta)+0.9
        elif 0.8<any_beta and any_beta<=0.9:
            any_B=(3*any_beta)+0.899
        elif 0.9<any_beta and any_beta<1:
            any_B=(4*any_beta)
        elif 1<=any_beta:
            any_B=4  
        else:
            print('**Error! The entered beta value is not in the allowable range!')
            return 0
        return round(any_B,3)
    @staticmethod
    def Ts(sd1,sds):
        """Calculates the Ts ratio (Ts = SD1/SDS)

        Args:
            sd1 (double): Design, 5% damped, spectral response acceleration parameter at a period of 1 s
            sds (double): Design, 5% damped, spectral response acceleration parameter at a short periods
        """        
        Ts=sd1/sds
        return round(Ts,3)
    @staticmethod
    def Hysteresis_loop_adjustment_factor(Ts,T1):
        """Calculates the hysteresis loop adjusment factor (qH)

        Args:
            Ts (double): SD1 to SDs ratio
            T1 (double): Preiod of the fundamental mode of the structure
        """        
        qh=0.67*(Ts/T1)
        if qh>1:
            qh=1
            return qh
        elif qh<0.5:
            qh=0.5
            return qh
        else:
            return round(qh,3)
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
    def Over_strength_factor(seismic_force_resisting_system_classification,type_number):
        """Returns the over strength factor for the chosen type of force-resisting system 
           present in the given structure (Omega_0)

        Args:
            seismic_force_resisting_system_classification (string): Designation pointed out in Table 12.2-1, ASCE7-16
                                                                    ('A','B','C','D','E','F','G','H')
            type_number (int): Row number under the chosen designation (for 'F' and 'H' this number is equal to 0)
        """
        a=[2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,3,3,2.5,2]
        b=[2,2,2,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2,2,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2]
        c=[3,3,3,3,3,3,3,3,3,3,3,3]
        d=[2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,3,3,2.5,2.5]
        e=[2.5,2.5,3,3,2.5,2.5,3,2.5]
        f=2.5
        g=[1.25,1.25,1.25,1.25,1.25,1.5]
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
    def Deflection_amplification_factor(seismic_force_resisting_system_classification,type_number):
        """Returns the deflection amplification factor for the chosen type of force-resisting system 
           present in the given structure (Cd)

        Args:
            seismic_force_resisting_system_classification (string): Designation pointed out in Table 12.2-1, ASCE7-16
                                                                    ('A','B','C','D','E','F','G','H')
            type_number (int): Row number under the chosen designation (for 'F' and 'H' this number is equal to 0)
        """
        a=[5,4,2,1.5,4,3,3.5,2.25,1.75,1.75,1.25,1.75,2,1.5,4,4,2,3.5]
        b=[4,5,3.25,5,4.5,2,1.5,4.5,4,4,4.5,3,5.5,5,4.5,4,4,2,2,1.25,1.75,4.5,4.5,2.5,5,6,]
        c=[5.5,5.5,4,3,5.5,4.5,2.5,5.5,4.5,5.5,2.5,3.5]
        d=[4,5.5,5.5,5,4,5,6,6,5,5,3.5,5,6.5,]
        e=[5,5,2.5,3,4.5,3,4.5,4.5]
        f=4
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
    def Total_vertical_load_per_unit_area(dead_load,live_load,dead_load_combo_coefficient,live_load_combo_coefficient):
        """Calculates the total combined vertical load per unit area

        Args:
            dead_load (double): Imposed dead load (per unit area)
            live_load (double): Imposed live load (per unit area)
            dead_load_combo_coefficient (double): Dead load coefficient in the desired load combination
            live_load_combo_coefficient (double): Live load coefficient in the desired load combination
        """        
        total_combined_load=(dead_load*dead_load_combo_coefficient)+(live_load*live_load_combo_coefficient)
        return round(total_combined_load,3)
    @staticmethod
    def Story_weight_values(combined_load_list,story_area_list):
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
    def Total_mass_of_the_structure(story_weight_values):
        """Calculates the total mass of the structure given the values of each story (W_bar)

        Args:
            story_weight_values (list): List containing story weight values from the first floor to the top roof,respectively
        """                
        total_mass=0
        for i in range(len(story_weight_values)):
            total_mass+=story_weight_values[i]
        return round(total_mass,3)
    @staticmethod
    def Element_hysteretic_damping(beta_I,qh,mu_D):
        """Calculates the element hysteretic damping value (Beta_HD)

        Args:
            beta_I (double): Inherent damping (<=3%)
            qh (double): Hysteresis loop adjustment factor
            mu_D (double): Effective ductility demand on the seismic force-resisting system
        """    
        if beta_I>0.03 or beta_I<0:
            print("**Error! Inherent damping value must be in range 0<=beta_I<=0.03! revise your input!")
            return 0
        elif not 1<=mu_D:
            print("**Error! Ductility demand value must be equal to or greater than 1! revise your input!")
            return 0
        else:        
            beta_HD=qh*(0.64-beta_I)*(1-(1/mu_D))
            return round(beta_HD,3)
class Fundamental_Mode_Parameters: # A class for the fundamental-mode-related calculations
    @staticmethod
    def Accurate_fundamental_period(story_weight_list,story_force_list,delta_1D):
        numerator=0
        denominator=0
        for i in range(len(delta_1D)):
            numerator+=story_weight_list[i]*power(delta_1D[i],2)
            denominator+=story_weight_list[i]*delta_1D[i]
        T1_accurate=2*pi*square_root((numerator)/(9.806*denominator))
        return round(T1_accurate,3)
    @staticmethod
    def Fundamental_effective_damping_of_dampers(W_1,damper_displacement_list,maximum_damper_force):
        """Calculates the eefctive damping of dampers for the fundamental mode (beta_v1)

        Args:
            W_1 (double): Modal strain energy for the fundamental mode
            damper_displacement_list (list): List containing the absolute damper displacement values for the fundamental mode
            maximum_damper_force (double): Maximum nominal yield force of dampers (provided by the device producer)
        """        
        numerator=0
        for i in range (len(damper_displacement_list)):
                numerator+= 4*damper_displacement_list[i]*maximum_damper_force
        beta_v1=numerator/(4*pi*W_1)
        return round(beta_v1,3)
    @staticmethod
    def Fundamental_damper_displacement_list(delta_1):
        """Calculates the absolute displacement for each damper at each story fo the fundamental mode

        Args:
            delta_1 (list): List containing the design deflection value for each story for the fundamental mode
        """        
        damp_disp=[0]*len(delta_1)
        for i in range(len(delta_1)):
            if i==0:
                damp_disp[i]=round(abs(delta_1[i]),3)
            else:
                damp_disp[i]=round(abs(delta_1[i]-delta_1[i-1]),3)
        return damp_disp
    @staticmethod
    def Fundamental_modal_strain_energy(F1,delta_1):
        """Calculates the modal strain energy for the fundamental mode (W_1)

        Args:
            F1 (list): List containing lateral force values for each story for the fundamental mode
            delta_1 (list): List containing the design deflection value for each story for the fundamental mode
        """        
        W_1=0
        for i in range(len(F1)):
            W_1+=0.5*(F1[i]*delta_1[i])
        return round(W_1,3)
    @staticmethod
    def Fundamental_lateral_force_lsit(w1,phi_1,Gamma_1,W_bar_1,V1):
        """Calculates the list containing lateral force values for each story for the fundamental mode (F1)

        Args:
            w1 (list): List containing story weight values
            phi_1 (list): Modal shape vector for the fundamental mode
            Gamma_1 (double): Mode participationfactro for the fundamental mode
            W_bar_1 (double): Modal mass of the structure for the fundamental mode
            V1 (double): Base shear for the fundamental mode
        """        
        F1=[0]*len(phi_1)
        for i in range(len(phi_1)):
            F1[i]=round(((w1[i]*phi_1[i]*Gamma_1*V1)/W_bar_1),3)
        return F1
    @staticmethod
    def Fundamental_design_deflection(D1D,phi_1):
        """Calculates the design deflection of each story at its
           center of rigidity for the fisrt mode (delta_1D)

        Args:
            D1D (double): Design displacement of the center of rigidity at the last roof level for the fundamental mode
            phi_1 (list): Modal shape vector of the structure for the fisrt mode
        """        
        delta_1D=[0]*len(phi_1)
        for i in range(len(phi_1)):
            delta_1D[i]=round((D1D*phi_1[i]),3)
        return delta_1D
    @staticmethod
    def Fundamental_design_roof_displacement(Gamma_1,SDS,T1D,T1,B1D,B1E,SD1):
        """Calculates the design displacement of the center of rigidity at the last roof level for the fundamental mode (D1D)

        Args:
            Gamma_1 (double): Participation factor for the fundamental mode
            SDS (double): Design, 5% damped, spectral response acceleration parameter for short periods
            T1D (double): Damped period for the fundamental mode
            T1 (double): Approximate period for the fundamental mode
            B1D (double): B value for damped performance assumption for the fundamental mode
            B1E (double): B value for the effective damping (beta_E)
            SD1 (double): Design, 5% damped, spectral response acceleration parameter at a period of 1 s
        """   
        Ts=SD1/SDS     
        if T1D<Ts:
            a = (9.81*Gamma_1*SDS*power(T1D,2))/(4*B1D*power(pi,2))
            b = (9.81*Gamma_1*SD1*power(T1,2))/(4*B1E*power(pi,2))
            if b<=a:
                D1D=a
            else:
                D1D=b
        else:
            a = (9.81*Gamma_1*SD1*T1D)/(4*B1D*power(pi,2))
            b = (9.81*Gamma_1*SD1*T1)/(4*B1E*power(pi,2))
            if a<=b:
                D1D=a
            else:
                D1D=b
        return round(D1D,3)
    @staticmethod
    def Fundamental_design_base_shear(Cs1,W_bar_1):
        """Calculates the base shear design value for the fundamental mode (V1)

        Args:
            Cs1 (double): Seismic response coefficient for the fundamental mode
            W_bar_1 (double): Total modal mass of the structure for the fundamental mode
        """        
        V1=Cs1*W_bar_1
        return round(V1,3)
    @staticmethod
    def Fundamental_response_coefficient(Ts,T1D,R,Cd,Omega_0,SD1,B1D):
        """Calculates the seismic response coefficient for the fundamental mode (Cs1)

        Args:
            Ts (double): Ts ratio (= sd1/sds)
            T1D (double): Damped period of the structure for the fundamental mode
            R (double): Response modification coefficient
            Cd (double): Deflection amplification factor
            Omega_0 (double): Over strength factor
            SD1 (double): Design, 5% damped, spectral response acceleration parameter at a period of 1 s 
            B1D (double): B coefficient for the damped fundamental mode assumption (Table 18-7.2, ASCE7-16)
        """        
        if T1D<Ts:
            Cs1=(R*SD1)/(Cd*Omega_0*B1D)
        else:
            Cs1=(R*SD1)/(Cd*T1D*Omega_0*B1D)
        return round(Cs1,3)
    @staticmethod
    def Fundamental_damped_period(approxiamte_period,mu_D):
        """Calculates the damped period of the structure for the fundamental mode (T1_D)

        Args:
            approxiamte_period (double): Approximate value of the period of fundamental mode for the given structure
            mu_D (double): Effective ductility demand on the seismic force-resisting system
        """        
        T1_D=approxiamte_period*square_root(mu_D)
        return round(T1_D,3)
    @staticmethod
    def Fundamental_approxiamte_period(structure_type,total_structure_height):
        """"Calculates the approximate value of the period of fundamental mode for the given structure
            based on Table 12.8-2 of ASCE7-16 (T1)

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
    def Fundamental_modal_shape_vector(story_height_list):
        """Calculates the modal shape vector of the fundamental mode (Phi_i1) for the given structure according to
           procedures described in chapter 18, ASCE7-16 (phi_i1)

        Args:
            story_height_list (list): An array containing the story height values from the base floor to the top roof, respectively
        """       
        phi_i1=[0]*(len(story_height_list)-1) 
        total_height = story_height_list[len(story_height_list)-1]-story_height_list[0]
        story_height=0
        for i in range (len(story_height_list)-1):
            story_height+=story_height_list[i+1]-story_height_list[i]
            phi_i1[i]=round((story_height/total_height),3)
        return phi_i1
    @staticmethod
    def Fundamental_modal_mass(story_weight_list,modal_shape_vector_1):
        """Calculates the total modal mass of the structure for the fundamental mode (W_bar_1)

        Args:
            story_weight_list (list): List containing weight values for each story from the first floor to the top roof, respectively
            modal_shape_vector_1 (list): The modal shape vector of the fundamental mode for the given structure from the first floor to the top roof, respectively
        """ 
        numerator=0
        denominator=0
        for i in range(len(story_weight_list)):
            numerator+=story_weight_list[i]*modal_shape_vector_1[i]
        numerator=power(numerator,2)
        for i in range (len(story_weight_list)):
            denominator+=power(modal_shape_vector_1[i],2)*story_weight_list[i]
        total_modal_mass=numerator/denominator
        return total_modal_mass
    @staticmethod
    def Gamma_1(total_modal_mass,story_weight_list,modal_shape_vector_1):
        """Calculates the mode participation factor for the fundamental mode (Gamma_1)

        Args:
            total_modal_mass (double):The total modal mass of the structure
            story_weight_list (list):List containing weight values for each story from the first floor to the top roof, respectively
            modal_shape_vector_1 (list): The modal shape vector of the fundamental mode for the given structure from the first floor to the top roof, respectively
        """        
        denominator=0
        for i in range(len(story_weight_list)):
            denominator+=modal_shape_vector_1[i]*story_weight_list[i]
        Gamma_1=total_modal_mass/denominator
        return round(Gamma_1,3)
    @staticmethod
    def Fundamental_total_damping(beta_I,beta_v1,mu_D,beta_HD):
        """Calculates to total effective damping value for the fundamental mode (beta_1D)

        Args:
            beta_I (double): Inherent damping (<=3%)
            beta_v1 (double): Damping of dampers (nearly any kind)
            mu_D (double): Effective ductility demand on the seismic force-resisting system
            beta_HD (double): Hysteretic damping
        """   
        if beta_I>0.03 or beta_I<0:
            print("**Error! Inherent damping value must be in range 0<=beta_I<=0.03! revise your input!")
            return 0
        elif not 1<=mu_D:
            print("**Error! Ductility demand value must be equal to or greater than 1! revise your input!")
            return 0
        else:
            beta_1D=beta_I+(beta_v1*square_root(mu_D))+beta_HD
            return round(beta_1D,3)
    @staticmethod
    def Effective_damping(beta_I,beta_v1):
        """Calculates the effective damping value (beta_E)

        Args:
            beta_I (double): Inherent damping value (<=3%)
            beta_v1 (double): Damping of dampers (nearly any kind)
        """        
        beta_E=beta_I+beta_v1
        return round(beta_E,3)
class Residual_Mode_Parameters: # A class for the residual-mode-related calculations
    @staticmethod
    def Residual_effective_damping_of_dampers(W_R,damper_displacement_list,maximum_damper_force):
        """Calculates the eefctive damping of dampers for the residual mode (beta_vr)

        Args:
            W_R (double): Modal strain energy for the residual mode
            damper_displacement_list (list): List containing the absolute damper displacement values for the residual mode
            maximum_damper_force (double): Maximum nominal yield force of dampers (provided by the device producer)
        """
        if W_R==0:
            beta_v1=0
        else:        
            numerator=0
            for i in range (len(damper_displacement_list)):
                    numerator += (4*damper_displacement_list[i]*maximum_damper_force)
            beta_v1=numerator/(4*pi*W_R)
        return round(beta_v1,3)
    @staticmethod
    def Residual_damper_displacement_list(delta_R):
        """Calculates the absolute displacement for each damper at each story fo the residual mode

        Args:
            delta_R (list): List containing the design deflection value for each story for the residual mode
        """        
        damp_disp=[0]*len(delta_R)
        for i in range(len(delta_R)):
            if i==0:
                damp_disp[i]=round(abs(delta_R[i]),3)
            else:
                damp_disp[i]=round(abs(delta_R[i]-delta_R[i-1]),3)
        return damp_disp
    @staticmethod
    def Residual_modal_strain_energy(FR,delta_R):
        """Calculates the modal strain energy for the residual mode (W_R)

        Args:
            F1 (list): List containing lateral force values for each story for the residual mode
            delta_1 (list): List containing the design deflection value for each story for the residual mode
        """        
        W_R=0
        for i in range(len(FR)):
            W_R+=0.5*(FR[i]*delta_R[i])
        return round(W_R,3)
    @staticmethod
    def Residual_lateral_force_lsit(w1,phi_R,Gamma_R,W_bar_R,VR):
        """Calculates the list containing lateral force values for each story for the residual mode(FR)

        Args:
            w1 (list): List containing story weight values
            phi_R (list): Modal shape vector for the residual mode
            Gamma_R (double): Mode participationfactro for the residual mode
            W_bar_R (double): Modal mass of the structure for the residual mode
            VR (double): Base shear for the residual mode
        """        
        FR=[0]*len(phi_R)
        for i in range(len(phi_R)):
            FR[i]=round(((w1[i]*phi_R[i]*Gamma_R*VR)/W_bar_R),3)
        return FR
    @staticmethod
    def Residual_design_deflection(DRD,phi_R):
        """Calculates the design deflection of each story at its
           center of rigidity for the residual mode (delta_RD)

        Args:
            DRD (double): Design displacement of the center of rigidity at the last roof level for the residual mode
            phi_R (list): Modal shape vector of the structure for the residual mode
        """        
        delta_RD=[0]*len(phi_R)
        for i in range(len(phi_R)):
            delta_RD[i]=round((DRD*phi_R[i]),3)
        return delta_RD
    @staticmethod
    def Residual_design_roof_displacement(Gamma_R,SDS,TR,BR,SD1):
        """Calculates the design displacement of the center of rigidity at the last roof level for the residual mode (DRD)

        Args:
            Gamma_R (double): Participation factor for the residual mode
            SDS (double): Design, 5% damped, spectral response acceleration parameter for short periods
            Tr (double): Period for the residual mode
            BR (double): B value for the residual mode
            SD1 (double): Design, 5% damped, spectral response acceleration parameter at a period of 1 s
        """        
        a = (9.81*Gamma_R*SD1*TR)/(4*BR*power(pi,2))
        b = (9.81*Gamma_R*SDS*power(TR,2))/(4*BR*power(pi,2))
        if a<=b:
            DRD=b
        else:
            DRD=a
        return round(DRD,3)
    @staticmethod
    def Residual_design_base_shear(Csr,W_bar_R):
        """Calculates the base shear design value for the residual mode (Vr)

        Args:
            Csr (double): Seismic response coefficient for the residual mode
            W_bar_R (double): Total modal mass of the structure for the residual mode
        """        
        Vr=Csr*W_bar_R
        return round(Vr,3)
    @staticmethod
    def Residual_mode_response_coefficient(R,Cd,Omega_0,SDS,BR):
        """Calculates the seismic response coefficient for the residual mode (CsR)

        Args:
            R (double): Response modification coefficient
            Cd (double): Deflection amplification factor
            Omega_0 (double): Over strength factor
            SDS (double): Design, 5% damped, spectral response acceleration parameter for short periods
            BR (double): B coefficient for the redidual mode (Table 18-7.2, ASCE7-16)
        """        
        CsR=(R*SDS)/(Cd*Omega_0*BR)
        return round(CsR,3)
    @staticmethod
    def Residual_mode_period(fundamental_mode_period):
        """Calculates the period for the residual mode of a given structure

        Args:
            fundamental_mode_period (double): Period of the fundamental mode of the structure
        """        
        Tr=0.4*fundamental_mode_period
        return round(Tr,3)
    @staticmethod
    def Gamma_R(gamma_1):
        """Calculates the mode participation factor (Gamma_R) for the residual mode

        Args:
            gamma_1 (double): The mode participation factor for the fist mode
        """       
        gamma_r=1-gamma_1
        return round(gamma_r,3)
    @staticmethod
    def Residual_modal_shape_vector(modal_shape_vector_1,gamma_1):
        """Calculates the modal shape vector of the residual mode (Phi_ir) for the given structure according to
           procedures described in chapter 18, ASCE7-16

        Args:
            modal_shape_vector_1 (list): The modal shape vector of the fundamental mode for the given structure from the first floor to the top roof, respectively
            gamma_1 (double): The mode participation factor for the fist mode
        """       
        phi_ir=[0]*(len(modal_shape_vector_1)) 
        for i in range (len(modal_shape_vector_1)):
            temp=(1-gamma_1*modal_shape_vector_1[i])/(1-gamma_1)
            phi_ir[i]=round(temp,3)
        return phi_ir
    @staticmethod
    def Residual_modal_mass(total_mass,modal_mass_1):
        """Calculates the total modal mass of the structure for the residual mode

        Args:
            total_mass (double): The total mass of the structure
            modal_mass_1 (double): The total modal mass of the structure for the fundamental mode
        """        
        modal_mass_r=total_mass-modal_mass_1
        return modal_mass_r
    @staticmethod
    def Residual_mode_total_damping(beta_I,beta_vr,beta_HD):
        """Calculates to total effective damping value for the residual mode (beta_R)

        Args:
            beta_I (double): Inherent damping e (<=3%)
            beta_vr (double): Damping of dampers (nearly any kind) for the residual mode
            beta_HD (double): Hysteretic damping
        """   
        if beta_I>0.03 or beta_I<0:
            print("**Error! Inherent damping value must be in range 0<=beta_I<=0.03! revise your input!")
            return 0
        else:
            beta_r=beta_I+beta_vr+beta_HD
            return round(beta_r,3)
class Combinatory_Calculations: # A class for conducting calculations related to both the fundamental, and the residual modes
    @staticmethod
    def Effective_yield_displacement(Omega_0,Cd,R,Gamma_1,Cs1,T1):
        """Calculates the effective yield displacement of the center of rigidity of the last roof (Dy)

        Args:
            Omega_0 (double): Over strength factor
            Cd (double): Deflection amplification factor
            R (double): Response modification factor
            Gamma_1 (double): Participation factor for the fundamental mode
            Cs1 (double): Seismic response coefficient for the fundamental mode
            T1 (double): Approximate period of the structure for the fundamental mode
        """        
        Dy=(9.81*Omega_0*Cd*Gamma_1*Cs1*power(T1,2))/(4*R*power(pi,2))
        return round(Dy,6)
    @staticmethod
    def Existing_ductility_demand(Dy,D1D):
        """Calculates the existing ductility demand of the structure (mu_D)

        Args:
            Dy (double): Effective yield displacement of the center of rigidity of the last roof
            D1D (double): Design displacement of the center of rigidity at the last roof level for the fundamental mode
        """        
        mu_D=D1D/Dy
        return round(mu_D,3)
    @staticmethod
    def Maximum_ductility_demand(Ts,T1D,R,Omega_0,Ie):
        """Calculates the maximum achieveable ductility demand under the conditions of the strutcure (mu_max)

        Args:
            Ts (double): Ts ratio (=SD1/SDS)
            T1D (double): Damped period for the fundamental mode
            R (double): Response modification factor
            Omega_0 (double): Over strength factor
            Ie (double): Importance factor
        """        
        if T1D<=Ts:
            mu_max=0.5*(power((R/(Omega_0*Ie)),2)+1)
        else:
            mu_max=R/(Omega_0*Ie)
        return round(mu_max,3)
    @staticmethod
    def Combined_base_shear(V1,VR):
        """Calculates the combined base shear using SRSS method (V)

        Args:
            V1 (double): Base shear for the fundamental mode
            VR (double): Base shear for the residual mode
        """        
        V=square_root(power(V1,2)+power(VR,2))
        return round(V,3)
    @staticmethod
    def Total_damper_displacement_vector(damp_disp_1,damp_disp_r):
        """Calculates the total displacement vector for dampers

        Args:
            damp_disp_1 (double): Damper displacement vector in fundamental mode
            damp_disp_r (double): Damper displacement vector in residual mode
        """        
        total_damper_displacement=[0]*len(damp_disp_1)
        for i in range(len(damp_disp_1)):
            total_damper_displacement[i]=round((square_root(power(damp_disp_r[i],2)+power(damp_disp_1[i],2))),3)
        return total_damper_displacement
    @staticmethod
    def Fixed_base_shear(R,SDS,SD1,Ie,T1,W):
        """Calculates the fixed base shearing force

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
    def Minimum_base_shear(V_Fixed):
        """Calculates the minimum base shearing force

        Args:
            V_Fixed (double): Fixed base shearing force
        """        
        V_min = 0.75*V_Fixed
        return round(V_min,3)
