import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



def get_model_nr(df):
    '''Gets the model number'''
    return df.iloc[:,0].to_numpy()

def get_age(df):
    '''Gets the age of the star'''
    return df.iloc[:,1].to_numpy()
    
def get_mass(df):
    '''Gets the total mass of the star'''
    return df.iloc[:,5].to_numpy()

def get_he_core_mass(df):
    '''Gets the hydrogen exhausted core mass (X<0.1 per STARS default definition)'''
    return df.iloc[:,6].to_numpy()

def get_co_core_mass(df):
    '''Gets the carbon-oxygen core mass'''
    return df.iloc[:,7].to_numpy()


def get_m_conv(df):
    '''Gets the mass coordinates of convective-radiative boundaries'''
    return df.iloc[:,11:22].to_numpy()

def get_conv_env_mass(df):
    '''Gets the mass of the convective envelope'''
    return df.iloc[:,70].to_numpy()


def read_plot_file(filename):
        
    spec = ([6,16]                               # I6, E16.9
    + [10 for _ in range(24)]               # 24F10.5
    + [13, 13, 13]                          # 3E13.6
    # should the i in the following LC be a _?
    + [12 for _ in range(18)]# for i in (1,12)]#18(1X,E12.5)
    + [9 for _ in range(52)])               # 52F9.5

    df = pd.read_fwf(filename,widths=spec, header=None)
    return df


class Zone:
    def __init__(self, x, y_min, y_max):
        self.x = [x]
        self.y_min = [y_min]
        self.y_max = [y_max]
        self.mid = (y_max + y_min)/2
        
    def append(self, x, y_min, y_max):
        self.x.append(x)
        self.y_min.append(y_min)
        self.y_max.append(y_max)
        self.mid = (y_max + y_min)/2
        
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('name')
    parser.add_argument('filename')
    args = parser.parse_args()
    
    
    df = read_plot_file(args.filename)
    conv_reg = get_m_conv(df)
    top_mass = get_conv_env_mass(df)
    model_numbers = get_model_nr(df)
    total_mass = get_mass(df)
    zones = []
    conv_envelope_tolerance = 0.001

    for inter in range(1, len(conv_reg)):

        #boundaries = conv_reg.iloc[i,:][conv_reg.iloc[i,:] > 0].to_numpy()
        mass = total_mass[inter]
        N = model_numbers[inter]
        # skip if the next model number is the same
        if inter < len(model_numbers)-1:
            if N == model_numbers[inter+1]:
                continue
            
        sorted_by_mass = np.array(sorted(conv_reg[inter,:][conv_reg[inter,:] != 0], key=lambda x: abs(x)))
        conv = False
        semi_conv = False
        czone_start = 0
        czone_end = 0
        sczone_start = 0
        sczone_end = 0
        
        sorted_by_mass = sorted_by_mass[np.abs(sorted_by_mass) < mass-0.0001]

        if len(sorted_by_mass) > 0:
            if sorted_by_mass[0] < 0:
                sorted_by_mass = np.concatenate([[1e-5], sorted_by_mass])

            # odd number of boundaries
            if np.sum(sorted_by_mass > 0) % 2:
                # Check where negative values are.
                tmp = np.where(sorted_by_mass < 0)[0]
                # take a zone before the first negative value
                new_sorted_by_mass = sorted_by_mass.copy()  
                if len(tmp) != 0 and tmp[0] == 0:
                    sorted_by_mass = np.concatenate([[1e-5], sorted_by_mass])
                
                elif len(tmp) != 0 and tmp[0] == 2 and sorted_by_mass[0] > 0 and sorted_by_mass[1] > 0:
                    sorted_by_mass = np.concatenate([[1e-5], sorted_by_mass])
                    


        for i in sorted_by_mass:
            if i > 0 and conv == False:
                conv = True
                czone_start = i
            elif i > 0 and conv == True:
                conv = False
                czone_end = i
                if not np.isclose(czone_start - czone_end, 0 ):
                    
                    # find closest zone
                    exists = False
                    for z in zones:
                        if np.isclose(z.y_min[-1], czone_start, atol=0.5) and np.isclose(z.y_max[-1], czone_end, atol=0.5) and (N - z.x[-1] < 100):
                            z.append(N, czone_start, czone_end)
                            exists = True
                            break
                    if exists == False:
                        z = Zone(N, czone_start, czone_end)
                        zones.append(z)
                
            elif i < 0 and conv == True and semi_conv == False:
                semi_conv = True
                sczone_start = np.abs(i)
            elif i < 0 and conv == True and semi_conv == True:
                semi_conv = False
                sczone_end = np.abs(i)
            else:
                print('unexpected state')
                print(i, conv, semi_conv)
                print(sorted_by_mass)
                print(mass)
                print(conv_reg[inter,:])
                break
        if conv:
            exists = False
            if top_mass[inter] > conv_envelope_tolerance:
                czone_end = mass
            
                for z in zones:
                    if np.isclose(z.y_min[-1], czone_start, atol=0.5) and np.isclose(z.y_max[-1], czone_end, atol=0.5) and (N - z.x[-1] < 100):
                        z.append(N, czone_start, czone_end)
                        exists = True
                        break
                if exists == False:
                    z = Zone(N, czone_start, czone_end)
                    zones.append(z)
        
      
    
    plt.figure(figsize=(7,5))
    for z in zones:
        plt.fill_between(z.x, z.y_min, z.y_max, color='#d2f8d2')
    
    # plot the main stellar properties
    plt.plot(model_numbers, total_mass, color='black', lw=2)
    
    # bSTARS will return a zero value for the helium core mass
    # if there's no envelope left.
    he_core_mass = get_he_core_mass(df)
    
    # check where the helium core mass is zero
    mask = np.where(he_core_mass > 0)[0]
    if len(mask) > 0:
        mask2 = np.where(he_core_mass[mask[-1]:] == 0)[0]
        if len(mask2) > 0:
            he_core_mass[mask[-1]+1:] = total_mass[mask[-1]+1:]
        
    
    # plot the helium core mass
    plt.plot(model_numbers, get_he_core_mass(df), color='#1f77b4', lw=2)
    # plot the carbon-oxygen core mass
    plt.plot(model_numbers, get_co_core_mass(df),  color='#b40424', lw=2)
    
    plt.ylim(0, np.max(total_mass)+2)
    plt.xlim(model_numbers[0], model_numbers[-1])
    plt.savefig(f'{args.name}.png')
    
        
            
        
        # zones = calculate_zones(df)
        # star = get_cores(df)
        # N = get_index(df)
        # plt.figure(figsize=(15,15))

        # for z in zones:
        #     plt.fill_between(z.x, z.y_min, z.y_max, color='#d2f8d2')

        # plt.plot(N, star[0], lw=5, color='black')
        # if (star[1] != 0).any():
        #     mask = np.where(star[1] != 0)[0][-1]
        #     star[1][mask:] = star[0][mask:]
        # plt.plot(N, star[1], lw=5, color='#1f77b4')
        # plt.plot(N, star[2], lw=5, color='#b40424')
        # plt.ylim(0,star[0][0]+2)
        # plt.xlim(N[0], N[-1])
        # plt.savefig(f'{args.name}.png')