import numpy as np
import matplotlib.pyplot as plt
from skplatform.geodesy import Geodetic

def method_one(va, vb):
    nva = va / np.linalg.norm(va)
    nvb = vb / np.linalg.norm(vb)
    angle = np.arccos(np.dot(nva, nvb))
    cross = np.cross(nva, nvb)
    vertical = np.array([0.0, 0.0, 1.0])
    return angle*np.sign(np.dot(vertical, cross))

def method_two(va, vb, vertical):
    nva = va / np.linalg.norm(va)
    nvb = vb / np.linalg.norm(vb)
    cross = np.cross(nva, nvb)
    # vertical = np.array([0.0, 0.0, 1.0])
    return np.arctan2(np.dot(cross, vertical), np.dot(nva, nvb))

def comapare_angles():
    N = 100_000
    angle = np.zeros(N)
    a1 = np.zeros(N)
    a2 = np.zeros(N)
    for ai , a in enumerate(np.linspace(-np.pi, np.pi, N)):
        angle[ai] = np.rad2deg(-a)
        # 2D case
        axis = np.array([1.0, 0.0, 0.0])
        rotated = np.array([np.cos(a), np.sin(a), 0.0])
        
        a1[ai] = np.rad2deg(method_one(rotated, axis))
        a2[ai] = np.rad2deg(method_two(rotated, axis))
        # ( angle, np.rad2deg(a1), np.rad2deg(a2))
    
    print(np.max(np.abs(angle - a1)), np.max(np.abs(angle - a2)))
    
    fig, ax = plt.subplots(2,1, figsize=(16.0, 9.0), dpi=150)
    ax[0].plot(angle)
    ax[0].plot(a1)
    ax[0].plot(a2)
    ax[0].legend(['Angle', 'Method 1', 'Method 2'])
    
    a1_diff = angle - a1
    a2_diff = angle - a2
    # a1_diff[a1_diff == 0.0] = np.nan
    # a2_diff[a2_diff == 0.0] = np.nan
    ax[1].plot(np.log10(a1_diff))
    ax[1].plot(np.log10(a2_diff))
    ax[1].legend(['Method 1 Error', "Method 2 Error"])
    
    # fig.show()
    fig.savefig('angles.png')
    

def calculate_angle_deg(sat_lla, tp_lla):
    sat_geo = Geodetic()
    sat_xyz = sat_geo.xyz_from_llh(sat_lla)
    # looking north and further west, expect a negative angle (CCW)
    tp_geo = Geodetic()
    tp_xyz = tp_geo.xyz_from_llh(tp_lla)
    
    north, east, down = sat_geo.xyz_north_east_down(tp_xyz)
    look = tp_xyz - sat_xyz
    return np.rad2deg(method_two(look, north, -down))
    
    
def test_angles():
    print('========== West Sat (Negative Longitude) ==========')
    west_sat = np.array([-7.2, -166.0, 490_000.0])
    nw_tp = np.array([-13.2, -169.0, 490_000.0])
    sw_tp = np.array([ 13.2, -169.0, 490_000.0])
    ne_tp = np.array([-13.2, -160.0, 490_000.0])
    se_tp = np.array([ 13.2, -160.0, 490_000.0])
    
    nw_angle = calculate_angle_deg(west_sat, nw_tp)
    print(f'NW of West Sat: {nw_angle}. Should be negative')
    sw_angle = calculate_angle_deg(west_sat, sw_tp)
    print(f'SW of West Sat: {sw_angle}. Should be negative')
    
    ne_angle = calculate_angle_deg(west_sat, ne_tp)
    print(f'NE of West Sat: {ne_angle}. Should be positive')
    se_angle = calculate_angle_deg(west_sat, se_tp)
    print(f'SE of West Sat: {se_angle}. Should be positive')

    print('')
    
    print('========== East Sat (Positive Longitude) ==========')
    east_sat = np.array([-7.2, 166.0, 490_000.0])
    ne_tp = np.array([-13.2, 169.0, 490_000.0])
    se_tp = np.array([ 13.2, 169.0, 490_000.0])
    nw_tp = np.array([-13.2, 160.0, 490_000.0])
    sw_tp = np.array([ 13.2, 160.0, 490_000.0])
    
    nw_angle = calculate_angle_deg(east_sat, nw_tp)
    print(f'NW of East Sat: {nw_angle}. Should be negative')
    sw_angle = calculate_angle_deg(east_sat, sw_tp)
    print(f'SW of East Sat: {sw_angle}. Should be negative')
    
    ne_angle = calculate_angle_deg(east_sat, ne_tp)
    print(f'NE of East Sat: {ne_angle}. Should be positive')
    se_angle = calculate_angle_deg(east_sat, se_tp)
    print(f'SE of East Sat: {se_angle}. Should be positive')
    
    
if __name__ == "__main__":
    # comapare_angles()
    test_angles()