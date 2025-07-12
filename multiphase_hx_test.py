import thermopt as th
import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt
import ht

def interface_capturing(fluid,h_in,p_in,h_out,p_out,invert=False):
    """
    Approximate the phase change points using a linear relationship between
    enthalpy and pressure.
    """
    #internally algorithm works assuming cooling, swap if input is heating
    if h_out>h_in:
        invert = not(invert)
        h_in,p_in,h_out,p_out=h_out,p_out,h_in,p_in
    eps=1e-5 #margin to ensure Ph point is inside dome
    p=[p_in]
    h=[h_in]
    #first guess, pahse change at intermediate pressure between inlet and outlet
    p_av=0.5*(p_in+p_out)
    vap_sat=fluid.get_state(th.PQ_INPUTS,p_av,1.0-eps)
    x_vap_sat=(vap_sat["h"]-h_in)/(h_out-h_in)
    p_av=0.5*(p_in+p_out)
    liq_sat=fluid.get_state(th.PQ_INPUTS,p_av,0.0+eps)
    x_liq_sat=(liq_sat["h"]-h_in)/(h_out-h_in)
    #if phase transition is within process, add points.
    #Pressure is estimated from previous guess, then PQ call used to compute enthalpy
    if 0<x_vap_sat<1:
        p_vap_sat=p_in+x_vap_sat*(p_out-p_in)
        vap_sat=fluid.get_state(th.PQ_INPUTS,p_vap_sat,1.0-eps)
        h.append(vap_sat["h"])
        p.append(vap_sat["p"])
    if 0<x_liq_sat<1:
        p_liq_sat=p_in+x_vap_sat*(p_out-p_in)
        liq_sat=fluid.get_state(th.PQ_INPUTS,p_vap_sat,0.0+eps)
        h.append(liq_sat["h"])
        p.append(liq_sat["p"])
    h.append(h_out)
    p.append(p_out)
    if invert:
        h=h[::-1]
        p=p[::-1]
    return h,p 

def hx_correl(i,qq,G,fld,liq,gas,geom,heating=True):
    D_h=geom["hydraulic_diameter"]
    two_phase=((0<fld[i]["Q"]<1) and (0<fld[i+1]["Q"]<1))
    if two_phase:
        A_channel_flow=geom["cross_section_area"]
        m = G*A_channel_flow
        Q=0.5*(fld[i]["Q"]+fld[i+1]["Q"])
        rho_l=0.5*(liq[i]["rhomass"]+liq[i+1]["rhomass"])
        rho_g=0.5*(gas[i]["rhomass"]+gas[i+1]["rhomass"])
        mu_l=0.5*(liq[i]["viscosity"]+liq[i+1]["viscosity"])
        mu_g=0.5*(gas[i]["viscosity"]+gas[i+1]["viscosity"])
        k_l=0.5*(liq[i]["conductivity"]+liq[i+1]["conductivity"])
        Hvap=0.5*(gas[i]["hmass"]+gas[i+1]["hmass"]-liq[i]["hmass"]-liq[i+1]["hmass"])
        cp_l=0.5*(liq[i]["cpmass"]+liq[i+1]["cpmass"])
        if heating:
            ht_coeff=ht.boiling_plate.h_boiling_Han_Lee_Kim(m, Q, D_h, rho_l, rho_g, mu_l, k_l, \
            Hvap, cp_l, qq, A_channel_flow, geom["wavelength"], chevron_angle=geom["chevron_angle"])
        else:
            ht_coeff=ht.condensation.Cavallini_Smith_Zecchin(m, Q, D_h, rho_l, rho_g, mu_l, mu_g, k_l, cp_l)
        #TODO: replace gas pressure loss with Lockart-Martinelli or similar correlations
        Re_g = G*D_h/mu_g
        fr = ht.conv_plate.friction_plate_Martin_1999(Re_g, geom["chevron_angle"])
        dPdx=G**2 / (2*rho_g) * fr / D_h
    else:
        rho = 0.5*(fld[i]["rhomass"]+fld[i+1]["rhomass"])
        cp = 0.5*(fld[i]["cpmass"]+fld[i+1]["cpmass"])
        mu = 0.5*(fld[i]["viscosity"]+fld[i+1]["viscosity"])
        k = 0.5*(fld[i]["conductivity"]+fld[i+1]["conductivity"])
        Re = G*D_h/mu
        Pr = cp*mu/k
        Nu = ht.conv_plate.Nu_plate_Martin(Re, Pr, geom["chevron_angle"])
        fr = ht.conv_plate.friction_plate_Martin_1999(Re, geom["chevron_angle"])
        ht_coeff=Nu*k/D_h
        dPdx = G**2 / (2*rho) * fr / D_h
    return ht_coeff,dPdx

def heat_exchanger_segmenting(fluid_hot,fluid_cold,h_hot,p_hot,h_cold,p_cold,num_steps=10):
    #merges the list of points where phase transition happens.
    #the points in x are computed as fraction of total heat transfered
    x_hot =[(hi-h_hot[0] )/(h_hot[-1] -h_hot[0] ) for hi in h_hot]
    x_cold=[(hi-h_cold[0])/(h_cold[-1]-h_cold[0]) for hi in h_cold]
    x_merged = list(set(x_hot+x_cold))
    x_merged.sort()
    #subdivide each element in steps
    x_refined=[]
    for i in range(len(x_merged)-1):
        x_refined.append(x_merged[i])
        intermediate_points=np.linspace(x_merged[i], x_merged[i+1], num_steps)[1:-1]
        x_refined.extend(intermediate_points)
    x_refined.append(x_merged[-1])

    h_hot_ref = np.interp(x_refined, x_hot, h_hot)
    p_hot_ref = np.interp(x_refined, x_hot, p_hot)
    h_cold_ref = np.interp(x_refined, x_cold, h_cold)
    p_cold_ref = np.interp(x_refined, x_cold, p_cold)
    #compute properties
    hot  = []; hot_l  = []; hot_g  = []
    cold = []; cold_l = []; cold_g = []
    for i in range(len(h_hot_ref)):
        f=fluid_hot
        hot.append(f.get_state(th.HmassP_INPUTS,h_hot_ref[i],p_hot_ref[i]))
        if 0<hot[-1]["Q"]<1:
            hot_l.append(f.get_state(th.PQ_INPUTS,p_hot_ref[i],0.0))
            hot_g.append(f.get_state(th.PQ_INPUTS,p_hot_ref[i],1.0))
        else:
            hot_l.append(None)
            hot_g.append(None)
        f=fluid_cold
        cold.append(f.get_state(th.HmassP_INPUTS,h_cold_ref[i],p_cold_ref[i]))
        if 0<cold[-1]["Q"]<1:
            cold_l.append(f.get_state(th.PQ_INPUTS,p_cold_ref[i],0.0))
            cold_g.append(f.get_state(th.PQ_INPUTS,p_cold_ref[i],1.0))
        else:
            cold_l.append(None)
            cold_g.append(None)

    return x_refined,hot,hot_l,hot_g,cold,cold_l,cold_g

def hx_segment(i,hot,hot_l,hot_g,cold,cold_l,cold_g,G_hot,G_cold,geom):
    DT_DEGENERATE=1e-3 #avoid errors with logarithm in degenerate cases
    #check if two phase correlation
    tho_phase_hot  = ((0<hot[i]["Q"]<1) and (0<hot[i+1]["Q"]<1))
    two_phase_cold = ((0<cold[i]["Q"]<1) and (0<cold[i+1]["Q"]<1))
    two_phase = tho_phase_hot or two_phase_cold
    dT1=hot[i]["T"]-cold[i]["T"]
    dT2=hot[i+1]["T"]-cold[i+1]["T"]
    #LMTD
    if dT1<DT_DEGENERATE or dT2<DT_DEGENERATE: 
        LMTD=DT_DEGENERATE
    elif np.isclose(dT1,dT2):
        LMTD=(dT1+dT2)*0.5
    else:
        LMTD=(dT1-dT2)/np.log(dT1/dT2)
    #iterative solution if two phase, otherwise direct solution
    if two_phase: 
        fun=lambda qq: qq-LMTD/(1/hx_correl(i,qq,G_cold,cold,cold_l,cold_g,geom,heating=True)[0] \
                               +1/hx_correl(i,qq,G_hot ,hot ,hot_l ,hot_g ,geom,heating=False)[0] \
                               +geom["wall_resistance"])
        try: 
            q = brentq(fun,1e2,1e6)
        except:
            if fun(1e2)>0:
                q = 1e2
            else:
                q = 1e6
    else:
        qq=np.nan #not used for single phase
        q=                LMTD/(1/hx_correl(i,qq,G_cold,cold,cold_l,cold_g,geom,heating=True)[0] \
                               +1/hx_correl(i,qq,G_hot ,hot ,hot_l ,hot_g ,geom,heating=False)[0] \
                               +geom["wall_resistance"])
    Area=m_dot_cold*(cold[i+1]["h"]-cold[i]["h"])/q
    Length=Area/geom["cross_section_area"]*geom["hydraulic_diameter"]/4 #by defintion of hydraulic diameter
    dP_cold=Length*hx_correl(i,q,G_cold,cold,cold_l,cold_g,geom,heating=True)[1]
    dP_hot =Length*hx_correl(i,q,G_hot ,hot ,hot_l ,hot_g ,geom,heating=False)[1]
    U=q/LMTD
    return Area,dP_cold,dP_hot,LMTD,U

def composite_heat_exchanger(
    fluid_hot,
    h_in_hot,
    h_out_hot,
    p_in_hot,
    p_out_hot,
    fluid_cold,
    h_in_cold,
    h_out_cold,
    p_in_cold,
    p_out_cold,
    m_dot_hot=0.1,
    m_dot_cold=0.1,
    num_steps=10,
    counter_current=True,
    geom=    {"topology":"plate",
              "cross_section_area":0.001,
              "hydraulic_diameter":0.003,
              "wavelength":0.01,
              "wall_resistance":0.0001,#K/(W/m2)
              "chevron_angle":45}
):
    #compute a list of enthalpy and pressure which split the hx according to the phase of either fluid
    h_hot ,p_hot  = interface_capturing(fluid_hot,h_in_hot ,p_in_hot ,h_out_hot ,p_out_hot ,invert=counter_current)
    h_cold,p_cold = interface_capturing(fluid_cold,h_in_cold,p_in_cold,h_out_cold,p_out_cold)
    x,hot,hot_l,hot_g,cold,cold_l,cold_g = heat_exchanger_segmenting(fluid_hot,fluid_cold,h_hot,p_hot,h_cold,p_cold,num_steps=num_steps)
    G_hot = m_dot_hot/geom["cross_section_area"]
    G_cold = m_dot_cold/geom["cross_section_area"]
    Area=0;dP_cold_tot=0;dP_hot_tot=0
    for i in range(len(hot)-1):
        dArea,dP_cold,dP_hot,LMTD,U=hx_segment(i,hot,hot_l,hot_g,cold,cold_l,cold_g,G_hot,G_cold,geom)
        Area+=dArea
        dP_cold_tot+=dP_cold
        dP_hot_tot+=dP_hot
        #TODO: for debug, remove and implement as residuals
        plt.figure(1)
        plt.plot(x[i],hot[i]["T"],'r*')
        plt.plot(x[i],cold[i]["T"],'b*')
        plt.figure(2)
        plt.plot(x[i:i+2],[U,U],'k')
        plt.figure(3)
        plt.plot(x[i],1-dP_hot_tot/(p_in_cold-p_out_cold),'*r')
        plt.plot(x[i],dP_cold_tot/(p_in_cold-p_out_cold),'*b')
    print("HX area:",Area)
    print("Cold pressure loss factor",dP_cold_tot/(p_in_cold-p_out_cold))
    print("Hot pressure loss factor",dP_hot_tot/(p_in_hot-p_out_hot))
    print("Heat balance factor",m_dot_hot*(h_in_hot-h_out_hot)/(m_dot_cold*(h_out_cold-h_in_cold)))


if __name__=="__main__":
    fluid_hot=th.Fluid("R245fa",backend="HEOS")
    p_in_hot=2.2e6
    p_out_hot=2.15e6
    h_in_hot=550e3
    h_out_hot=330e3
    m_dot_hot=2.0
    #fluid_cold=th.Fluid("R245fa",backend="HEOS")
    #p_in_cold=1.5e6
    #p_out_cold=1.45e6
    #h_in_cold=280e3
    #h_out_cold=530e3
    #m_dot_cold=2.0*0.88
    fluid_cold=th.Fluid("Water",backend="HEOS")
    p_in_cold=1e7
    p_out_cold=0.98e7
    h_in_cold=280e3
    h_out_cold=320e3
    m_dot_cold=2.0*99/18

    geom=  {"topology":"plate",
            "cross_section_area":0.008,
            "hydraulic_diameter":0.003,
            "wavelength":0.01,
            "wall_resistance":0.0003,#K/(W/m2)
            "chevron_angle":45}

    composite_heat_exchanger(
    fluid_hot,
    h_in_hot,
    h_out_hot,
    p_in_hot,
    p_out_hot,
    fluid_cold,
    h_in_cold,
    h_out_cold,
    p_in_cold,
    p_out_cold,
    m_dot_hot=m_dot_hot,
    m_dot_cold=m_dot_cold,
    num_steps=5,
    counter_current=True,
    geom=geom)
    plt.show()
