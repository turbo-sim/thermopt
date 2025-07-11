import CoolProp as cp

fluid_cold=cp.AbstractState("REFPROP","R1233zd(E)")
fluid_cold.update(cp.HmassP_INPUTS,393327.87288754666,1477334.4254224906)
print("Hmass HP inputs (in):", 93327.87288754666)
print("Hmass HP inputs (out):", fluid_cold.hmass())
fluid_cold.update(cp.PQ_INPUTS,1477334.4254224906,0.08152285)
print("Hmass PQ inputs:", fluid_cold.hmass())


fluid_cold=cp.AbstractState("REFPROP","R1233zd(E)")
print("Phase envelope data:", fluid_cold.get_phase_envelope_data().T)
