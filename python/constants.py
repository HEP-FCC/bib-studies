"""
    List of constants
"""

C_MM_NS = 299.792   # speed of light in mm/ns
b_to_GB = 1.25e-10  # bit to GB conversion factor
MHz_to_Hz = 1e6    # MHz to Hz conversion factor
cm2_to_mm2 = 1e2    # cm^2 to mm^2 conversion factor

# Example rate assumptions, used to scale the hit densities per event to rates
IPC_BX_RATE_MHZ = 52.         # bunch crossing rate [MHz]
IPC_PIXELS_PER_CLUSTER = 5.   # pixels fired per hit (cluster size)
