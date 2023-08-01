from . import gro

def load_system(gro_file:str, top_file:str):
    system = gro.GroSystem(gro_file, top_file)
    return system