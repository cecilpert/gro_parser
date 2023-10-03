from . import gro
from .base_logger import logger
import logging

def load_system(gro_file:str, top_file:str):
    system = gro.GroSystem(gro_file, top_file)
    return system

def set_logger_info():
    logger.setLevel(logging.INFO)
    
def set_logger_debug():
    logger.setLevel(logging.DEBUG)


