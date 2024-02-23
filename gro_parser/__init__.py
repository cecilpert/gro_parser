from . import system
from .base_logger import logger
import logging

def load_system(gro_file:str, top_file:str = None, itp_file: str = None):
    return system.GroSystem(gro_file, top_file, itp_file)


def set_logger_info():
    logger.setLevel(logging.INFO)
    
def set_logger_debug():
    logger.setLevel(logging.DEBUG)


