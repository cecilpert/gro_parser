from . import system
from .base_logger import logger
import logging

def load_system(gro_file:str, top_file:str = None, itp_file: str = None, name = None):
    logger.info(f'gro : {gro_file}')
    if itp_file:
        logger.info(f'itp : {itp_file}')
    return system.GroSystem(gro_file, top_file, itp_file, name)


def set_logger_info():
    logger.setLevel(logging.INFO)
    
def set_logger_debug():
    logger.setLevel(logging.DEBUG)


