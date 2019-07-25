
import logging
logger = logging.getLogger(__name__)

logger.warn("This module is deprecated and just a wrapper for netgen.occ, import netgen.occ instead")

from .occ import *
