import logging
progress = {}
logger = logging.Logger('Chemviewer')
formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
ch = logging.StreamHandler()
ch.setLevel('DEBUG')
ch.setFormatter(formatter)
logger.addHandler(ch)
